#!/usr/bin/env Rscript
# ==============================================================================
# 06f_validate_synthetic_pseudobulks_MuSiC.R
# FASE 02 — AUDITORIA DE RECUPERAÇÃO DE COMPOSIÇÕES CONHECIDAS
# ==============================================================================
# HIPÓTESE
#   O MuSiC consegue recuperar misturas de composição CONHECIDA quando os
#   perfis que compõem a mistura são de doadores EXCLUÍDOS da referência?
#
# IMPORTANTE
#   - Usa somente os dados Mathys (não há bulk ROSMAP/GEO neste teste).
#   - A referência de treino exclui COMPLETAMENTE o doador que gera o teste.
#   - Perfis de teste são pseudobulks SINTÉTICOS ESPERADOS: a média por núcleo
#     de cada tipo do doador excluído é multiplicada por um número conhecido de
#     núcleos virtuais, e os perfis são somados e arredondados para counts.
#   - Isto NÃO equivale a sequenciar nova amostra de bulk RNA-seq: é um teste
#     interno (in-domain) de identificabilidade/separabilidade.
#   - Números de núcleos virtuais não significam 1.000 núcleos reais capturados.
#   - Os perfis top300 foram selecionados originalmente usando os 24 doadores:
#     isso gera vazamento de informações na seleção das FEATURES. A análise
#     primária aqui é allgenes, sem seleção supervisionada de marcadores; top300
#     é sensibilidade e NÃO constitui validação estritamente independente.
#   - Não substitui nem altera outputs de 06/06c/06d/06e.
#
# ENTRADA
#   data/reference/mathys2019/phase02_prepared/reference/
#     mathys2019_pathology_free_geneid_counts_sparse.rds
#     mathys2019_pathology_free_cell_metadata_aligned.tsv
#   data/reference/mathys2019/phase02_prepared/markers/
#     mathys2019_top300_markers.tsv
#
# SAÍDA
#   results/deconvolution/MuSiC_Synthetic_Validation/
#     06f_donor_selection.tsv
#     06f_mixture_design.tsv
#     06f_fit_audit.tsv
#     06f_proportions_long.tsv
#     06f_metrics_by_type.tsv
#     06f_metrics_by_run.tsv
#     06f_status.tsv
#     06f_sessionInfo.txt
#     plots/06f_observed_vs_expected.png
#     plots/06f_mean_absolute_error.png
#
# EXECUÇÃO
#   Rscript scripts/phase02_deconvolution/06f_validate_synthetic_pseudobulks_MuSiC.R
# ============================================================================== 

# ==============================================================================
# 1. CONFIGURAÇÕES
# ==============================================================================
CELL_TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
TEST_DONORS <- 4L    # leave-one-donor-out em quatro doadores elegíveis
SEED <- 62026L
VIRTUAL_NUCLEI_PER_MIX <- 1000L
GENE_SETS <- c("allgenes", "top300")
ZERO_TOL <- 1e-10
SUM_TOL <- 1e-6

reference_dir <- file.path("data", "reference", "mathys2019", "phase02_prepared", "reference")
markers_dir <- file.path("data", "reference", "mathys2019", "phase02_prepared", "markers")
ref_file <- file.path(reference_dir, "mathys2019_pathology_free_geneid_counts_sparse.rds")
meta_file <- file.path(reference_dir, "mathys2019_pathology_free_cell_metadata_aligned.tsv")
marker_file <- file.path(markers_dir, "mathys2019_top300_markers.tsv")
result_dir <- file.path("results", "deconvolution", "MuSiC_Synthetic_Validation")
plot_dir <- file.path(result_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

pkgs <- c("data.table", "Matrix", "MuSiC", "SingleCellExperiment", "S4Vectors")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Pacotes não instalados: ", paste(missing_pkgs, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(MuSiC)
  library(SingleCellExperiment)
  library(S4Vectors)
})

for (f in c(ref_file, meta_file, marker_file)) {
  if (!file.exists(f)) stop("Arquivo de entrada ausente: ", f)
}

# Perfis são proporções de NÚCLEOS VIRTUAIS, não frações de RNA.
compositions <- rbind(
  equal = c(Ex=.125, In=.125, Ast=.125, Oli=.125,
            Opc=.125, Mic=.125, End=.125, Per=.125),
  cortical_like = c(Ex=.38, In=.17, Ast=.12, Oli=.20,
                   Opc=.05, Mic=.04, End=.02, Per=.02),
  inhibitory_enriched = c(Ex=.22, In=.30, Ast=.12, Oli=.20,
                         Opc=.06, Mic=.05, End=.03, Per=.02),
  opc_enriched = c(Ex=.20, In=.12, Ast=.13, Oli=.22,
                   Opc=.20, Mic=.07, End=.03, Per=.03),
  vascular_enriched = c(Ex=.23, In=.12, Ast=.15, Oli=.18,
                        Opc=.06, Mic=.06, End=.10, Per=.10)
)

if (!identical(colnames(compositions), CELL_TYPES)) stop("Ordem das células nas misturas incorreta.")
if (any(compositions <= 0) || max(abs(rowSums(compositions) - 1)) > 1e-12) {
  stop("Misturas precisam ter oito frações positivas e somar 1.")
}
virtual_n <- compositions * VIRTUAL_NUCLEI_PER_MIX
if (any(abs(virtual_n - round(virtual_n)) > 1e-8)) stop("Número de núcleos virtuais não inteiro.")

mix_table <- data.table::melt(
  as.data.table(compositions, keep.rownames = "mixture"),
  id.vars = "mixture", variable.name = "cell_type", value.name = "true_fraction"
)
mix_table[, virtual_nuclei := as.integer(round(true_fraction * VIRTUAL_NUCLEI_PER_MIX))]
fwrite(mix_table, file.path(result_dir, "06f_mixture_design.tsv"), sep = "\t")

# ==============================================================================
# 2. LER E AUDITAR REFERENCE
# ==============================================================================
cat("Carregando referência Mathys...\n")
ref <- readRDS(ref_file)
if (!inherits(ref, "dgCMatrix")) ref <- as(ref, "dgCMatrix")
meta <- fread(meta_file, na.strings = c("NA", ""))
if (!all(c("TAG", "donor_id", "cell_type") %in% names(meta))) {
  stop("Metadata sem TAG/donor_id/cell_type.")
}
if (is.null(rownames(ref)) || is.null(colnames(ref)) || anyDuplicated(rownames(ref)) ||
    anyDuplicated(colnames(ref)) || anyDuplicated(meta$TAG)) {
  stop("Dimnames GeneID/TAG ausentes ou duplicadas.")
}
idx <- match(colnames(ref), as.character(meta$TAG))
if (anyNA(idx)) stop("Falha no alinhamento das colunas de núcleos.")
meta <- meta[idx]
if (!identical(as.character(meta$TAG), colnames(ref))) stop("TAG desalinhados.")
if (uniqueN(meta$donor_id) != 24L || ncol(ref) != 35140L || nrow(ref) != 17394L) {
  stop("A referência não corresponde ao conjunto main24 auditado.")
}
if (!setequal(unique(meta$cell_type), CELL_TYPES)) stop("Tipos celulares inesperados.")

marker_dt <- fread(marker_file)
if (!all(c("gene_id", "cell_type") %in% names(marker_dt))) stop("Marcadores sem GeneID/cell_type.")
marker_geneids <- unique(as.character(marker_dt$gene_id))
if (!all(marker_geneids %in% rownames(ref))) stop("Marcadores fora do universo da referência.")

# Exigir que os doadores de TESTE tenham pelo menos um núcleo observado em cada
# tipo; caso contrário, a simulação de uma mistura dos oito tipos seria inventada.
support <- meta[, .N, by = .(donor_id, cell_type)]
support_wide <- data.table::dcast(support, donor_id ~ cell_type, value.var = "N", fill = 0L)
for (ct in CELL_TYPES) if (!(ct %in% names(support_wide))) stop("Tipo não observado em nenhum doador: ", ct)
eligible <- support_wide[
  rowSums(as.matrix(support_wide[, ..CELL_TYPES]) > 0L) == length(CELL_TYPES)
]
if (nrow(eligible) < TEST_DONORS) {
  stop("Menos de ", TEST_DONORS, " doadores com oito tipos observados.")
}
set.seed(SEED)
selected_rows <- sample(seq_len(nrow(eligible)), size = TEST_DONORS, replace = FALSE)
selected <- eligible[selected_rows]
selected[, holdout_index := seq_len(.N)]
selected[, holdout_label := sprintf("holdout_%02d", holdout_index)]
# Não exportamos donor_id: só índices e contagens por tipo celular.
donor_selection_out <- selected[, c("holdout_index", "holdout_label", CELL_TYPES), with = FALSE]
fwrite(donor_selection_out, file.path(result_dir, "06f_donor_selection.tsv"), sep = "\t")

# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================
standardize <- function(raw, sample_ids, label) {
  if (is.null(raw)) stop(label, ": campo de proporções ausente.")
  p <- as.matrix(raw)
  if (!is.null(rownames(p)) && all(sample_ids %in% rownames(p))) {
    p <- p[sample_ids, , drop = FALSE]
  } else if (!is.null(colnames(p)) && all(sample_ids %in% colnames(p))) {
    p <- t(p[, sample_ids, drop = FALSE])
  } else stop(label, ": amostras não identificadas nas proporções.")
  if (!setequal(colnames(p), CELL_TYPES)) stop(label, ": oito tipos não encontrados.")
  p <- p[, CELL_TYPES, drop = FALSE]
  storage.mode(p) <- "double"
  if (anyNA(p) || any(!is.finite(p)) || any(p < -1e-8)) stop(label, ": proporções inválidas.")
  if (max(abs(rowSums(p) - 1)) > SUM_TOL) stop(label, ": proporções não somam 1.")
  p
}

# Gerar perfis por tipo como MÉDIA de counts por núcleo do doador retido.
# Misturas são EXPECTATIVAS determinísticas e não requerem reamostragem de
# 1.000 núcleos reais, que não existem nos tipos raros.
make_synthetic_counts <- function(reference_counts, reference_meta, donor) {
  means <- matrix(NA_real_, nrow = nrow(reference_counts), ncol = length(CELL_TYPES),
                  dimnames = list(rownames(reference_counts), CELL_TYPES))
  for (ct in CELL_TYPES) {
    ix <- which(reference_meta$donor_id == donor & reference_meta$cell_type == ct)
    if (!length(ix)) stop("Doação ausente para ", donor, " / ", ct)
    means[, ct] <- as.numeric(Matrix::rowMeans(reference_counts[, ix, drop = FALSE]))
  }
  synthetic <- means %*% t(virtual_n)
  synthetic <- round(synthetic)
  colnames(synthetic) <- rownames(compositions)
  storage.mode(synthetic) <- "double"
  if (anyNA(synthetic) || any(synthetic < 0)) stop("Counts sintéticos inválidos.")
  synthetic
}

# ==============================================================================
# 4. LEAVE-ONE-DONOR-OUT: TREINO vs TESTE
# ==============================================================================
proportion_results <- list()
fit_audit <- list()
fit_i <- 1L
for (holdout_i in seq_len(nrow(selected))) {
  donor <- selected$donor_id[holdout_i]
  label <- selected$holdout_label[holdout_i]
  cat("\n", label, " | gerando mistura de doador não presente no treino...\n", sep = "")
  ref_train_idx <- which(meta$donor_id != donor)
  ref_test_idx <- which(meta$donor_id == donor)
  if (length(intersect(ref_train_idx, ref_test_idx))) stop("Vazamento de núcleos treino/teste.")
  train_meta <- meta[ref_train_idx]
  if (uniqueN(train_meta$donor_id) != 23L ||
      !setequal(unique(train_meta$cell_type), CELL_TYPES)) {
    stop("Treinamento deveria incluir 23 doadores e oito tipos.")
  }
  test_mix <- make_synthetic_counts(ref, meta, donor)
  sample_ids <- paste(label, colnames(test_mix), sep = "_")
  colnames(test_mix) <- sample_ids
  true <- compositions
  rownames(true) <- sample_ids

  for (gene_set in GENE_SETS) {
    geneids <- if (gene_set == "top300") {
      rownames(ref)[rownames(ref) %in% marker_geneids]
    } else rownames(ref)

    sc <- ref[geneids, ref_train_idx, drop = FALSE]
    bulk <- test_mix[geneids, , drop = FALSE]
    genes_keep <- Matrix::rowSums(sc) > 0 & rowSums(bulk) > 0
    genes_used <- geneids[genes_keep]
    if (length(genes_used) < 100L) stop("Poucos genes em ", label, "/", gene_set)
    sc <- sc[genes_used, , drop = FALSE]
    bulk <- bulk[genes_used, , drop = FALSE]
    coldata <- S4Vectors::DataFrame(
      cellType = factor(train_meta$cell_type, levels = CELL_TYPES),
      sampleID = factor(train_meta$donor_id),
      row.names = as.character(train_meta$TAG)
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = sc), colData = coldata
    )
    rownames(sce) <- genes_used
    colnames(sce) <- as.character(train_meta$TAG)
    if (!identical(rownames(sce), rownames(bulk))) stop("GeneIDs desalinhados.")

    cat("  ", gene_set, " | ", length(genes_used), " genes...\n", sep = "")
    begin <- Sys.time()
    fit <- MuSiC::music_prop(
      bulk.mtx = bulk,
      sc.sce = sce,
      clusters = "cellType",
      samples = "sampleID",
      select.ct = CELL_TYPES,
      markers = NULL,
      cell_size = NULL,
      ct.cov = FALSE,
      iter.max = 1000,
      nu = 0.0001,
      eps = 0.01,
      centered = FALSE,
      normalize = FALSE,
      verbose = FALSE
    )
    seconds <- as.numeric(difftime(Sys.time(), begin, units = "secs"))
    for (method in c("MuSiC_weighted", "NNLS_unweighted")) {
      raw <- if (method == "MuSiC_weighted") fit$Est.prop.weighted else fit$Est.prop.allgene
      estimates <- standardize(raw, sample_ids, paste(label, gene_set, method))
      x <- as.data.table(estimates, keep.rownames = "sample_id")
      long <- data.table::melt(x, id.vars = "sample_id", variable.name = "cell_type",
                               value.name = "estimated_fraction")
      long[, cell_type := as.character(cell_type)]
      truth <- as.data.table(true, keep.rownames = "sample_id")
      truth <- data.table::melt(truth, id.vars = "sample_id", variable.name = "cell_type",
                                value.name = "true_fraction")
      truth[, cell_type := as.character(cell_type)]
      long <- merge(long, truth, by = c("sample_id", "cell_type"), sort = FALSE)
      long[, `:=`(holdout_index = holdout_i, mixture = sub(paste0("^", label, "_"), "", sample_id),
                  gene_set = gene_set, method = method)]
      long[, `:=`(error_pp = 100 * (estimated_fraction - true_fraction),
                  absolute_error_pp = 100 * abs(estimated_fraction - true_fraction),
                  estimated_zero = estimated_fraction <= ZERO_TOL)]
      setcolorder(long, c("holdout_index", "mixture", "gene_set", "method", "sample_id",
                          "cell_type", "true_fraction", "estimated_fraction", "error_pp",
                          "absolute_error_pp", "estimated_zero"))
      proportion_results[[length(proportion_results) + 1L]] <- long
    }
    fit_audit[[fit_i]] <- data.table(
      holdout_index = holdout_i, gene_set = gene_set,
      training_donors = uniqueN(train_meta$donor_id),
      training_nuclei = length(ref_train_idx), test_nuclei_observed = length(ref_test_idx),
      synthetic_samples = ncol(bulk), genes_used = length(genes_used),
      mean_library_size_synthetic = mean(colSums(bulk)),
      elapsed_seconds = seconds,
      feature_selection_leakage = gene_set == "top300"
    )
    fit_i <- fit_i + 1L
    fwrite(rbindlist(proportion_results), file.path(result_dir, "06f_proportions_partial.tsv"), sep = "\t")
    fwrite(rbindlist(fit_audit), file.path(result_dir, "06f_fit_audit_partial.tsv"), sep = "\t")
    rm(sc, bulk, sce, fit)
    invisible(gc())
  }
  rm(test_mix, true)
  invisible(gc())
}

# ==============================================================================
# 5. MÉTRICAS DE RECUPERAÇÃO E QC
# ==============================================================================
p <- rbindlist(proportion_results, use.names = TRUE, fill = TRUE)
fit_dt <- rbindlist(fit_audit, use.names = TRUE, fill = TRUE)
by_type <- p[, .(
  n_mixtures = .N,
  true_fraction_min = min(true_fraction), true_fraction_max = max(true_fraction),
  mean_true_fraction = mean(true_fraction), mean_estimated_fraction = mean(estimated_fraction),
  mean_absolute_error_pp = mean(absolute_error_pp),
  median_absolute_error_pp = median(absolute_error_pp),
  median_bias_pp = median(error_pp),
  fraction_estimated_zero = mean(estimated_zero)
), by = .(gene_set, method, cell_type)]
by_run <- p[, .(
  n_predictions = .N,
  mean_absolute_error_pp = mean(absolute_error_pp),
  root_mean_squared_error_pp = sqrt(mean(error_pp^2)),
  max_absolute_error_pp = max(absolute_error_pp),
  fraction_estimated_zero_for_positive_truth = mean(estimated_zero & true_fraction > 0)
), by = .(gene_set, method)]
status <- data.table(
  check = c("donors_held_out", "all_eight_types_in_heldout", "all_mixture_proportions_positive",
            "all_expected_fits_completed", "no_train_test_nucleus_overlap",
            "finite_estimates", "fraction_sums_one", "marker_leakage_documented"),
  status = c(
    if (nrow(selected) == TEST_DONORS) "PASS" else "FAIL",
    if (all(as.matrix(selected[, ..CELL_TYPES]) > 0L)) "PASS" else "FAIL",
    if (all(p$true_fraction > 0)) "PASS" else "FAIL",
    if (nrow(fit_dt) == TEST_DONORS * length(GENE_SETS)) "PASS" else "FAIL",
    "PASS", # verificação executada dentro de cada fold antes de executar o ajuste
    if (all(is.finite(p$estimated_fraction))) "PASS" else "FAIL",
    if (max(abs(p[, .(s = sum(estimated_fraction)),
                  by = .(holdout_index, mixture, gene_set, method)]$s - 1)) <= SUM_TOL) "PASS" else "FAIL",
    if (all(fit_dt[gene_set == "top300", feature_selection_leakage])) "PASS" else "FAIL"
  ),
  note = c(
    "Os doadores usados para criar o bulk sintético não são usados no treino do fold.",
    "End e Per foram realmente observados nos doadores de teste.",
    "A verdade conhecida é a composição de NÚCLEOS VIRTUAIS, não composição bulk externa.",
    "Quatro doadores; allgenes e top300 para cada doador.",
    "Nenhum núcleo do doador de teste está na referência de treino do fold.",
    NA_character_, NA_character_,
    "Marcadores top300 foram descobertos com todos os 24 doadores: avaliar allgenes como teste principal."
  )
)
fwrite(p, file.path(result_dir, "06f_proportions_long.tsv"), sep = "\t", na = "NA")
fwrite(fit_dt, file.path(result_dir, "06f_fit_audit.tsv"), sep = "\t", na = "NA")
fwrite(by_type, file.path(result_dir, "06f_metrics_by_type.tsv"), sep = "\t", na = "NA")
fwrite(by_run, file.path(result_dir, "06f_metrics_by_run.tsv"), sep = "\t", na = "NA")
fwrite(status, file.path(result_dir, "06f_status.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(result_dir, "06f_sessionInfo.txt"))

# ==============================================================================
# 6. GRÁFICOS
# ==============================================================================
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  colors <- c("MuSiC_weighted" = "#00ADA2", "NNLS_unweighted" = "#7B2CBF")
  g1 <- ggplot(p, aes(x = 100 * true_fraction, y = 100 * estimated_fraction, color = method)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#64748B") +
    geom_point(alpha = 0.65, size = 1.7) +
    facet_grid(cell_type ~ gene_set, scales = "free") +
    scale_color_manual(values = colors) +
    labs(title = "06f — recuperação de misturas sintéticas Mathys",
         subtitle = "Doador de teste fora da referência; linha pontilhada = recuperação perfeita",
         x = "Composição conhecida (%)", y = "Estimativa MuSiC (%)",
         caption = "Misturas ideais da mesma plataforma: não validam desempenho em bulk real") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(plot_dir, "06f_observed_vs_expected.png"), g1, width = 13, height = 11, dpi = 180)
  g2 <- ggplot(by_type, aes(x = cell_type, y = mean_absolute_error_pp, fill = method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ gene_set) +
    scale_fill_manual(values = colors) +
    labs(title = "06f — erro absoluto médio por tipo celular",
         subtitle = "Quatro doadores excluídos individualmente; 5 misturas por doador",
         x = "Tipo celular", y = "Erro médio absoluto (pontos percentuais)",
         caption = "Pior desempenho de In/Opc na simulação reforçaria hipótese de não separabilidade") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.title = element_blank())
  ggsave(file.path(plot_dir, "06f_mean_absolute_error.png"), g2, width = 11, height = 5.5, dpi = 180)
} else message("ggplot2 não instalado; tabelas foram geradas sem gráficos.")

cat("\n=== 06f: VALIDAÇÃO SINTÉTICA CONCLUÍDA ===\n")
print(fit_dt[, .(holdout_index, gene_set, training_donors, genes_used, elapsed_seconds)])
cat("\nErros globais:\n")
print(by_run)
cat("\nResultados de In/Opc:\n")
print(by_type[cell_type %in% c("In", "Opc")])
cat("\nStatus:\n")
print(status[, .(check, status)])
cat("\nSaídas em: ", result_dir, "\n", sep = "")
if (any(status$status == "FAIL")) stop("06f possui FAIL; revisar os outputs.")
