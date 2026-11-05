#!/usr/bin/env Rscript
# =============================================================================
# 06c_audit_MuSiC_vs_MuSiC2.R
# =============================================================================
# Fase 02 | Auditoria diagnóstica do MuSiC2 com MuSiC convencional.
#
# PROBLEMA
#   In e Opc são estimados como zero em quase todas as amostras do MuSiC2,
#   inclusive em CONTROLES. Adicionalmente, o MuSiC2 não converge
#   completamente em nenhuma das cinco execuções.
#
# OBJETIVO
#   (a) Reexecutar APENAS MuSiC convencional nas MESMAS entradas top300,
#       referência Mathys e amostras usadas pelo MuSiC2.
#   (b) Executar separadamente Control e Case, tal como a inicialização do
#       próprio music2_prop(). Isso permite comparar o componente inicial
#       (MuSiC convencional) ao retorno final MuSiC2, sem reexecutar MuSiC2.
#   (c) Inspecionar zeros em In/Opc, alterações por tipo, fit e qualidade da
#       referência. Não remover amostras, células ou genes por conveniência.
#
# DOCUMENTAÇÃO: MuSiC::music_prop() retorna $Est.prop.weighted.
#   xuranw.github.io/MuSiC/reference/music_prop.html
#   github.com/xuranw/MuSiC/blob/master/R/music2.R
#
# EXECUTAR DA RAIZ DO REPOSITÓRIO:
#   Rscript scripts/phase02_deconvolution/06c_audit_MuSiC_vs_MuSiC2.R
#
# PRÉ-REQUISITO: etapas 03, 05 e 06 já executadas. Nenhum arquivo anterior
# é substituído: novos resultados vão para results/deconvolution/MuSiC_Audit/.
#
# PARÂMETROS DE INÍCIO:
#   DATASETS = ROSMAP e GSE125583 para comparação principal/validação.
#   Posteriormente, adicionar GSE64810 e GSE68719.
#   RUN_ROSMAP_NO_OVERLAP21 = FALSE (opcional, após diagnóstico inicial).
#   O universo de entrada é estritamente top300, igual ao script 06.
#
# IMPORTANTE: isto NÃO é um estudo de acurácia com ground truth celular. É
# uma auditoria de consistência entre métodos, não validação externa.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. CONFIGURAÇÃO
# -----------------------------------------------------------------------------
DATASETS <- c("ROSMAP", "GSE125583")
# Para completar a auditoria nas quatro coortes, alterar para:
# DATASETS <- c("ROSMAP", "GSE125583", "GSE64810", "GSE68719")
RUN_ROSMAP_NO_OVERLAP21 <- FALSE

CELL_TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
SEED <- 1234L
EXACT_ZERO <- 0
NEAR_ZERO <- 0.001  # 0.1% em escala de proporção
SUM_TOL <- 1e-6
MATCH_CONTROLS_TOL <- 1e-6

PKGS <- c("data.table", "Matrix", "SingleCellExperiment", "S4Vectors", "MuSiC")
not_installed <- PKGS[!vapply(PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(not_installed)) {
  stop("Pacotes ausentes: ", paste(not_installed, collapse = ", "))
}
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(MuSiC)
})

OUT <- file.path("results", "deconvolution", "MuSiC_Audit")
MUSIC2_DIR <- file.path("results", "deconvolution", "MuSiC2")
BULK_DIR <- file.path("data", "processed", "deconvolution")
REF_DIR <- file.path("data", "reference", "mathys2019", "phase02_prepared", "reference")
REF_COUNTS <- file.path(REF_DIR, "mathys2019_pathology_free_geneid_counts_sparse.rds")
REF_META <- file.path(REF_DIR, "mathys2019_pathology_free_cell_metadata_aligned.tsv")
MUSIC2_LONG <- file.path(MUSIC2_DIR, "06_proportions_long_all.tsv")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. FUNÇÕES AUXILIARES
# -----------------------------------------------------------------------------
check_exists <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Arquivos ausentes:\n", paste(missing, collapse = "\n"))
}

num_mat <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}

# music_prop retorna amostras x tipos em Est.prop.weighted. A função lida
# também com saída transposta em caso de variação de versões do pacote.
standardize_prop <- function(result, sample_ids, label) {
  p <- result$Est.prop.weighted
  if (is.null(p)) stop(label, ": music_prop não retornou Est.prop.weighted.")
  p <- num_mat(p)
  if (!is.null(rownames(p)) && all(sample_ids %in% rownames(p))) {
    p <- p[sample_ids, , drop = FALSE]
  } else if (!is.null(colnames(p)) && all(sample_ids %in% colnames(p))) {
    p <- t(p[, sample_ids, drop = FALSE])
  } else {
    stop(label, ": impossível alinhar o Est.prop.weighted às amostras.")
  }
  missing_ct <- setdiff(CELL_TYPES, colnames(p))
  if (length(missing_ct)) {
    stop(label, ": MuSiC retornou menos de 8 tipos (",
         paste(missing_ct, collapse = ", "), "). Não serão imputados zeros.")
  }
  p <- p[, CELL_TYPES, drop = FALSE]
  if (any(!is.finite(p)) || any(p < -1e-8)) {
    stop(label, ": proporções não finitas ou negativas.")
  }
  deviation <- max(abs(rowSums(p) - 1))
  if (deviation > SUM_TOL) {
    warning(label, ": max |soma-1| = ", signif(deviation, 5),
            ". Valores preservados, sem renormalização.")
  }
  p
}

# Extrai as estimativas MuSiC2 diretamente do resumo detalhado da etapa 06;
# uma única linha por amostra x tipo x versão da referência.
read_previous_music2 <- function(data, dataset, reference_version) {
  ds_target <- dataset
  ref_target <- reference_version
  x <- data[dataset == ds_target & reference_version == ref_target &
              marker_set == "top300"]
  if (!nrow(x)) stop("Sem MuSiC2 prévio: ", dataset, " / ", reference_version)
  if (anyDuplicated(x[, .(sample_id, cell_type)])) {
    stop("Dados MuSiC2 duplicados: ", dataset)
  }
  x[, .(sample_id, group, cell_type, music2 = proportion)]
}

# Avalia não só zeros, mas diferenças nas mesmas amostras.
summarize_method <- function(comp) {
  ans <- comp[, .(
    n_samples = .N,
    music_mean = mean(music),
    music2_mean = mean(music2),
    music_zero_n = sum(music == EXACT_ZERO),
    music2_zero_n = sum(music2 == EXACT_ZERO),
    music_zero_fraction = mean(music == EXACT_ZERO),
    music2_zero_fraction = mean(music2 == EXACT_ZERO),
    music_near_zero_fraction = mean(music < NEAR_ZERO),
    music2_near_zero_fraction = mean(music2 < NEAR_ZERO),
    mean_absolute_difference_pp = 100 * mean(abs(music - music2)),
    median_difference_pp = 100 * median(music2 - music),
    spearman = if (sd(music) > 0 && sd(music2) > 0) {
      suppressWarnings(cor(music, music2, method = "spearman"))
    } else {
      NA_real_
    }
  ), by = .(dataset, reference_version, group, cell_type)]
  ans
}

# -----------------------------------------------------------------------------
# 2. PRÉ-VALIDAÇÃO
# -----------------------------------------------------------------------------
check_exists(c(REF_COUNTS, REF_META, MUSIC2_LONG))
if (!exists("music_prop", envir = asNamespace("MuSiC"), inherits = FALSE)) {
  stop("O pacote MuSiC instalado não contém music_prop().")
}
required_args <- c("bulk.mtx", "sc.sce", "clusters", "samples", "select.ct",
                   "markers", "cell_size", "ct.cov", "iter.max", "nu", "eps",
                   "centered", "normalize")
actual_args <- names(formals(MuSiC::music_prop))
if (!all(required_args %in% actual_args)) {
  stop("API MuSiC::music_prop diferente da esperada. Verifique versão instalada.")
}

music2_all <- fread(MUSIC2_LONG, na.strings = c("NA", ""))
required_cols <- c("dataset", "reference_version", "marker_set", "sample_id",
                   "group", "cell_type", "proportion")
if (!all(required_cols %in% names(music2_all))) {
  stop("06_proportions_long_all.tsv sem colunas obrigatórias.")
}

reference_counts <- readRDS(REF_COUNTS)
if (!inherits(reference_counts, "dgCMatrix")) reference_counts <- as(reference_counts, "dgCMatrix")
ref_metadata <- fread(REF_META, na.strings = c("NA", ""))
needed_ref_cols <- c("TAG", "donor_id", "cell_type", "overlap_with_bulk_ROSMAP")
if (!all(needed_ref_cols %in% names(ref_metadata))) stop("Metadados Mathys incompletos.")
if (anyDuplicated(ref_metadata$TAG) || anyDuplicated(colnames(reference_counts))) {
  stop("TAGs duplicadas na referência.")
}
ref_i <- match(colnames(reference_counts), ref_metadata$TAG)
if (anyNA(ref_i)) stop("Referência: núcleos sem metadata.")
ref_metadata <- ref_metadata[ref_i]
stopifnot(identical(as.character(ref_metadata$TAG), colnames(reference_counts)))
if (!setequal(ref_metadata$cell_type, CELL_TYPES)) stop("Tipos celulares Mathys inesperados.")

fwrite(data.table(
  parameter = c("MuSiC_version", "R_version", "top300_only", "match_music2_settings",
                "iter.max", "nu", "eps", "ct.cov", "centered", "normalize",
                "datasets", "RUN_ROSMAP_NO_OVERLAP21"),
  value = c(as.character(packageVersion("MuSiC")), as.character(getRversion()),
            "TRUE", "TRUE", "1000", "0.0001", "0.01", "FALSE", "FALSE", "FALSE",
            paste(DATASETS, collapse = ";"), as.character(RUN_ROSMAP_NO_OVERLAP21))
), file.path(OUT, "06c_parameters.tsv"), sep = "\t", na = "NA")

# -----------------------------------------------------------------------------
# 3. EXECUÇÕES DIAGNÓSTICAS
# -----------------------------------------------------------------------------
plan <- data.table(dataset = DATASETS, reference_version = "main24")
if (RUN_ROSMAP_NO_OVERLAP21 && "ROSMAP" %in% DATASETS) {
  plan <- rbind(plan, data.table(dataset = "ROSMAP",
                                reference_version = "no_ROSMAP_overlap21"))
}

all_comp <- list()
all_summary <- list()
all_status <- list()
all_support <- list()
all_similarity <- list()

for (i in seq_len(nrow(plan))) {
  ds <- plan$dataset[i]
  ref_version <- plan$reference_version[i]
  disease <- switch(ds, ROSMAP = "AD", GSE125583 = "AD",
                    GSE64810 = "HD", GSE68719 = "PD",
                    stop("Dataset desconhecido: ", ds))

  cat("\n==================================================\n")
  cat("MuSiC diagnóstico: ", ds, " / ", ref_version, " / top300\n", sep = "")

  cfile <- file.path(BULK_DIR, ds, "bulk_counts_top300.rds")
  mfile <- file.path(BULK_DIR, ds, "metadata_deconvolution.tsv")
  check_exists(c(cfile, mfile))
  bulk <- num_mat(readRDS(cfile))
  bulk_meta <- fread(mfile, na.strings = c("NA", ""))
  if (!all(c("sample_id", "group") %in% names(bulk_meta))) {
    stop(ds, ": metadata sem sample_id/group.")
  }
  idx <- match(colnames(bulk), bulk_meta$sample_id)
  if (anyNA(idx)) stop(ds, ": bulk com samples sem metadata.")
  bulk_meta <- bulk_meta[idx]
  if (anyDuplicated(bulk_meta$sample_id)) stop(ds, ": sample_ids duplicados.")
  if (!all(bulk_meta$group %in% c("Control", disease))) {
    stop(ds, ": grupos diferentes de Control/", disease)
  }

  if (ref_version == "main24") {
    keep_cell <- rep(TRUE, ncol(reference_counts))
  } else {
    ov <- toupper(as.character(ref_metadata$overlap_with_bulk_ROSMAP))
    keep_cell <- !(ov %in% c("TRUE", "T", "1"))
    if (any(is.na(ov) | ov == "NA")) stop("Overlap Mathys com valores ausentes.")
  }
  sc_meta <- ref_metadata[keep_cell]
  n_donors <- uniqueN(sc_meta$donor_id)
  if ((ref_version == "main24" && n_donors != 24L) ||
      (ref_version == "no_ROSMAP_overlap21" && n_donors != 21L)) {
    stop(ds, ": quantidade inesperada de doadores em ", ref_version)
  }

  # Interseção e filtro exatamente como em 06_prepare_run_inputs().
  g <- rownames(bulk)[rownames(bulk) %in% rownames(reference_counts)]
  if (length(g) < 100L) stop(ds, ": poucos genes em comum.")
  sc <- reference_counts[g, keep_cell, drop = FALSE]
  bulk <- bulk[g, , drop = FALSE]
  keep_gene <- rowSums(bulk) > 0 & Matrix::rowSums(sc) > 0
  g <- g[keep_gene]
  bulk <- bulk[g, , drop = FALSE]
  sc <- sc[g, , drop = FALSE]
  if (length(g) < 100L) stop(ds, ": poucos genes após all-zero QC.")
  cat("Genes após QC: ", length(g), " | Núcleos: ", ncol(sc),
      " | Doadores: ", n_donors, "\n", sep = "")

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = sc),
    colData = S4Vectors::DataFrame(
      cellType = factor(sc_meta$cell_type, levels = CELL_TYPES),
      sampleID = factor(sc_meta$donor_id),
      row.names = sc_meta$TAG
    )
  )
  colnames(sce) <- as.character(sc_meta$TAG)
  rownames(sce) <- g

  # A auditoria da referência independe da inferência e ajuda a explicar zeros.
  support <- sc_meta[, .(n_nuclei = .N, n_donors = uniqueN(donor_id)), by = cell_type]
  support[, `:=`(dataset = ds, reference_version = ref_version)]
  all_support[[i]] <- support

  # Correlação entre perfis médios (indicativo de similaridade de assinatura,
  # NÃO avaliação de exatidão das proporções). Registra abundância total média.
  mean_counts <- vapply(CELL_TYPES, function(ct) {
    Matrix::rowMeans(sc[, sc_meta$cell_type == ct, drop = FALSE])
  }, numeric(nrow(sc)))
  rownames(mean_counts) <- g
  cell_corr <- suppressWarnings(stats::cor(log1p(mean_counts), method = "spearman"))
  sim <- as.data.table(as.table(cell_corr))
  setnames(sim, c("cell_type_1", "cell_type_2", "spearman_logmean"))
  sim[, `:=`(dataset = ds, reference_version = ref_version)]
  all_similarity[[i]] <- sim

  # Comparação direta com inicialização MuSiC dentro de music2_prop():
  # O MuSiC2 chama music_prop() separadamente nos controles e nos casos.
  result_by_group <- list()
  for (grp in c("Control", disease)) {
    ids <- bulk_meta[group == grp, sample_id]
    if (length(ids) < 2L) stop("Poucas amostras: ", ds, " / ", grp)
    b <- bulk[, ids, drop = FALSE]
    cat("  Rodando MuSiC convencional — ", grp, " (", length(ids), ")...\n", sep = "")
    tic <- Sys.time()
    set.seed(SEED + i)
    # Mesma configuração usada internamente no music2_prop() para a
    # estimativa inicial. Não executa o ciclo iterativo de MuSiC2.
    fit <- MuSiC::music_prop(
      bulk.mtx = b,
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
    elapsed <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
    result_by_group[[grp]] <- standardize_prop(fit, ids, paste(ds, grp))
    cat("    MuSiC finalizou em ", round(elapsed, 1), " s\n", sep = "")
    rm(fit)
    invisible(gc())
  }

  music <- rbind(result_by_group[["Control"]], result_by_group[[disease]])
  music2 <- read_previous_music2(music2_all, ds, ref_version)
  md <- melt(as.data.table(music, keep.rownames = "sample_id"),
             id.vars = "sample_id", variable.name = "cell_type", value.name = "music")
  md[, cell_type := as.character(cell_type)]
  md <- merge(md, bulk_meta[, .(sample_id, group)], by = "sample_id", all.x = TRUE)
  comparison <- merge(md, music2, by = c("sample_id", "group", "cell_type"),
                      all = FALSE, sort = FALSE)
  expect_n <- nrow(bulk_meta) * length(CELL_TYPES)
  if (nrow(comparison) != expect_n || anyDuplicated(comparison[, .(sample_id, cell_type)])) {
    stop(ds, ": junção MuSiC/MuSiC2 incompleta ou duplicada. Não usar resultado.")
  }
  comparison[, `:=`(
    dataset = ds,
    reference_version = ref_version,
    difference = music2 - music,
    absolute_difference_pp = 100 * abs(music2 - music)
  )]
  all_comp[[i]] <- comparison
  all_summary[[i]] <- summarize_method(comparison)

  # Controle: no music2_prop() os controles provêm do MuSiC convencional.
  ctrl <- comparison[group == "Control"]
  max_control_diff <- max(abs(ctrl$music2 - ctrl$music))
  nonzero_mu <- comparison[, .(n_nonzero_music = sum(music > 0),
                               n_nonzero_music2 = sum(music2 > 0)),
                           by = .(group, cell_type)]
  all_status[[i]] <- data.table(
    dataset = ds,
    reference_version = ref_version,
    gene_count = length(g),
    reference_donors = n_donors,
    bulk_samples = nrow(bulk_meta),
    control_max_absolute_diff = max_control_diff,
    control_matches_music2_initialization = max_control_diff <= MATCH_CONTROLS_TOL,
    In_nonzero_music = sum(comparison[cell_type == "In", music > 0]),
    In_nonzero_music2 = sum(comparison[cell_type == "In", music2 > 0]),
    Opc_nonzero_music = sum(comparison[cell_type == "Opc", music > 0]),
    Opc_nonzero_music2 = sum(comparison[cell_type == "Opc", music2 > 0]),
    note = if (max_control_diff <= MATCH_CONTROLS_TOL) {
      "MuSiC de controle reproduz MuSiC2: verificar a etapa inicial antes da iteração."
    } else {
      "MuSiC e MuSiC2 diferem nos controles: AUDITAR entradas/API antes de interpretar."
    }
  )
  cat("  Máx. diferença nos controles: ", signif(max_control_diff, 4), "\n", sep = "")
  if (max_control_diff > MATCH_CONTROLS_TOL) {
    warning(ds, ": MuSiC controle NÃO reproduz MuSiC2. Examine inputs e versões.")
  }
  run_dir <- file.path(OUT, ds, ref_version, "top300")
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(comparison, file.path(run_dir, "06c_comparison_by_sample.tsv"), sep = "\t", na = "NA")
  fwrite(all_summary[[i]], file.path(run_dir, "06c_celltype_summary.tsv"), sep = "\t", na = "NA")
  fwrite(nonzero_mu, file.path(run_dir, "06c_nonzero_counts.tsv"), sep = "\t", na = "NA")

  rm(bulk, sce, sc, sc_meta, mean_counts, cell_corr, music,
     result_by_group, comparison)
  invisible(gc())
}

# -----------------------------------------------------------------------------
# 4. CONSOLIDAR E SALVAR
# -----------------------------------------------------------------------------
comp_all <- rbindlist(all_comp, fill = TRUE, use.names = TRUE)
summary_all <- rbindlist(all_summary, fill = TRUE, use.names = TRUE)
status_all <- rbindlist(all_status, fill = TRUE, use.names = TRUE)
support_all <- unique(rbindlist(all_support, fill = TRUE, use.names = TRUE))
sim_all <- unique(rbindlist(all_similarity, fill = TRUE, use.names = TRUE))

fwrite(comp_all, file.path(OUT, "06c_comparison_long_all.tsv"), sep = "\t", na = "NA")
fwrite(summary_all, file.path(OUT, "06c_celltype_summary_all.tsv"), sep = "\t", na = "NA")
fwrite(status_all, file.path(OUT, "06c_status.tsv"), sep = "\t", na = "NA")
fwrite(support_all, file.path(OUT, "06c_reference_support.tsv"), sep = "\t", na = "NA")
fwrite(sim_all, file.path(OUT, "06c_reference_profile_similarity.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(OUT, "06c_sessionInfo.txt"))

# -----------------------------------------------------------------------------
# 5. VISUALIZAÇÕES DA AUDITORIA (opcionais, não modificam as estimativas)
# -----------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  vis_dir <- file.path(OUT, "plots")
  dir.create(vis_dir, recursive = TRUE, showWarnings = FALSE)
  assay_colors <- c("MuSiC" = "#277DA1", "MuSiC2" = "#D55E80")
  group_colors <- c("Control" = "#94A3B8", "AD" = "#7B2CBF",
                    "HD" = "#277DA1", "PD" = "#E75480")
  common_theme <- theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", colour = "#141937"),
          panel.grid.minor = element_blank())

  # Proporção de resultados exatamente zero: cada método, cada tipo, cada grupo.
  z <- summary_all[cell_type %in% c("In", "Opc")]
  z <- melt(z,
            id.vars = c("dataset", "reference_version", "group", "cell_type"),
            measure.vars = c("music_zero_fraction", "music2_zero_fraction"),
            variable.name = "method", value.name = "zero_fraction")
  z[, method := factor(ifelse(method == "music_zero_fraction", "MuSiC", "MuSiC2"),
                       levels = c("MuSiC", "MuSiC2"))]
  z[, condition_label := paste(dataset, reference_version, group, sep = " / ")]
  p_zero <- ggplot(z, aes(x = condition_label, y = 100 * zero_fraction,
                         fill = method)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    scale_fill_manual(values = assay_colors) +
    coord_flip() + facet_wrap(~cell_type, ncol = 2) +
    labs(title = "MuSiC convencional vs MuSiC2 — estimativas zero",
         x = NULL, y = "Estimativas iguais a zero (%)", fill = "Método",
         caption = "Comparação técnica usando a mesma referência, top300 e amostras.") +
    common_theme
  ggsave(file.path(vis_dir, "06c_01_zeros_MuSiC_vs_MuSiC2.png"),
         p_zero, width = 12, height = 5.8, dpi = 180)

  # Concordância em cada coorte: uma amostra representa um ponto.
  for (ds in unique(comp_all$dataset)) {
    dd <- comp_all[dataset == ds & reference_version == "main24"]
    if (!nrow(dd)) next
    p_cmp <- ggplot(dd, aes(x = 100 * music, y = 100 * music2,
                            colour = group)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  colour = "#64748B", linewidth = 0.4) +
      geom_point(alpha = 0.55, size = 1.5) +
      scale_colour_manual(values = group_colors) +
      facet_wrap(~cell_type, scales = "free", ncol = 4) +
      labs(title = paste0("Concordância MuSiC inicial vs MuSiC2: ", ds),
           subtitle = "Cada ponto corresponde à mesma amostra e tipo celular",
           x = "MuSiC convencional (%)", y = "MuSiC2 (%)", colour = "Grupo",
           caption = "Diagonal = estimativas idênticas. Escalas livres por tipo; não são estimativas validadas.") +
      common_theme
    ggsave(file.path(vis_dir, paste0("06c_02_concordancia_", ds, ".png")),
           p_cmp, width = 12.5, height = 7.1, dpi = 180)
  }
} else {
  message("ggplot2 não instalado: relatórios tabulares gerados, sem gráficos 06c.")
}

cat("\n=======================================================\n")
cat("Auditoria MuSiC convencional concluída\n")
cat("=======================================================\n\n")
print(status_all)
cat("\nResumo In / Opc por coorte e grupo:\n")
print(summary_all[cell_type %in% c("In", "Opc"), .(
  dataset, reference_version, group, cell_type, n_samples,
  music_zero_fraction, music2_zero_fraction,
  mean_absolute_difference_pp
)])
cat("\nResultados salvos em: ", OUT, "\n", sep = "")
cat("OBS: MuSiC não possui o mesmo campo de convergência global do MuSiC2.\n")
