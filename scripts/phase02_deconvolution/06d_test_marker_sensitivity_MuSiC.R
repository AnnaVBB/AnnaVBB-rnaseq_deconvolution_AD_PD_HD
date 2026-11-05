#!/usr/bin/env Rscript
# =============================================================================
# 06d_test_marker_sensitivity_MuSiC.R
# FASE 02 — TESTE DIAGNÓSTICO DA SELEÇÃO DE MARCADORES
# =============================================================================
# OBJETIVO:
# Investigar se proporções quase zero de In e Opc já na inicialização do MuSiC
# dependem do tamanho do conjunto de genes (top100 / top300 / top500).
#
# - Executa MuSiC convencional nos CONTROLES de ROSMAP e GSE125583 para
#   top100 e top500, mantendo a mesma referência e parâmetros do 06c.
# - Reutiliza as estimativas já salvas em 06 para top300 (controles), sem nova
#   execução: 06c comprovou identidade numérica com a inicialização do MuSiC2.
# - Preserva todos os oito tipos celulares. Não imputa genes nem proporções.
# - Não muda a escolha de coorte nem valida biologicamente as estimativas.
#
# EXECUTAR A PARTIR DA RAIZ DO PROJETO:
# Rscript scripts/phase02_deconvolution/06d_test_marker_sensitivity_MuSiC.R
#
# OUTPUT:
# results/deconvolution/MuSiC_Marker_Sensitivity/06d_*.tsv
# results/deconvolution/MuSiC_Marker_Sensitivity/plots/
# =============================================================================

DATASETS <- c("ROSMAP", "GSE125583")
MARKER_SETS_TO_RUN <- c("top100", "top500")
REFERENCE_VERSION <- "main24"
CELL_TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
NEAR_ZERO <- 0.001           # proporção 0.001 = 0,1%
ZERO_TOL <- 0                # zero exatamente igual a zero
REFERENCE_DONORS_EXPECTED <- 24L
REFERENCE_NUCLEI_EXPECTED <- 35140L

required <- c("data.table", "Matrix", "MuSiC", "S4Vectors", "SingleCellExperiment")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Instale pacotes ausentes: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(S4Vectors)
  library(SingleCellExperiment)
  library(MuSiC)
})

out_dir <- file.path("results", "deconvolution", "MuSiC_Marker_Sensitivity")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

ref_dir <- file.path("data", "reference", "mathys2019", "phase02_prepared", "reference")
ref_counts_path <- file.path(ref_dir, "mathys2019_pathology_free_geneid_counts_sparse.rds")
ref_meta_path <- file.path(ref_dir, "mathys2019_pathology_free_cell_metadata_aligned.tsv")
bulk_root <- file.path("data", "processed", "deconvolution")
mu2_path <- file.path("results", "deconvolution", "MuSiC2", "06_proportions_long_all.tsv")

check_paths <- c(ref_counts_path, ref_meta_path, mu2_path)
if (!all(file.exists(check_paths))) {
  stop("Arquivos ausentes: ", paste(check_paths[!file.exists(check_paths)], collapse = "; "))
}

message("Carregando referência Mathys 2019...")
ref <- readRDS(ref_counts_path)
if (!inherits(ref, "dgCMatrix")) ref <- as(ref, "dgCMatrix")
ref_meta <- fread(ref_meta_path, na.strings = c("NA", ""))
if (!all(c("TAG", "donor_id", "cell_type") %in% names(ref_meta))) {
  stop("Metadados da referência sem TAG, donor_id ou cell_type.")
}
if (anyDuplicated(ref_meta$TAG) || anyDuplicated(colnames(ref))) {
  stop("TAG duplicado na referência.")
}
ridx <- match(colnames(ref), ref_meta$TAG)
if (anyNA(ridx)) stop("Colunas da referência sem metadados correspondentes.")
ref_meta <- ref_meta[ridx]
if (!identical(as.character(ref_meta$TAG), colnames(ref))) {
  stop("Falha no alinhamento da referência.")
}
if (uniqueN(ref_meta$donor_id) != REFERENCE_DONORS_EXPECTED ||
    ncol(ref) != REFERENCE_NUCLEI_EXPECTED) {
  stop("Números de doadores/núcleos diferem da referência main24 esperada.")
}
if (!setequal(unique(ref_meta$cell_type), CELL_TYPES)) {
  stop("Tipos celulares diferentes dos oito previstos.")
}

previous <- fread(mu2_path, na.strings = c("NA", ""))
needed_previous <- c("dataset", "reference_version", "marker_set", "group",
                     "sample_id", "cell_type", "proportion")
if (!all(needed_previous %in% names(previous))) stop("MuSiC2 consolidated sem colunas esperadas.")

standardize_music <- function(fit, sample_ids, label) {
  p <- fit$Est.prop.weighted
  if (is.null(p)) stop(label, ": music_prop não retornou Est.prop.weighted.")
  p <- as.matrix(p)
  if (!is.null(rownames(p)) && all(sample_ids %in% rownames(p))) {
    p <- p[sample_ids, , drop = FALSE]
  } else if (!is.null(colnames(p)) && all(sample_ids %in% colnames(p))) {
    p <- t(p[, sample_ids, drop = FALSE])
  } else {
    stop(label, ": não foi possível alinhar as amostras na saída.")
  }
  if (!setequal(colnames(p), CELL_TYPES)) stop(label, ": saída sem os oito tipos.")
  p <- p[, CELL_TYPES, drop = FALSE]
  storage.mode(p) <- "double"
  if (anyNA(p) || any(!is.finite(p)) || any(p < -1e-8)) stop(label, ": proporções inválidas.")
  if (max(abs(rowSums(p) - 1)) > 1e-6) stop(label, ": soma das proporções fora da tolerância.")
  p
}

calc_stats <- function(x) {
  x[, .(
    n_samples = .N,
    mean_proportion = mean(proportion),
    median_proportion = median(proportion),
    n_zero = sum(proportion == ZERO_TOL),
    fraction_zero = mean(proportion == ZERO_TOL),
    n_near_zero = sum(proportion < NEAR_ZERO),
    fraction_near_zero = mean(proportion < NEAR_ZERO),
    n_nonzero = sum(proportion > 0)
  ), by = .(dataset, marker_set, cell_type)]
}

all_prop <- list()
all_inputs <- list()
row_counter <- 1L

for (ds in DATASETS) {
  disease <- "AD"
  bulk_meta_path <- file.path(bulk_root, ds, "metadata_deconvolution.tsv")
  if (!file.exists(bulk_meta_path)) stop("Metadata não encontrado: ", bulk_meta_path)
  bulk_meta <- fread(bulk_meta_path, na.strings = c("NA", ""))
  if (!all(c("sample_id", "group") %in% names(bulk_meta))) stop(ds, ": metadata incompleto.")
  if (anyDuplicated(bulk_meta$sample_id)) stop(ds, ": sample IDs duplicados.")
  control_ids <- bulk_meta[group == "Control", sample_id]
  if (length(control_ids) < 2L) stop(ds, ": controles insuficientes.")
  message("\nDataset: ", ds, " | controles: ", length(control_ids))

  # MuSiC2 usa MuSiC convencional para os controles. A auditoria 06c confirmou
  # identidade numérica, portanto reutilizamos este top300 como controle-base.
  baseline <- previous[
    dataset == ds & reference_version == REFERENCE_VERSION &
      marker_set == "top300" & group == "Control",
    .(sample_id, cell_type, proportion)
  ]
  if (nrow(baseline) != length(control_ids) * length(CELL_TYPES) ||
      anyDuplicated(baseline[, .(sample_id, cell_type)]) ||
      !setequal(baseline$sample_id, control_ids)) {
    stop(ds, ": top300 controle de MuSiC2 ausente/inconsistente.")
  }
  baseline[, `:=`(dataset = ds, marker_set = "top300")]
  all_prop[[row_counter]] <- baseline
  row_counter <- row_counter + 1L

  for (mk in MARKER_SETS_TO_RUN) {
    message("  MuSiC convencional ", mk, ": iniciando...")
    bulk_path <- file.path(bulk_root, ds, paste0("bulk_counts_", mk, ".rds"))
    if (!file.exists(bulk_path)) stop("Bulk não encontrado: ", bulk_path)
    bulk <- as.matrix(readRDS(bulk_path))
    storage.mode(bulk) <- "double"
    if (is.null(rownames(bulk)) || is.null(colnames(bulk))) stop("Matriz bulk sem dimnames.")
    if (anyDuplicated(rownames(bulk)) || anyDuplicated(colnames(bulk))) {
      stop("Matriz bulk com dimnames duplicados.")
    }
    if (anyNA(bulk) || any(!is.finite(bulk)) || any(bulk < 0)) stop("Counts bulk inválidos.")
    if (!all(control_ids %in% colnames(bulk))) stop(ds, ": controles não encontrados no bulk.")
    common <- rownames(bulk)[rownames(bulk) %in% rownames(ref)]
    if (length(common) < 100L) stop(ds, ": poucos genes comuns.")
    b <- bulk[common, control_ids, drop = FALSE]
    sc <- ref[common, , drop = FALSE]
    # Para comparabilidade com 06: o filtro all-zero no bulk é calculado usando
    # TODAS as amostras Control+AD, e não só os controles rodados neste 06d.
    b_full <- bulk[common, , drop = FALSE]
    keep <- rowSums(b_full) > 0 & Matrix::rowSums(sc) > 0
    g <- common[keep]
    if (length(g) < 100L) stop(ds, " / ", mk, ": <100 genes após QC.")
    b <- b[g, , drop = FALSE]
    sc <- sc[g, , drop = FALSE]
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = sc),
      colData = S4Vectors::DataFrame(
        cellType = factor(ref_meta$cell_type, levels = CELL_TYPES),
        sampleID = factor(ref_meta$donor_id),
        row.names = ref_meta$TAG
      )
    )
    rownames(sce) <- g
    colnames(sce) <- ref_meta$TAG
    started <- Sys.time()
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
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    p <- standardize_music(fit, control_ids, paste(ds, mk))
    long <- melt(as.data.table(p, keep.rownames = "sample_id"), id.vars = "sample_id",
                 variable.name = "cell_type", value.name = "proportion")
    long[, cell_type := as.character(cell_type)]
    long[, `:=`(dataset = ds, marker_set = mk)]
    all_prop[[row_counter]] <- long
    row_counter <- row_counter + 1L
    all_inputs[[length(all_inputs) + 1L]] <- data.table(
      dataset = ds, marker_set = mk, n_controls = length(control_ids),
      genes_used = length(g), reference_donors = uniqueN(ref_meta$donor_id),
      reference_nuclei = ncol(sc), elapsed_seconds = elapsed
    )
    message("  ", mk, " concluído: ", length(g), " genes / ", round(elapsed, 1), " segundos")
    rm(bulk, b, b_full, sc, sce, fit, p, long)
    invisible(gc())
  }
}

proportions <- rbindlist(all_prop, use.names = TRUE, fill = TRUE)
setcolorder(proportions, c("dataset", "marker_set", "sample_id", "cell_type", "proportion"))
if (anyDuplicated(proportions[, .(dataset, marker_set, sample_id, cell_type)])) {
  stop("Combinações duplicadas na tabela final.")
}
expected_rows <- sum(vapply(DATASETS, function(ds) {
  n_samples <- uniqueN(proportions[dataset == ds & marker_set == "top300", sample_id])
  as.numeric(n_samples * length(CELL_TYPES) * (length(MARKER_SETS_TO_RUN) + 1L))
}, numeric(1)))
if (nrow(proportions) != expected_rows) stop("Número de linhas final inesperado.")

summary_dt <- calc_stats(proportions)
setorder(summary_dt, dataset, cell_type, marker_set)

baseline <- proportions[marker_set == "top300", .(
  dataset, sample_id, cell_type, baseline_top300 = proportion
)]
compare <- merge(
  proportions[marker_set != "top300"], baseline,
  by = c("dataset", "sample_id", "cell_type"),
  all.x = TRUE, sort = FALSE
)
if (anyNA(compare$baseline_top300)) stop("Pareamento com top300 incompleto.")
comparison_summary <- compare[, .(
  n_samples = .N,
  mean_absolute_difference_pp = 100 * mean(abs(proportion - baseline_top300)),
  median_difference_pp = 100 * median(proportion - baseline_top300),
  spearman = if (sd(proportion) > 0 && sd(baseline_top300) > 0) {
    suppressWarnings(cor(proportion, baseline_top300, method = "spearman"))
  } else {
    NA_real_
  }
), by = .(dataset, marker_set, cell_type)]

input_summary <- if (length(all_inputs)) rbindlist(all_inputs) else data.table()
fwrite(proportions, file.path(out_dir, "06d_proportions_controls_long.tsv"), sep = "\t", na = "NA")
fwrite(summary_dt, file.path(out_dir, "06d_celltype_summary.tsv"), sep = "\t", na = "NA")
fwrite(comparison_summary, file.path(out_dir, "06d_comparison_to_top300.tsv"), sep = "\t", na = "NA")
fwrite(input_summary, file.path(out_dir, "06d_input_audit.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "06d_sessionInfo.txt"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  colors <- c(top100 = "#00ADA2", top300 = "#64748B", top500 = "#7B2CBF")
  df <- summary_dt[cell_type %in% c("Ex", "In", "Oli", "Opc")]
  df[, marker_set := factor(marker_set, levels = c("top100", "top300", "top500"))]
  p <- ggplot(df, aes(x = marker_set, y = 100 * fraction_zero, fill = marker_set)) +
    geom_col(width = 0.75) +
    facet_grid(cell_type ~ dataset) +
    scale_fill_manual(values = colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(
      title = "MuSiC convencional: frequência de zeros por conjunto de marcadores",
      subtitle = "Apenas controles; top300 obtido da inicialização do MuSiC2, validada em 06c",
      x = "Conjunto", y = "Proporções iguais a zero (%)", fill = "Marcadores",
      caption = "Diagnóstico de robustez, não estimativa de acurácia da composição celular"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  ggplot2::ggsave(file.path(plot_dir, "06d_zero_fraction_top100_300_500.png"),
                  p, width = 10.6, height = 8.2, dpi = 180)
}

cat("\n=== AUDITORIA MuSiC POR CONJUNTO DE MARCADORES CONCLUÍDA ===\n")
print(summary_dt[cell_type %in% c("Ex", "In", "Oli", "Opc"), .(
  dataset, marker_set, cell_type, n_samples, n_nonzero,
  fraction_zero, fraction_near_zero, mean_proportion
)])
cat("\nArquivos em: ", out_dir, "\n", sep = "")
cat("Atenção: não houve validação externa e nenhuma coorte/tipo foi excluída.\n")
