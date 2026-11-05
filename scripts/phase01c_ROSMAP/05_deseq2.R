# ==============================================================================
# 05_deseq2_ROSMAP.R
# ==============================================================================
# Fase 01c — ROSMAP
# Expressão diferencial com DESeq2
#
# Entradas:
#   results/candidate_ROSMAP/qc/counts_ROSMAP_GeneID_filtered.tsv.gz
#   results/candidate_ROSMAP/metadata/metadata_ROSMAP.tsv
#   results/candidate_ROSMAP/qc/gene_annotation_ROSMAP_filtered.tsv
#   results/candidate_ROSMAP/qc/vst_ROSMAP_filtered.rds
#
# Modelos:
#   Principal: ~ group
#   Ajustado:  ~ sequencingBatch + msex + RIN + pmi + group
#
# Implementação do modelo ajustado:
#   ~ sequencingBatch_factor + msex_factor + RIN_centered +
#     pmi_centered + group_factor
#
# Critérios pré-definidos:
#   Principal: BH-FDR < 0,05 e |log2FC| > 0,58
#   Estrito:   BH-FDR < 0,05 e |log2FC| > 1
#   Nominal:   p < 0,05 e |log2FC| > 1
#
# Control é a referência:
#   log2FC > 0 → maior expressão em AD
#   log2FC < 0 → menor expressão em AD
#
# A matriz VST é usada somente para visualização/heatmap.
# O DESeq2 recebe raw counts filtrados.
# ==============================================================================


# ==============================================================================
# 0. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table",
  "DESeq2",
  "SummarizedExperiment",
  "ggplot2",
  "pheatmap"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Pacotes necessários não instalados: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
  library(pheatmap)
})


# ==============================================================================
# 1. PARÂMETROS
# ==============================================================================

EXPECTED_MAIN_N <- 144L
EXPECTED_AD_N <- 77L
EXPECTED_CONTROL_N <- 67L

PRIMARY_ALPHA <- 0.05
PRIMARY_LFC <- 0.58

STRICT_ALPHA <- 0.05
STRICT_LFC <- 1.00

NOMINAL_P <- 0.05
NOMINAL_LFC <- 1.00

HEATMAP_TOP_N <- 50L
VOLCANO_LABEL_N <- 12L


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

candidate_dir <- file.path("results", "candidate_ROSMAP")
qc_dir <- file.path(candidate_dir, "qc")
metadata_dir <- file.path(candidate_dir, "metadata")
deseq2_dir <- file.path(candidate_dir, "deseq2")
plot_dir <- file.path(deseq2_dir, "plots")

dir.create(deseq2_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

counts_input <- file.path(
  qc_dir,
  "counts_ROSMAP_GeneID_filtered.tsv.gz"
)

metadata_input <- file.path(
  metadata_dir,
  "metadata_ROSMAP.tsv"
)

annotation_input <- file.path(
  qc_dir,
  "gene_annotation_ROSMAP_filtered.tsv"
)

vst_input <- file.path(
  qc_dir,
  "vst_ROSMAP_filtered.rds"
)

main_results_output <- file.path(
  deseq2_dir,
  "ROSMAP_deseq2_results.tsv"
)

adjusted_results_output <- file.path(
  deseq2_dir,
  "ROSMAP_deseq2_results_adjusted.tsv"
)

main_primary_output <- file.path(
  deseq2_dir,
  "ROSMAP_DEG_primary.tsv"
)

main_strict_output <- file.path(
  deseq2_dir,
  "ROSMAP_DEG_strict.tsv"
)

main_nominal_output <- file.path(
  deseq2_dir,
  "ROSMAP_DEG_nominal_comparative.tsv"
)

adjusted_primary_output <- file.path(
  deseq2_dir,
  "ROSMAP_DEG_primary_adjusted.tsv"
)

adjusted_strict_output <- file.path(
  deseq2_dir,
  "ROSMAP_DEG_strict_adjusted.tsv"
)

adjusted_nominal_output <- file.path(
  deseq2_dir,
  "ROSMAP_DEG_nominal_comparative_adjusted.tsv"
)

model_summary_output <- file.path(
  deseq2_dir,
  "ROSMAP_deseq2_model_summary.tsv"
)

stability_gene_output <- file.path(
  deseq2_dir,
  "ROSMAP_model_stability_gene_level.tsv"
)

stability_summary_output <- file.path(
  deseq2_dir,
  "ROSMAP_model_stability_summary.tsv"
)

design_matrix_output <- file.path(
  deseq2_dir,
  "ROSMAP_adjusted_design_matrix.tsv"
)

normalized_counts_output <- file.path(
  deseq2_dir,
  "ROSMAP_normalized_counts_main.tsv.gz"
)

cooks_output <- file.path(
  deseq2_dir,
  "ROSMAP_Cooks_distance_by_sample.tsv"
)

top_genes_output <- file.path(
  deseq2_dir,
  "ROSMAP_top_genes_for_heatmap.tsv"
)

top_genes_adjusted_output <- file.path(
  deseq2_dir,
  "ROSMAP_top_genes_for_heatmap_adjusted.tsv"
)

dds_main_output <- file.path(
  deseq2_dir,
  "ROSMAP_dds_main.rds"
)

dds_adjusted_output <- file.path(
  deseq2_dir,
  "ROSMAP_dds_adjusted.rds"
)

session_info_output <- file.path(
  deseq2_dir,
  "ROSMAP_deseq2_sessionInfo.txt"
)

report_output <- file.path(
  deseq2_dir,
  "05_deseq2_ROSMAP.txt"
)

ma_plot_output <- file.path(
  plot_dir,
  "MA_ROSMAP.png"
)

volcano_plot_output <- file.path(
  plot_dir,
  "volcano_ROSMAP.png"
)

volcano_adjusted_plot_output <- file.path(
  plot_dir,
  "volcano_ROSMAP_adjusted.png"
)

heatmap_plot_output <- file.path(
  plot_dir,
  "heatmap_ROSMAP_top50_primary.png"
)

heatmap_adjusted_plot_output <- file.path(
  plot_dir,
  "heatmap_ROSMAP_top50_primary_adjusted.png"
)

stability_plot_output <- file.path(
  plot_dir,
  "log2FC_main_vs_adjusted_ROSMAP.png"
)

dispersion_main_output <- file.path(
  plot_dir,
  "dispersion_ROSMAP_main.png"
)

dispersion_adjusted_output <- file.path(
  plot_dir,
  "dispersion_ROSMAP_adjusted.png"
)


# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================

assert_file <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo não encontrado:\n", path)
  }
  invisible(TRUE)
}

assert_columns <- function(dt, required, object_name) {
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop(
      object_name,
      ": colunas obrigatórias ausentes: ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}

fmt_int <- function(x) {
  format(
    round(as.numeric(x), 0),
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

fmt_num <- function(x, digits = 4) {
  format(
    round(as.numeric(x), digits),
    nsmall = digits,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

count_fractional_values <- function(x) {
  sum(abs(x - round(x)) > 1e-8, na.rm = TRUE)
}

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

write_tsv_gzip <- function(dt, path_gz) {
  if (!grepl("\\.gz$", path_gz)) {
    stop("write_tsv_gzip exige caminho terminado em .gz.")
  }

  path_plain <- sub("\\.gz$", "", path_gz)

  fwrite(
    dt,
    path_plain,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  gzip_bin <- Sys.which("gzip")

  if (identical(gzip_bin, "")) {
    stop(
      "gzip não encontrado. Arquivo sem compressão salvo em:\n",
      path_plain
    )
  }

  status <- system2(
    gzip_bin,
    args = c("-f", path_plain)
  )

  if (status != 0 || !file.exists(path_gz)) {
    stop("Falha ao criar: ", path_gz)
  }

  invisible(path_gz)
}


# ------------------------------------------------------------------------------
# Classificar resultados segundo os três critérios definidos.
# ------------------------------------------------------------------------------

classify_results <- function(x) {

  x[
    ,
    significant_FDR :=
      !is.na(padj) &
      padj < PRIMARY_ALPHA
  ]

  x[
    ,
    DEG :=
      significant_FDR &
      abs(log2FoldChange) > PRIMARY_LFC
  ]

  x[
    ,
    direction := fifelse(
      DEG & log2FoldChange > PRIMARY_LFC,
      "Up",
      fifelse(
        DEG & log2FoldChange < -PRIMARY_LFC,
        "Down",
        "Not_DEG"
      )
    )
  ]

  x[
    ,
    DEG_strict :=
      !is.na(padj) &
      padj < STRICT_ALPHA &
      abs(log2FoldChange) > STRICT_LFC
  ]

  x[
    ,
    direction_strict := fifelse(
      DEG_strict & log2FoldChange > STRICT_LFC,
      "Up",
      fifelse(
        DEG_strict & log2FoldChange < -STRICT_LFC,
        "Down",
        "Not_DEG"
      )
    )
  ]

  x[
    ,
    DEG_nominal_comparative :=
      !is.na(pvalue) &
      pvalue < NOMINAL_P &
      abs(log2FoldChange) > NOMINAL_LFC
  ]

  x[
    ,
    direction_nominal := fifelse(
      DEG_nominal_comparative & log2FoldChange > NOMINAL_LFC,
      "Up",
      fifelse(
        DEG_nominal_comparative & log2FoldChange < -NOMINAL_LFC,
        "Down",
        "Not_DEG"
      )
    )
  ]

  x
}


summarize_model <- function(x, model_name) {

  min_p <- if (all(is.na(x$pvalue))) {
    NA_real_
  } else {
    min(x$pvalue, na.rm = TRUE)
  }

  min_q <- if (all(is.na(x$padj))) {
    NA_real_
  } else {
    min(x$padj, na.rm = TRUE)
  }

  data.table(
    model = model_name,
    n_genes = nrow(x),
    n_pvalue_nonmissing = sum(!is.na(x$pvalue)),
    n_pvalue_missing = sum(is.na(x$pvalue)),
    n_padj_nonmissing = sum(!is.na(x$padj)),
    n_padj_missing = sum(is.na(x$padj)),
    FDR_lt_0.05 = sum(x$significant_FDR, na.rm = TRUE),
    primary_total = sum(x$DEG, na.rm = TRUE),
    primary_up = sum(x$direction == "Up", na.rm = TRUE),
    primary_down = sum(x$direction == "Down", na.rm = TRUE),
    strict_total = sum(x$DEG_strict, na.rm = TRUE),
    strict_up = sum(x$direction_strict == "Up", na.rm = TRUE),
    strict_down = sum(x$direction_strict == "Down", na.rm = TRUE),
    nominal_total = sum(x$DEG_nominal_comparative, na.rm = TRUE),
    nominal_up = sum(x$direction_nominal == "Up", na.rm = TRUE),
    nominal_down = sum(x$direction_nominal == "Down", na.rm = TRUE),
    min_pvalue = min_p,
    min_padj = min_q
  )
}


# ==============================================================================
# 4. VALIDAR ENTRADAS
# ==============================================================================

invisible(
  lapply(
    c(
      counts_input,
      metadata_input,
      annotation_input,
      vst_input
    ),
    assert_file
  )
)


# ==============================================================================
# 5. METADATA E VARIÁVEIS DO MODELO
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat("05 — DESeq2 ROSMAP\n")
cat(
  "============================================================\n\n"
)

metadata <- fread(
  metadata_input,
  na.strings = c("NA", "")
)

assert_columns(
  metadata,
  c(
    "sample_id",
    "individual_id",
    "group",
    "cogdx",
    "sequencingBatch",
    "msex",
    "RIN",
    "pmi"
  ),
  "metadata_ROSMAP.tsv"
)

metadata[
  ,
  sample_id := trimws(as.character(sample_id))
]

metadata[
  ,
  individual_id := trimws(as.character(individual_id))
]

metadata[
  ,
  group := trimws(as.character(group))
]

metadata[
  ,
  sequencingBatch := trimws(as.character(sequencingBatch))
]

if (nrow(metadata) != EXPECTED_MAIN_N) {
  stop("Metadata principal deveria conter 144 amostras.")
}

if (anyDuplicated(metadata$sample_id)) {
  stop("metadata contém sample_id duplicado.")
}

if (anyDuplicated(metadata$individual_id)) {
  stop("metadata contém individual_id duplicado.")
}

group_counts <- metadata[, .N, by = group]

n_ad <- group_counts[group == "AD", N]
n_control <- group_counts[group == "Control", N]

if (length(n_ad) != 1 || n_ad != EXPECTED_AD_N) {
  stop("Número inesperado de amostras AD.")
}

if (length(n_control) != 1 || n_control != EXPECTED_CONTROL_N) {
  stop("Número inesperado de controles.")
}


# ------------------------------------------------------------------------------
# Control como referência.
# ------------------------------------------------------------------------------

metadata[
  ,
  group_factor := factor(
    group,
    levels = c("Control", "AD")
  )
]


# ------------------------------------------------------------------------------
# Batch técnico como variável categórica.
# Os rótulos oficiais são preservados sem interpretação adicional.
# ------------------------------------------------------------------------------

metadata[
  ,
  sequencingBatch_factor := factor(sequencingBatch)
]


# ------------------------------------------------------------------------------
# msex como variável categórica 0/1.
# Não traduzir os níveis antes de confirmar o codebook.
# ------------------------------------------------------------------------------

metadata[
  ,
  msex_factor := factor(as.character(msex))
]


# ------------------------------------------------------------------------------
# RIN e PMI como covariáveis contínuas centralizadas.
# ------------------------------------------------------------------------------

metadata[
  ,
  RIN_centered :=
    as.numeric(RIN) -
    mean(as.numeric(RIN), na.rm = TRUE)
]

metadata[
  ,
  pmi_centered :=
    as.numeric(pmi) -
    mean(as.numeric(pmi), na.rm = TRUE)
]

adjusted_vars <- c(
  "group_factor",
  "sequencingBatch_factor",
  "msex_factor",
  "RIN_centered",
  "pmi_centered"
)

if (!all(complete.cases(metadata[, ..adjusted_vars]))) {
  stop("Existem dados ausentes nas covariáveis do modelo ajustado.")
}


# ==============================================================================
# 6. COUNTS FILTRADOS
# ==============================================================================

cat("Lendo counts filtrados...\n")

counts_dt <- fread(
  counts_input,
  na.strings = c("NA", "")
)

assert_columns(
  counts_dt,
  "GeneID",
  "counts_ROSMAP_GeneID_filtered.tsv.gz"
)

counts_dt[
  ,
  GeneID := trimws(as.character(GeneID))
]

if (anyDuplicated(counts_dt$GeneID)) {
  stop("GeneID duplicado na matriz filtrada.")
}

if (
  anyNA(counts_dt$GeneID) ||
  any(counts_dt$GeneID == "")
) {
  stop("GeneID ausente/vazio na matriz filtrada.")
}

if (
  !identical(
    names(counts_dt)[-1],
    metadata$sample_id
  )
) {
  stop("Ordem counts ↔ metadata não corresponde.")
}

count_matrix <- as.matrix(
  counts_dt[, -1, with = FALSE]
)

rownames(count_matrix) <- counts_dt$GeneID
storage.mode(count_matrix) <- "numeric"

if (
  sum(is.na(count_matrix)) > 0 ||
  sum(!is.finite(count_matrix), na.rm = TRUE) > 0 ||
  sum(count_matrix < 0, na.rm = TRUE) > 0 ||
  count_fractional_values(count_matrix) > 0
) {
  stop("A matriz de entrada do DESeq2 apresenta problema de integridade.")
}

if (any(rowSums(count_matrix) == 0)) {
  stop("Genes all-zero permaneceram após o filtro.")
}

if (any(colSums(count_matrix) == 0)) {
  stop("Existe amostra all-zero.")
}

if (max(count_matrix) > .Machine$integer.max) {
  stop("Existe count individual acima do limite integer do R.")
}

storage.mode(count_matrix) <- "integer"


# ==============================================================================
# 7. ANOTAÇÃO
# ==============================================================================

annotation <- fread(
  annotation_input,
  na.strings = c("NA", "")
)

assert_columns(
  annotation,
  "GeneID",
  "gene_annotation_ROSMAP_filtered.tsv"
)

annotation[
  ,
  GeneID := as.character(GeneID)
]

if (anyDuplicated(annotation$GeneID)) {
  stop("A anotação contém GeneID duplicado.")
}

annotation <- annotation[
  match(rownames(count_matrix), GeneID)
]

if (
  anyNA(annotation$GeneID) ||
  !identical(annotation$GeneID, rownames(count_matrix))
) {
  stop("A anotação não corresponde exatamente aos GeneIDs da matriz.")
}


# ==============================================================================
# 8. colData E MATRIZ DE DESIGN
# ==============================================================================

coldata <- as.data.frame(metadata)
rownames(coldata) <- metadata$sample_id

if (
  !identical(
    rownames(coldata),
    colnames(count_matrix)
  )
) {
  stop("rownames(colData) != colnames(count_matrix).")
}

adjusted_formula <- ~
  sequencingBatch_factor +
  msex_factor +
  RIN_centered +
  pmi_centered +
  group_factor

design_matrix <- model.matrix(
  adjusted_formula,
  data = coldata
)

design_rank <- qr(design_matrix)$rank
design_ncol <- ncol(design_matrix)

if (design_rank != design_ncol) {
  stop(
    "Design ajustado sem posto completo. rank=",
    design_rank,
    "; ncol=",
    design_ncol
  )
}

fwrite(
  as.data.table(
    design_matrix,
    keep.rownames = "sample_id"
  ),
  design_matrix_output,
  sep = "\t"
)


# ==============================================================================
# 9. MODELO PRINCIPAL: ~ group
# ==============================================================================

cat("Executando modelo principal ~ group...\n")

dds_main <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ group_factor
)

dds_main <- DESeq(
  dds_main,
  quiet = TRUE
)

if (
  !"group_factor_AD_vs_Control" %in%
  resultsNames(dds_main)
) {
  stop(
    "Coeficiente AD vs Control não encontrado no modelo principal."
  )
}

res_main <- results(
  dds_main,
  contrast = c(
    "group_factor",
    "AD",
    "Control"
  ),
  alpha = PRIMARY_ALPHA
)


# ==============================================================================
# 10. MODELO AJUSTADO
# ==============================================================================

cat("Executando modelo ajustado...\n")

dds_adjusted <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = adjusted_formula
)

dds_adjusted <- DESeq(
  dds_adjusted,
  quiet = TRUE
)

if (
  !"group_factor_AD_vs_Control" %in%
  resultsNames(dds_adjusted)
) {
  stop(
    "Coeficiente AD vs Control não encontrado no modelo ajustado."
  )
}

res_adjusted <- results(
  dds_adjusted,
  contrast = c(
    "group_factor",
    "AD",
    "Control"
  ),
  alpha = PRIMARY_ALPHA
)


# ==============================================================================
# 11. RESULTADOS + ANOTAÇÃO
# ==============================================================================

res_main_dt <- as.data.table(
  as.data.frame(res_main),
  keep.rownames = "GeneID"
)

res_adjusted_dt <- as.data.table(
  as.data.frame(res_adjusted),
  keep.rownames = "GeneID"
)

res_main_dt[
  ,
  GeneID := as.character(GeneID)
]

res_adjusted_dt[
  ,
  GeneID := as.character(GeneID)
]

res_main_dt <- merge(
  res_main_dt,
  annotation,
  by = "GeneID",
  all.x = TRUE,
  sort = FALSE
)

res_adjusted_dt <- merge(
  res_adjusted_dt,
  annotation,
  by = "GeneID",
  all.x = TRUE,
  sort = FALSE
)

res_main_dt <- res_main_dt[
  match(rownames(count_matrix), GeneID)
]

res_adjusted_dt <- res_adjusted_dt[
  match(rownames(count_matrix), GeneID)
]

if (
  anyNA(res_main_dt$GeneID) ||
  anyNA(res_adjusted_dt$GeneID)
) {
  stop("Falha ao restaurar a ordem dos resultados.")
}

res_main_dt <- classify_results(res_main_dt)
res_adjusted_dt <- classify_results(res_adjusted_dt)

if ("symbol" %in% names(res_main_dt)) {

  res_main_dt[
    ,
    plot_label := fifelse(
      !is.na(symbol) & symbol != "",
      symbol,
      GeneID
    )
  ]

  res_adjusted_dt[
    ,
    plot_label := fifelse(
      !is.na(symbol) & symbol != "",
      symbol,
      GeneID
    )
  ]

} else {

  res_main_dt[, plot_label := GeneID]
  res_adjusted_dt[, plot_label := GeneID]
}


# ==============================================================================
# 12. SALVAR RESULTADOS E LISTAS
# ==============================================================================

fwrite(
  res_main_dt,
  main_results_output,
  sep = "\t",
  na = "NA"
)

fwrite(
  res_adjusted_dt,
  adjusted_results_output,
  sep = "\t",
  na = "NA"
)

fwrite(
  res_main_dt[
    DEG == TRUE
  ][
    order(padj, -abs(log2FoldChange))
  ],
  main_primary_output,
  sep = "\t",
  na = "NA"
)

fwrite(
  res_main_dt[
    DEG_strict == TRUE
  ][
    order(padj, -abs(log2FoldChange))
  ],
  main_strict_output,
  sep = "\t",
  na = "NA"
)

fwrite(
  res_main_dt[
    DEG_nominal_comparative == TRUE
  ][
    order(pvalue, -abs(log2FoldChange))
  ],
  main_nominal_output,
  sep = "\t",
  na = "NA"
)

fwrite(
  res_adjusted_dt[
    DEG == TRUE
  ][
    order(padj, -abs(log2FoldChange))
  ],
  adjusted_primary_output,
  sep = "\t",
  na = "NA"
)

fwrite(
  res_adjusted_dt[
    DEG_strict == TRUE
  ][
    order(padj, -abs(log2FoldChange))
  ],
  adjusted_strict_output,
  sep = "\t",
  na = "NA"
)

fwrite(
  res_adjusted_dt[
    DEG_nominal_comparative == TRUE
  ][
    order(pvalue, -abs(log2FoldChange))
  ],
  adjusted_nominal_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 13. RESUMO DOS MODELOS
# ==============================================================================

model_summary <- rbind(
  summarize_model(
    res_main_dt,
    "main_group_only"
  ),
  summarize_model(
    res_adjusted_dt,
    "adjusted_batch_msex_RIN_PMI"
  ),
  fill = TRUE
)

fwrite(
  model_summary,
  model_summary_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 14. ESTABILIDADE PRINCIPAL × AJUSTADO
# ==============================================================================

stability <- merge(
  res_main_dt[
    ,
    .(
      GeneID,
      symbol = if ("symbol" %in% names(res_main_dt)) symbol else NA_character_,
      log2FC_main = log2FoldChange,
      pvalue_main = pvalue,
      padj_main = padj,
      DEG_main = DEG,
      direction_main = direction,
      DEG_strict_main = DEG_strict,
      DEG_nominal_main = DEG_nominal_comparative
    )
  ],
  res_adjusted_dt[
    ,
    .(
      GeneID,
      log2FC_adjusted = log2FoldChange,
      pvalue_adjusted = pvalue,
      padj_adjusted = padj,
      DEG_adjusted = DEG,
      direction_adjusted = direction,
      DEG_strict_adjusted = DEG_strict,
      DEG_nominal_adjusted = DEG_nominal_comparative
    )
  ],
  by = "GeneID",
  all = FALSE,
  sort = FALSE
)

stability[
  ,
  primary_overlap_class := fifelse(
    DEG_main & DEG_adjusted,
    "Shared_primary",
    fifelse(
      DEG_main & !DEG_adjusted,
      "Main_only",
      fifelse(
        !DEG_main & DEG_adjusted,
        "Adjusted_only",
        "Neither"
      )
    )
  )
]

stability[
  ,
  same_effect_sign :=
    sign(log2FC_main) ==
    sign(log2FC_adjusted)
]

stability[
  ,
  sign_switch :=
    is.finite(log2FC_main) &
    is.finite(log2FC_adjusted) &
    sign(log2FC_main) != sign(log2FC_adjusted)
]

stability[
  ,
  shared_primary_same_direction :=
    DEG_main &
    DEG_adjusted &
    direction_main == direction_adjusted
]

stability[
  ,
  shared_primary_opposite_direction :=
    DEG_main &
    DEG_adjusted &
    direction_main != direction_adjusted
]

spearman_all <- safe_cor(
  stability$log2FC_main,
  stability$log2FC_adjusted,
  "spearman"
)

pearson_all <- safe_cor(
  stability$log2FC_main,
  stability$log2FC_adjusted,
  "pearson"
)

n_shared_primary <- stability[
  primary_overlap_class == "Shared_primary",
  .N
]

n_main_only <- stability[
  primary_overlap_class == "Main_only",
  .N
]

n_adjusted_only <- stability[
  primary_overlap_class == "Adjusted_only",
  .N
]

n_shared_same_direction <- stability[
  shared_primary_same_direction == TRUE,
  .N
]

n_shared_opposite_direction <- stability[
  shared_primary_opposite_direction == TRUE,
  .N
]

n_sign_switch_all <- stability[
  sign_switch == TRUE,
  .N
]

primary_union <- stability[
  DEG_main | DEG_adjusted,
  .N
]

primary_jaccard <- if (primary_union == 0) {
  NA_real_
} else {
  n_shared_primary / primary_union
}

stability_primary_union <- stability[
  DEG_main | DEG_adjusted
]

spearman_primary_union <- safe_cor(
  stability_primary_union$log2FC_main,
  stability_primary_union$log2FC_adjusted,
  "spearman"
)

pearson_primary_union <- safe_cor(
  stability_primary_union$log2FC_main,
  stability_primary_union$log2FC_adjusted,
  "pearson"
)

stability_summary <- data.table(
  metric = c(
    "genes_compared",
    "spearman_log2FC_all_genes",
    "pearson_log2FC_all_genes",
    "shared_primary_DEGs",
    "main_only_primary_DEGs",
    "adjusted_only_primary_DEGs",
    "primary_union",
    "primary_Jaccard",
    "shared_primary_same_direction",
    "shared_primary_opposite_direction",
    "all_genes_sign_switch",
    "spearman_log2FC_primary_union",
    "pearson_log2FC_primary_union"
  ),
  value = c(
    nrow(stability),
    spearman_all,
    pearson_all,
    n_shared_primary,
    n_main_only,
    n_adjusted_only,
    primary_union,
    primary_jaccard,
    n_shared_same_direction,
    n_shared_opposite_direction,
    n_sign_switch_all,
    spearman_primary_union,
    pearson_primary_union
  )
)

fwrite(
  stability,
  stability_gene_output,
  sep = "\t",
  na = "NA"
)

fwrite(
  stability_summary,
  stability_summary_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 15. NORMALIZED COUNTS E COOK'S DISTANCE
# ==============================================================================

normalized_counts <- counts(
  dds_main,
  normalized = TRUE
)

write_tsv_gzip(
  as.data.table(
    normalized_counts,
    keep.rownames = "GeneID"
  ),
  normalized_counts_output
)

cooks_main <- assays(dds_main)[["cooks"]]

if (is.null(cooks_main)) {

  cooks_sample <- data.table(
    sample_id = metadata$sample_id,
    max_Cooks_distance = NA_real_,
    median_Cooks_distance = NA_real_,
    group = metadata$group,
    sequencingBatch = metadata$sequencingBatch
  )

} else {

  cooks_sample <- data.table(
    sample_id = colnames(cooks_main),
    max_Cooks_distance = apply(
      cooks_main,
      2,
      max,
      na.rm = TRUE
    ),
    median_Cooks_distance = apply(
      cooks_main,
      2,
      median,
      na.rm = TRUE
    )
  )

  cooks_sample <- merge(
    cooks_sample,
    metadata[
      ,
      .(
        sample_id,
        group,
        sequencingBatch,
        RIN,
        pmi,
        msex
      )
    ],
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  cooks_sample <- cooks_sample[
    match(metadata$sample_id, sample_id)
  ]
}

fwrite(
  cooks_sample,
  cooks_output,
  sep = "\t",
  na = "NA"
)

saveRDS(dds_main, dds_main_output)
saveRDS(dds_adjusted, dds_adjusted_output)


# ==============================================================================
# 16. MA PLOT
# ==============================================================================

png(
  ma_plot_output,
  width = 1800,
  height = 1500,
  res = 220
)

plotMA(
  res_main,
  alpha = PRIMARY_ALPHA,
  ylim = c(-5, 5),
  main = "ROSMAP — MA plot | AD vs Control"
)

abline(
  h = c(-PRIMARY_LFC, PRIMARY_LFC),
  lty = 2
)

dev.off()


# ==============================================================================
# 17. VOLCANO PLOTS
# ==============================================================================
#
# Para manter o MESMO padrão visual adotado nas fases anteriores:
#
#   eixo X = log2FoldChange
#   eixo Y = -log10(p-value bruto)
#
# A CLASSIFICAÇÃO das cores, entretanto, continua obedecendo ao critério
# principal pré-definido:
#
#   BH-FDR < 0,05 e |log2FC| > 0,58
#
# Portanto:
#   - o p-value bruto é usado somente para a escala visual do eixo Y;
#   - a linha horizontal em p = 0,05 é apenas uma referência visual;
#   - quem é ou não DEG principal continua sendo decidido pelo BH-FDR + log2FC.
#
# Isso evita a compressão visual que ocorre quando se usa -log10(padj) em um
# conjunto cujo menor padj é relativamente moderado.
# ==============================================================================

make_volcano <- function(
  res_dt,
  output_file,
  model_label
) {

  volcano_dt <- copy(
    res_dt
  )

  positive_p <- volcano_dt$pvalue[
    !is.na(
      volcano_dt$pvalue
    ) &
      volcano_dt$pvalue > 0
  ]

  p_floor <- if (
    length(
      positive_p
    ) > 0
  ) {
    min(
      positive_p
    ) / 10
  } else {
    .Machine$double.xmin
  }

  volcano_dt[
    ,
    pvalue_for_plot := fifelse(
      is.na(
        pvalue
      ),
      1,
      fifelse(
        pvalue <= 0,
        p_floor,
        pvalue
      )
    )
  ]

  volcano_dt[
    ,
    minus_log10_pvalue := -log10(
      pvalue_for_plot
    )
  ]

  volcano_dt[
    ,
    volcano_class := factor(
      fifelse(
        DEG &
          direction == "Up",
        "Primary Up",
        fifelse(
          DEG &
            direction == "Down",
          "Primary Down",
          "Not primary DEG"
        )
      ),
      levels = c(
        "Not primary DEG",
        "Primary Down",
        "Primary Up"
      )
    )
  ]

  volcano_labels <- volcano_dt[
    DEG == TRUE &
      !is.na(
        padj
      )
  ][
    order(
      padj,
      -abs(
        log2FoldChange
      )
    )
  ]

  if (
    nrow(
      volcano_labels
    ) > VOLCANO_LABEL_N
  ) {
    volcano_labels <- volcano_labels[
      seq_len(
        VOLCANO_LABEL_N
      )
    ]
  }

  p <- ggplot(
    volcano_dt[
      is.finite(
        log2FoldChange
      ) &
        is.finite(
          minus_log10_pvalue
        )
    ],
    aes(
      x = log2FoldChange,
      y = minus_log10_pvalue,
      color = volcano_class
    )
  ) +
    geom_point(
      alpha = 0.60,
      size = 1.6
    ) +
    geom_vline(
      xintercept = c(
        -PRIMARY_LFC,
        PRIMARY_LFC
      ),
      linetype = 2,
      linewidth = 0.4
    ) +
    geom_hline(
      yintercept = -log10(
        NOMINAL_P
      ),
      linetype = 2,
      linewidth = 0.4
    ) +
    geom_text(
      data = volcano_labels,
      aes(
        label = plot_label
      ),
      check_overlap = TRUE,
      size = 3,
      vjust = -0.5,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Not primary DEG" = "grey70",
        "Primary Down" = "#2C7BB6",
        "Primary Up" = "#D7191C"
      ),
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "ROSMAP — Volcano plot | ",
        model_label
      ),
      subtitle = paste0(
        "Cores: BH-FDR < ",
        PRIMARY_ALPHA,
        " e |log2FC| > ",
        PRIMARY_LFC,
        " | linha horizontal: p bruto = ",
        NOMINAL_P
      ),
      x = "log2FoldChange (AD / Control)",
      y = "-log10(p-value)",
      color = "Classificação"
    ) +
    theme_bw(
      base_size = 13
    )

  ggsave(
    filename = output_file,
    plot = p,
    width = 8.5,
    height = 6.5,
    dpi = 300
  )

  invisible(
    p
  )
}


volcano_main_plot <- make_volcano(
  res_main_dt,
  volcano_plot_output,
  "modelo principal"
)

volcano_adjusted_plot <- make_volcano(
  res_adjusted_dt,
  volcano_adjusted_plot_output,
  "modelo ajustado"
)


# ==============================================================================
# 18. SCATTER PRINCIPAL × AJUSTADO
# ==============================================================================

stability_plot <- ggplot(
  stability,
  aes(
    x = log2FC_main,
    y = log2FC_adjusted,
    color = primary_overlap_class
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  geom_point(
    alpha = 0.45,
    size = 1.5
  ) +
  labs(
    title = "ROSMAP — Estabilidade do efeito AD",
    subtitle = paste0(
      "log2FC principal vs ajustado | Spearman rho = ",
      fmt_num(spearman_all, 3)
    ),
    x = "log2FC — modelo principal (~ group)",
    y = "log2FC — modelo ajustado",
    color = "Critério principal"
  ) +
  theme_bw(base_size = 13)

ggsave(
  stability_plot_output,
  stability_plot,
  width = 7.5,
  height = 6.5,
  dpi = 300
)


# ==============================================================================
# 19. HEATMAP DOS DEGs PRINCIPAIS
# ==============================================================================
#
# Correção metodológica/visual:
#
# A versão anterior completava o heatmap até 50 genes mesmo quando havia menos
# de 50 DEGs principais. No ROSMAP, o modelo principal possui apenas 15 DEGs
# pelo critério pré-definido. Portanto, adicionar 35 genes não-DEG ao mesmo
# heatmap poderia sugerir incorretamente que todos eram DEGs principais.
#
# Nesta versão:
#
#   - usamos SOMENTE os DEGs principais;
#   - usamos no máximo 50 genes;
#   - se houver 15 DEGs, o heatmap terá 15 genes;
#   - a matriz VST continua sendo usada apenas para visualização;
#   - scale = "row" gera z-score por gene;
#   - pheatmap restaura o padrão visual empregado nas fases anteriores;
#   - as colunas recebem anotação de Group e sequencingBatch;
#   - clustering é exploratório e não altera a análise diferencial.
# ==============================================================================

vst_matrix <- readRDS(
  vst_input
)

if (
  !identical(
    colnames(
      vst_matrix
    ),
    metadata$sample_id
  ) ||
    !identical(
      rownames(
        vst_matrix
      ),
      rownames(
        count_matrix
      )
    )
) {
  stop(
    "Matriz VST não corresponde aos counts/metadata atuais."
  )
}


make_primary_heatmap <- function(
  res_dt,
  output_file,
  gene_table_file,
  model_label
) {

  selected <- res_dt[
    DEG == TRUE &
      !is.na(
        padj
      )
  ][
    order(
      padj,
      -abs(
        log2FoldChange
      )
    )
  ]

  if (
    nrow(
      selected
    ) == 0
  ) {
    warning(
      "Nenhum DEG principal disponível para heatmap — ",
      model_label,
      "."
    )
    return(
      invisible(
        NULL
      )
    )
  }

  selected <- selected[
    seq_len(
      min(
        HEATMAP_TOP_N,
        .N
      )
    )
  ]

  selected_ids <- selected$GeneID

  heatmap_info <- selected[
    ,
    .(
      GeneID,
      plot_label,
      log2FoldChange,
      pvalue,
      padj,
      DEG,
      direction
    )
  ]

  fwrite(
    heatmap_info,
    gene_table_file,
    sep = "\t",
    na = "NA"
  )

  mat <- vst_matrix[
    selected_ids,
    ,
    drop = FALSE
  ]

  row_sd <- apply(
    mat,
    1,
    sd
  )

  keep_variable <- is.finite(
    row_sd
  ) &
    row_sd > 0

  mat <- mat[
    keep_variable,
    ,
    drop = FALSE
  ]

  if (
    nrow(
      mat
    ) < 2
  ) {
    warning(
      "Menos de 2 genes variáveis disponíveis para heatmap — ",
      model_label,
      "."
    )
    return(
      invisible(
        NULL
      )
    )
  }

  labels <- selected$plot_label
  names(
    labels
  ) <- selected$GeneID

  rownames(
    mat
  ) <- make.unique(
    labels[
      rownames(
        mat
      )
    ]
  )

  annotation_col <- data.frame(
    Group = factor(
      metadata$group,
      levels = c(
        "Control",
        "AD"
      )
    ),
    sequencingBatch = factor(
      metadata$sequencingBatch
    ),
    row.names = metadata$sample_id,
    check.names = FALSE
  )

  pheatmap::pheatmap(
    mat,
    scale = "row",
    annotation_col = annotation_col,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    border_color = NA,
    fontsize = 10,
    fontsize_row = 9,
    main = paste0(
      "ROSMAP — DEGs principais | ",
      model_label
    ),
    filename = output_file,
    width = 10,
    height = max(
      5.5,
      0.28 *
        nrow(
          mat
        ) +
        3
    )
  )

  invisible(
    selected_ids
  )
}


heatmap_main_ids <- make_primary_heatmap(
  res_main_dt,
  heatmap_plot_output,
  top_genes_output,
  "modelo principal"
)

heatmap_adjusted_ids <- make_primary_heatmap(
  res_adjusted_dt,
  heatmap_adjusted_plot_output,
  top_genes_adjusted_output,
  "modelo ajustado"
)


# ==============================================================================
# 20. DISPERSÕES
# ==============================================================================

png(
  dispersion_main_output,
  width = 1700,
  height = 1400,
  res = 220
)

plotDispEsts(dds_main)
title(main = "ROSMAP — Dispersões DESeq2 | modelo principal")
dev.off()

png(
  dispersion_adjusted_output,
  width = 1700,
  height = 1400,
  res = 220
)

plotDispEsts(dds_adjusted)
title(main = "ROSMAP — Dispersões DESeq2 | modelo ajustado")
dev.off()


# ==============================================================================
# 21. SESSION INFO
# ==============================================================================

sink(session_info_output)
print(sessionInfo())
sink()


# ==============================================================================
# 22. RELATÓRIO
# ==============================================================================

main_summary <- model_summary[
  model == "main_group_only"
]

adjusted_summary <- model_summary[
  model == "adjusted_batch_msex_RIN_PMI"
]

sink(report_output)

cat(
  "============================================================\n"
)
cat(
  "FASE 01c — 05_deseq2_ROSMAP.R\n"
)
cat(
  "EXPRESSÃO DIFERENCIAL — AD vs CONTROL\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat("1. ENTRADA\n")
cat(
  "------------------------------------------------------------\n"
)
cat(
  "GeneIDs testados: ",
  fmt_int(nrow(count_matrix)),
  "\n",
  sep = ""
)
cat(
  "Amostras: ",
  fmt_int(ncol(count_matrix)),
  "\n",
  sep = ""
)
cat(
  "AD: ",
  fmt_int(n_ad),
  "\n",
  sep = ""
)
cat(
  "Control: ",
  fmt_int(n_control),
  "\n",
  sep = ""
)
cat(
  "Referência: Control\n"
)
cat(
  "log2FC > 0 = maior expressão em AD.\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat("2. CRITÉRIOS PRÉ-DEFINIDOS\n")
cat(
  "------------------------------------------------------------\n"
)
cat(
  "Principal: BH-FDR < 0,05 e |log2FC| > 0,58\n"
)
cat(
  "Estrito:   BH-FDR < 0,05 e |log2FC| > 1\n"
)
cat(
  "Nominal:   p < 0,05 e |log2FC| > 1\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat("3. MODELOS\n")
cat(
  "------------------------------------------------------------\n"
)
cat(
  "Principal:\n  ~ group\n\n"
)
cat(
  "Ajustado:\n  ~ sequencingBatch + msex + RIN + pmi + group\n\n"
)
cat(
  "Implementação:\n"
)
cat(
  "  ~ sequencingBatch_factor + msex_factor + RIN_centered +\n"
)
cat(
  "    pmi_centered + group_factor\n\n"
)
cat(
  "Casos completos no ajustado: 144 / 144\n"
)
cat(
  "Design full rank: TRUE (rank ",
  fmt_int(design_rank),
  " / ",
  fmt_int(design_ncol),
  ")\n\n",
  sep = ""
)
cat(
  "sequencingBatch foi tratado como variável categórica com os rótulos\n"
)
cat(
  "originais do metadata, sem interpretação adicional de 9/NYGC2/NYGC3.\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat("4. MODELO PRINCIPAL (~ group)\n")
cat(
  "------------------------------------------------------------\n"
)
cat(
  "Genes com p-value: ",
  fmt_int(main_summary$n_pvalue_nonmissing),
  "\n",
  sep = ""
)
cat(
  "Genes com padj: ",
  fmt_int(main_summary$n_padj_nonmissing),
  "\n",
  sep = ""
)
cat(
  "FDR < 0,05: ",
  fmt_int(main_summary$FDR_lt_0.05),
  "\n\n",
  sep = ""
)
cat(
  "Principal: ",
  fmt_int(main_summary$primary_total),
  " (Up ",
  fmt_int(main_summary$primary_up),
  "; Down ",
  fmt_int(main_summary$primary_down),
  ")\n",
  sep = ""
)
cat(
  "Estrito: ",
  fmt_int(main_summary$strict_total),
  " (Up ",
  fmt_int(main_summary$strict_up),
  "; Down ",
  fmt_int(main_summary$strict_down),
  ")\n",
  sep = ""
)
cat(
  "Nominal: ",
  fmt_int(main_summary$nominal_total),
  " (Up ",
  fmt_int(main_summary$nominal_up),
  "; Down ",
  fmt_int(main_summary$nominal_down),
  ")\n",
  sep = ""
)
cat(
  "Menor p-value: ",
  format(
    main_summary$min_pvalue,
    scientific = TRUE,
    digits = 5
  ),
  "\n",
  sep = ""
)
cat(
  "Menor padj: ",
  format(
    main_summary$min_padj,
    scientific = TRUE,
    digits = 5
  ),
  "\n\n",
  sep = ""
)

cat(
  "------------------------------------------------------------\n"
)
cat("5. MODELO AJUSTADO\n")
cat(
  "------------------------------------------------------------\n"
)
cat(
  "Genes com p-value: ",
  fmt_int(adjusted_summary$n_pvalue_nonmissing),
  "\n",
  sep = ""
)
cat(
  "Genes com padj: ",
  fmt_int(adjusted_summary$n_padj_nonmissing),
  "\n",
  sep = ""
)
cat(
  "FDR < 0,05: ",
  fmt_int(adjusted_summary$FDR_lt_0.05),
  "\n\n",
  sep = ""
)
cat(
  "Principal: ",
  fmt_int(adjusted_summary$primary_total),
  " (Up ",
  fmt_int(adjusted_summary$primary_up),
  "; Down ",
  fmt_int(adjusted_summary$primary_down),
  ")\n",
  sep = ""
)
cat(
  "Estrito: ",
  fmt_int(adjusted_summary$strict_total),
  " (Up ",
  fmt_int(adjusted_summary$strict_up),
  "; Down ",
  fmt_int(adjusted_summary$strict_down),
  ")\n",
  sep = ""
)
cat(
  "Nominal: ",
  fmt_int(adjusted_summary$nominal_total),
  " (Up ",
  fmt_int(adjusted_summary$nominal_up),
  "; Down ",
  fmt_int(adjusted_summary$nominal_down),
  ")\n",
  sep = ""
)
cat(
  "Menor p-value: ",
  format(
    adjusted_summary$min_pvalue,
    scientific = TRUE,
    digits = 5
  ),
  "\n",
  sep = ""
)
cat(
  "Menor padj: ",
  format(
    adjusted_summary$min_padj,
    scientific = TRUE,
    digits = 5
  ),
  "\n\n",
  sep = ""
)

cat(
  "------------------------------------------------------------\n"
)
cat("6. ESTABILIDADE PRINCIPAL × AJUSTADO\n")
cat(
  "------------------------------------------------------------\n"
)
cat(
  "Genes comparados: ",
  fmt_int(nrow(stability)),
  "\n",
  sep = ""
)
cat(
  "Spearman log2FC — todos: ",
  fmt_num(spearman_all, 4),
  "\n",
  sep = ""
)
cat(
  "Pearson log2FC — todos: ",
  fmt_num(pearson_all, 4),
  "\n\n",
  sep = ""
)
cat(
  "Shared primary: ",
  fmt_int(n_shared_primary),
  "\n",
  sep = ""
)
cat(
  "Main only: ",
  fmt_int(n_main_only),
  "\n",
  sep = ""
)
cat(
  "Adjusted only: ",
  fmt_int(n_adjusted_only),
  "\n",
  sep = ""
)
cat(
  "União primary: ",
  fmt_int(primary_union),
  "\n",
  sep = ""
)
cat(
  "Jaccard: ",
  fmt_num(primary_jaccard, 4),
  "\n",
  sep = ""
)
cat(
  "Shared mesma direção: ",
  fmt_int(n_shared_same_direction),
  "\n",
  sep = ""
)
cat(
  "Shared direção oposta: ",
  fmt_int(n_shared_opposite_direction),
  "\n",
  sep = ""
)
cat(
  "Troca de sinal — todos genes: ",
  fmt_int(n_sign_switch_all),
  "\n\n",
  sep = ""
)
cat(
  "Spearman — união primary: ",
  fmt_num(spearman_primary_union, 4),
  "\n",
  sep = ""
)
cat(
  "Pearson — união primary: ",
  fmt_num(pearson_primary_union, 4),
  "\n\n",
  sep = ""
)
cat(
  "Estas métricas são descritivas; não há cutoff pós-hoc automático para\n"
)
cat(
  "declarar o ROSMAP adequado ou inadequado.\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat("7. COOK'S DISTANCE\n")
cat(
  "------------------------------------------------------------\n"
)
cat(
  "Cook's distance foi auditada por amostra.\n"
)
cat(
  "Nenhuma amostra foi removida manualmente neste script.\n"
)
cat(
  "p-value NA no principal: ",
  fmt_int(main_summary$n_pvalue_missing),
  "\n",
  sep = ""
)
cat(
  "padj NA no principal: ",
  fmt_int(main_summary$n_padj_missing),
  "\n\n",
  sep = ""
)

cat(
  "------------------------------------------------------------\n"
)
cat("8. FIGURAS\n")
cat(
  "------------------------------------------------------------\n"
)
cat(
  "MA: ",
  ma_plot_output,
  "\n",
  sep = ""
)
cat(
  "Volcano principal: ",
  volcano_plot_output,
  "\n",
  sep = ""
)
cat(
  "Volcano ajustado: ",
  volcano_adjusted_plot_output,
  "\n",
  sep = ""
)
cat(
  "Heatmap principal: ",
  heatmap_plot_output,
  "\n",
  sep = ""
)
cat(
  "Heatmap ajustado: ",
  heatmap_adjusted_plot_output,
  "\n",
  sep = ""
)
cat(
  "Principal × ajustado: ",
  stability_plot_output,
  "\n",
  sep = ""
)
cat(
  "Dispersões principal: ",
  dispersion_main_output,
  "\n",
  sep = ""
)
cat(
  "Dispersões ajustado: ",
  dispersion_adjusted_output,
  "\n\n",
  sep = ""
)

cat(
  "------------------------------------------------------------\n"
)
cat("9. CONCLUSÃO\n")
cat(
  "------------------------------------------------------------\n"
)
cat(
  "- Contraste: AD vs Control, com Control como referência.\n"
)
cat(
  "- Modelo principal: ~ group.\n"
)
cat(
  "- Modelo ajustado: sequencingBatch + msex + RIN + PMI + group.\n"
)
cat(
  "- Os mesmos 144 indivíduos foram usados nos dois modelos.\n"
)
cat(
  "- Os critérios de DEG foram definidos antes da análise do ROSMAP.\n"
)
cat(
  "- A estabilidade foi avaliada por correlação, sobreposição e direção.\n"
)
cat(
  "- Os heatmaps incluem somente DEGs que satisfazem o critério principal;\n"
)
cat(
  "  genes não significativos não são usados para completar artificialmente 50 linhas.\n"
)
cat(
  "- O número absoluto de DEGs não deve ser usado isoladamente para escolher\n"
)
cat(
  "  o dataset principal de AD.\n"
)
cat(
  "- Resultados prontos para enriquecimento e comparação com GSE125583.\n\n"
)

cat(
  "============================================================\n"
)
cat(
  "05_deseq2_ROSMAP.R CONCLUÍDO\n"
)
cat(
  "============================================================\n"
)

sink()


# ==============================================================================
# 23. SESSION INFO
# ==============================================================================

sink(session_info_output)
print(sessionInfo())
sink()


# ==============================================================================
# 24. RESUMO NO TERMINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "DESeq2 ROSMAP CONCLUÍDO\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Principal: ",
  fmt_int(main_summary$primary_total),
  " DEGs (Up ",
  fmt_int(main_summary$primary_up),
  "; Down ",
  fmt_int(main_summary$primary_down),
  ")\n",
  sep = ""
)

cat(
  "Ajustado: ",
  fmt_int(adjusted_summary$primary_total),
  " DEGs (Up ",
  fmt_int(adjusted_summary$primary_up),
  "; Down ",
  fmt_int(adjusted_summary$primary_down),
  ")\n\n",
  sep = ""
)

cat(
  "Spearman log2FC: ",
  fmt_num(spearman_all, 4),
  "\n",
  sep = ""
)

cat(
  "Shared primary: ",
  fmt_int(n_shared_primary),
  "\n",
  sep = ""
)

cat(
  "Main only: ",
  fmt_int(n_main_only),
  "\n",
  sep = ""
)

cat(
  "Adjusted only: ",
  fmt_int(n_adjusted_only),
  "\n",
  sep = ""
)

cat(
  "Shared mesma direção: ",
  fmt_int(n_shared_same_direction),
  "\n",
  sep = ""
)

cat(
  "Shared direção oposta: ",
  fmt_int(n_shared_opposite_direction),
  "\n\n",
  sep = ""
)

cat(
  "Relatório:\n  ",
  report_output,
  "\n\n",
  sep = ""
)

cat(
  "Próxima etapa após revisar os resultados:\n"
)
cat(
  "  06_pathway_enrichment_ROSMAP.R\n"
)
