#!/usr/bin/env Rscript

# ==============================================================================
# 07_compare_AD_datasets.R
#
# FASE 01b — COMPARAÇÃO FORMAL DOS DOIS DATASETS DE ALZHEIMER
#
# DATASET ORIGINAL:
#   GSE53697
#
# DATASET CANDIDATO:
#   GSE125583
#
# OBJETIVO GERAL
# ------------------------------------------------------------------------------
# Comparar formalmente o dataset de AD utilizado na análise original
# (GSE53697) com o novo candidato (GSE125583), documentando:
#
#   1. tamanho amostral;
#   2. universo gênico bruto;
#   3. universo testável após o mesmo princípio de pré-filtragem;
#   4. robustez da expressão diferencial;
#   5. concordância dos efeitos entre os dois estudos;
#   6. sobreposição das listas de DEGs;
#   7. compatibilidade do GSE125583 com os universos testáveis de HD e PD.
#
# A finalidade NÃO é escolher o novo dataset simplesmente porque ele produz
# "mais DEGs". A decisão deve combinar:
#
#   - adequação técnica;
#   - tamanho amostral;
#   - quantidade de genes testáveis;
#   - robustez após correção por múltiplos testes;
#   - estabilidade do efeito após ajuste por covariáveis;
#   - compatibilidade gênica com HD/PD;
#   - coerência funcional;
#   - limitações de tecido/região cerebral e desenho de estudo.
#
#
# CONTEXTO METODOLÓGICO IMPORTANTE
# ------------------------------------------------------------------------------
# GSE53697 e GSE125583 NÃO representam uma replicação perfeita um do outro:
#
#   GSE53697  -> dorsolateral prefrontal cortex (DLPFC)
#   GSE125583 -> fusiform gyrus
#
# Portanto, baixa concordância gene-a-gene NÃO deve ser automaticamente
# interpretada como falha de um dataset. Diferenças podem decorrer de:
#
#   - região cerebral;
#   - composição celular;
#   - tamanho amostral;
#   - características clínicas;
#   - processamento original;
#   - heterogeneidade biológica.
#
# Esta etapa avalia principalmente a USABILIDADE do novo dataset para substituir
# o GSE53697 na etapa transdiagnóstica/deconvolução.
#
#
# CRITÉRIOS DE DEG COMPARADOS
# ------------------------------------------------------------------------------
# Critério principal:
#   padj < 0,05 e |log2FC| > 0,58
#
# Sensibilidade:
#   padj < 0,05 e |log2FC| > 1
#
# Critério nominal comparativo ("article-like"):
#   pvalue < 0,05 e |log2FC| > 1
#
# O critério nominal é apenas comparativo e NÃO substitui FDR.
#
#
# ENTRADAS — ANÁLISE ORIGINAL
# ------------------------------------------------------------------------------
# results/counts/counts_AD.tsv
# results/qc/counts_AD_DE_filtered.tsv
# results/metadata/metadata_AD.tsv
# results/deseq2/AD_deseq2_results.tsv
#
# results/qc/counts_HD_DE_filtered.tsv
# results/qc/counts_PD_DE_filtered.tsv
# results/deseq2/HD_deseq2_results.tsv
# results/deseq2/PD_deseq2_results.tsv
#
#
# ENTRADAS — GSE125583
# ------------------------------------------------------------------------------
# results/candidate_GSE125583/counts/counts_GSE125583.tsv
# results/candidate_GSE125583/qc/counts_GSE125583_DE_filtered.tsv
# results/candidate_GSE125583/metadata/metadata_GSE125583.tsv
# results/candidate_GSE125583/deseq2/GSE125583_deseq2_results.tsv
#
#
# SAÍDAS
# ------------------------------------------------------------------------------
# results/candidate_GSE125583/comparison_AD/
# ├── AD_dataset_characteristics.tsv
# ├── AD_DEG_comparison.tsv
# ├── AD_testable_universe_comparison.tsv
# ├── AD_old_vs_new_common_genes.tsv
# ├── AD_old_vs_new_DEG_overlap.tsv
# ├── candidate_AD_HD_PD_testable_universe.tsv
# └── 07_compare_AD_datasets.txt
#
# results/candidate_GSE125583/plots/comparison_AD/
# ├── AD_dataset_sample_size.png
# ├── AD_dataset_testable_genes.png
# ├── AD_dataset_DEG_counts.png
# └── AD_old_vs_new_log2FC.png
#
# ==============================================================================


# ==============================================================================
# 1. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Pacotes necessários não instalados: ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

# ------------------------------------------------------------------------------
# GSE53697 — AD original
# ------------------------------------------------------------------------------

old_raw_counts_file <- file.path(
  "results",
  "counts",
  "counts_AD.tsv"
)

old_filtered_counts_file <- file.path(
  "results",
  "qc",
  "counts_AD_DE_filtered.tsv"
)

old_metadata_file <- file.path(
  "results",
  "metadata",
  "metadata_AD.tsv"
)

old_deseq_file <- file.path(
  "results",
  "deseq2",
  "AD_deseq2_results.tsv"
)


# ------------------------------------------------------------------------------
# GSE125583 — novo AD candidato
# ------------------------------------------------------------------------------

new_raw_counts_file <- file.path(
  "results",
  "candidate_GSE125583",
  "counts",
  "counts_GSE125583.tsv"
)

new_filtered_counts_file <- file.path(
  "results",
  "candidate_GSE125583",
  "qc",
  "counts_GSE125583_DE_filtered.tsv"
)

new_metadata_file <- file.path(
  "results",
  "candidate_GSE125583",
  "metadata",
  "metadata_GSE125583.tsv"
)

new_deseq_file <- file.path(
  "results",
  "candidate_GSE125583",
  "deseq2",
  "GSE125583_deseq2_results.tsv"
)


# ------------------------------------------------------------------------------
# HD / PD — universos testáveis já definidos na Fase 01 original
# ------------------------------------------------------------------------------

hd_filtered_counts_file <- file.path(
  "results",
  "qc",
  "counts_HD_DE_filtered.tsv"
)

pd_filtered_counts_file <- file.path(
  "results",
  "qc",
  "counts_PD_DE_filtered.tsv"
)


# ------------------------------------------------------------------------------
# Diretórios de saída
# ------------------------------------------------------------------------------

output_dir <- file.path(
  "results",
  "candidate_GSE125583",
  "comparison_AD"
)

plot_dir <- file.path(
  "results",
  "candidate_GSE125583",
  "plots",
  "comparison_AD"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

characteristics_file <- file.path(
  output_dir,
  "AD_dataset_characteristics.tsv"
)

deg_comparison_file <- file.path(
  output_dir,
  "AD_DEG_comparison.tsv"
)

testable_comparison_file <- file.path(
  output_dir,
  "AD_testable_universe_comparison.tsv"
)

common_genes_file <- file.path(
  output_dir,
  "AD_old_vs_new_common_genes.tsv"
)

deg_overlap_file <- file.path(
  output_dir,
  "AD_old_vs_new_DEG_overlap.tsv"
)

cross_disease_universe_file <- file.path(
  output_dir,
  "candidate_AD_HD_PD_testable_universe.tsv"
)

audit_file <- file.path(
  output_dir,
  "07_compare_AD_datasets.txt"
)


# ==============================================================================
# 3. CONSTANTES ANALÍTICAS
# ==============================================================================

ALPHA <- 0.05
PRIMARY_LFC <- 0.58
SENSITIVITY_LFC <- 1.00

DATASET_OLD <- "GSE53697"
DATASET_NEW <- "GSE125583"

TISSUE_OLD <- "dorsolateral prefrontal cortex"
TISSUE_NEW <- "fusiform gyrus"

# Paleta consistente para os dois datasets de AD em todos os gráficos.
DATASET_COLORS <- c(
  "GSE53697" = "#00A6A6",
  "GSE125583" = "#E67E22"
)

# Paleta para o scatter de concordância dos efeitos.
# O critério usado aqui é o nominal comparativo (article-like):
# p < 0,05 e |log2FC| > 1.
EFFECT_CLASS_COLORS <- c(
  "Nenhum" = "grey80",
  "GSE53697" = "#00A6A6",
  "GSE125583" = "#E67E22",
  "Ambos" = "#8E44AD"
)


# ==============================================================================
# 4. FUNÇÕES AUXILIARES
# ==============================================================================

validate_file <- function(path) {

  if (!file.exists(path)) {
    stop(
      "Arquivo não encontrado: ",
      path
    )
  }

  info <- file.info(path)

  if (is.na(info$size) || info$size == 0) {
    stop(
      "Arquivo vazio ou tamanho indisponível: ",
      path
    )
  }

  invisible(TRUE)
}


format_int <- function(x) {

  format(
    x,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


format_num <- function(
  x,
  digits = 3
) {

  format(
    round(
      x,
      digits
    ),
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


format_p <- function(x) {

  if (is.na(x)) {
    return("NA")
  }

  format(
    x,
    scientific = TRUE,
    digits = 4
  )
}


# ------------------------------------------------------------------------------
# Ler GeneIDs de matriz de counts
# ------------------------------------------------------------------------------

read_gene_ids_from_counts <- function(
  path,
  dataset_name
) {

  dt <- fread(
    path,
    select = "GeneID",
    showProgress = FALSE
  )

  if (!"GeneID" %in% names(dt)) {
    stop(
      dataset_name,
      ": coluna GeneID não encontrada em ",
      path
    )
  }

  ids <- as.character(
    dt$GeneID
  )

  if (
    any(
      is.na(ids) |
        trimws(ids) == ""
    )
  ) {
    stop(
      dataset_name,
      ": GeneIDs ausentes."
    )
  }

  if (anyDuplicated(ids)) {
    stop(
      dataset_name,
      ": GeneIDs duplicados."
    )
  }

  ids
}


# ------------------------------------------------------------------------------
# Padronizar resultados DESeq2 antigos/novos
# ------------------------------------------------------------------------------
#
# A Fase 01 original usa "gene_id".
# A Fase 01b usa "GeneID".
#
# Esta função converte os dois formatos para um padrão único.

standardize_deseq <- function(
  dt,
  dataset_name
) {

  id_col <- if ("GeneID" %in% names(dt)) {

    "GeneID"

  } else if ("gene_id" %in% names(dt)) {

    "gene_id"

  } else {

    stop(
      dataset_name,
      ": não foi encontrada coluna GeneID/gene_id."
    )
  }

  required_stats <- c(
    "log2FoldChange",
    "pvalue",
    "padj"
  )

  missing_stats <- setdiff(
    required_stats,
    names(dt)
  )

  if (length(missing_stats) > 0) {
    stop(
      dataset_name,
      ": colunas DESeq2 ausentes: ",
      paste(
        missing_stats,
        collapse = ", "
      )
    )
  }

  out <- data.table(
    GeneID = as.character(
      dt[[id_col]]
    ),
    log2FoldChange = as.numeric(
      dt$log2FoldChange
    ),
    pvalue = as.numeric(
      dt$pvalue
    ),
    padj = as.numeric(
      dt$padj
    )
  )

  if ("stat" %in% names(dt)) {
    out[
      ,
      stat := as.numeric(
        dt$stat
      )
    ]
  } else {
    out[
      ,
      stat := NA_real_
    ]
  }

  if ("Symbol" %in% names(dt)) {
    out[
      ,
      Symbol := as.character(
        dt$Symbol
      )
    ]
  } else {
    out[
      ,
      Symbol := NA_character_
    ]
  }

  if (anyDuplicated(out$GeneID)) {
    stop(
      dataset_name,
      ": GeneIDs duplicados nos resultados padronizados."
    )
  }

  # Recalcular as classificações explicitamente, garantindo que os dois
  # datasets sejam comparados pelos MESMOS thresholds.
  out[
    ,
    fdr_significant := (
      !is.na(padj) &
        padj < ALPHA
    )
  ]

  out[
    ,
    primary := (
      !is.na(padj) &
        padj < ALPHA &
        !is.na(log2FoldChange) &
        abs(log2FoldChange) > PRIMARY_LFC
    )
  ]

  out[
    ,
    sensitivity := (
      !is.na(padj) &
        padj < ALPHA &
        !is.na(log2FoldChange) &
        abs(log2FoldChange) > SENSITIVITY_LFC
    )
  ]

  out[
    ,
    article_like := (
      !is.na(pvalue) &
        pvalue < ALPHA &
        !is.na(log2FoldChange) &
        abs(log2FoldChange) > SENSITIVITY_LFC
    )
  ]

  out[
    ,
    primary_direction := fcase(
      primary &
        log2FoldChange > 0,
      "Up",

      primary &
        log2FoldChange < 0,
      "Down",

      default = "Not_DEG"
    )
  ]

  out[
    ,
    article_direction := fcase(
      article_like &
        log2FoldChange > 0,
      "Up",

      article_like &
        log2FoldChange < 0,
      "Down",

      default = "Not_DEG"
    )
  ]

  out
}


# ------------------------------------------------------------------------------
# Resumir resultados de DE
# ------------------------------------------------------------------------------

summarize_deseq <- function(
  dt,
  dataset_name
) {

  data.table(
    dataset = dataset_name,

    n_tested = nrow(dt),

    n_pvalue_nonNA = sum(
      !is.na(
        dt$pvalue
      )
    ),

    n_padj_nonNA = sum(
      !is.na(
        dt$padj
      )
    ),

    n_FDR_lt_005 = sum(
      dt$fdr_significant,
      na.rm = TRUE
    ),

    n_primary = sum(
      dt$primary,
      na.rm = TRUE
    ),

    n_primary_up = sum(
      dt$primary_direction == "Up",
      na.rm = TRUE
    ),

    n_primary_down = sum(
      dt$primary_direction == "Down",
      na.rm = TRUE
    ),

    n_sensitivity = sum(
      dt$sensitivity,
      na.rm = TRUE
    ),

    n_article_like = sum(
      dt$article_like,
      na.rm = TRUE
    ),

    min_pvalue = if (
      all(
        is.na(
          dt$pvalue
        )
      )
    ) {
      NA_real_
    } else {
      min(
        dt$pvalue,
        na.rm = TRUE
      )
    },

    min_padj = if (
      all(
        is.na(
          dt$padj
        )
      )
    ) {
      NA_real_
    } else {
      min(
        dt$padj,
        na.rm = TRUE
      )
    }
  )
}


# ------------------------------------------------------------------------------
# Tema gráfico
# ------------------------------------------------------------------------------

comparison_theme <- function() {

  theme_bw(
    base_size = 12
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 15
      ),
      plot.subtitle = element_text(
        size = 11
      ),
      axis.title = element_text(
        face = "bold"
      ),
      legend.title = element_text(
        face = "bold"
      )
    )
}


# ==============================================================================
# 5. VALIDAR TODAS AS ENTRADAS
# ==============================================================================

input_files <- c(
  old_raw_counts_file,
  old_filtered_counts_file,
  old_metadata_file,
  old_deseq_file,
  new_raw_counts_file,
  new_filtered_counts_file,
  new_metadata_file,
  new_deseq_file,
  hd_filtered_counts_file,
  pd_filtered_counts_file
)

invisible(
  lapply(
    input_files,
    validate_file
  )
)


# ==============================================================================
# 6. LER METADATA DOS DOIS DATASETS DE AD
# ==============================================================================

metadata_old <- fread(
  old_metadata_file,
  na.strings = c(
    "NA",
    ""
  ),
  showProgress = FALSE
)

metadata_new <- fread(
  new_metadata_file,
  na.strings = c(
    "NA",
    ""
  ),
  showProgress = FALSE
)

if (!"group" %in% names(metadata_old)) {
  stop(
    "Metadata original sem coluna group."
  )
}

if (!"group" %in% names(metadata_new)) {
  stop(
    "Metadata GSE125583 sem coluna group."
  )
}

n_old_samples <- nrow(
  metadata_old
)

n_new_samples <- nrow(
  metadata_new
)

n_old_ad <- metadata_old[
  group == "AD",
  .N
]

n_old_control <- metadata_old[
  group == "Control",
  .N
]

n_new_ad <- metadata_new[
  group == "AD",
  .N
]

n_new_control <- metadata_new[
  group == "Control",
  .N
]


# ==============================================================================
# 7. LER UNIVERSOS GÊNICOS BRUTOS E TESTÁVEIS
# ==============================================================================

old_raw_genes <- read_gene_ids_from_counts(
  old_raw_counts_file,
  DATASET_OLD
)

new_raw_genes <- read_gene_ids_from_counts(
  new_raw_counts_file,
  DATASET_NEW
)

old_testable_genes <- read_gene_ids_from_counts(
  old_filtered_counts_file,
  paste0(
    DATASET_OLD,
    " filtered"
  )
)

new_testable_genes <- read_gene_ids_from_counts(
  new_filtered_counts_file,
  paste0(
    DATASET_NEW,
    " filtered"
  )
)

hd_testable_genes <- read_gene_ids_from_counts(
  hd_filtered_counts_file,
  "HD"
)

pd_testable_genes <- read_gene_ids_from_counts(
  pd_filtered_counts_file,
  "PD"
)


# ==============================================================================
# 8. COMPARAR UNIVERSOS AD ANTIGO vs AD NOVO
# ==============================================================================

raw_intersection <- intersect(
  old_raw_genes,
  new_raw_genes
)

raw_union <- union(
  old_raw_genes,
  new_raw_genes
)

testable_intersection <- intersect(
  old_testable_genes,
  new_testable_genes
)

testable_union <- union(
  old_testable_genes,
  new_testable_genes
)

old_testable_only <- setdiff(
  old_testable_genes,
  new_testable_genes
)

new_testable_only <- setdiff(
  new_testable_genes,
  old_testable_genes
)

testable_jaccard <- (
  length(
    testable_intersection
  ) /
    length(
      testable_union
    )
)


testable_comparison <- data.table(
  metric = c(
    "raw_GSE53697",
    "raw_GSE125583",
    "raw_intersection",
    "raw_union",
    "testable_GSE53697",
    "testable_GSE125583",
    "testable_intersection",
    "testable_union",
    "testable_only_GSE53697",
    "testable_only_GSE125583"
  ),

  n_genes = c(
    length(
      old_raw_genes
    ),
    length(
      new_raw_genes
    ),
    length(
      raw_intersection
    ),
    length(
      raw_union
    ),
    length(
      old_testable_genes
    ),
    length(
      new_testable_genes
    ),
    length(
      testable_intersection
    ),
    length(
      testable_union
    ),
    length(
      old_testable_only
    ),
    length(
      new_testable_only
    )
  )
)

fwrite(
  testable_comparison,
  testable_comparison_file,
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 9. LER E PADRONIZAR RESULTADOS DESeq2
# ==============================================================================

old_deseq_raw <- fread(
  old_deseq_file,
  na.strings = c(
    "NA",
    ""
  ),
  showProgress = FALSE
)

new_deseq_raw <- fread(
  new_deseq_file,
  na.strings = c(
    "NA",
    ""
  ),
  showProgress = FALSE
)

old_deseq <- standardize_deseq(
  old_deseq_raw,
  DATASET_OLD
)

new_deseq <- standardize_deseq(
  new_deseq_raw,
  DATASET_NEW
)


# ==============================================================================
# 10. VALIDAR QUE RESULTADOS DESeq2 CORRESPONDEM AOS UNIVERSOS TESTÁVEIS
# ==============================================================================

if (!setequal(
  old_deseq$GeneID,
  old_testable_genes
)) {
  stop(
    "GSE53697: resultados DESeq2 não correspondem exatamente aos genes ",
    "da matriz filtrada."
  )
}

if (!setequal(
  new_deseq$GeneID,
  new_testable_genes
)) {
  stop(
    "GSE125583: resultados DESeq2 não correspondem exatamente aos genes ",
    "da matriz filtrada."
  )
}


# ==============================================================================
# 11. RESUMO DOS RESULTADOS DE EXPRESSÃO DIFERENCIAL
# ==============================================================================

deg_comparison <- rbind(
  summarize_deseq(
    old_deseq,
    DATASET_OLD
  ),
  summarize_deseq(
    new_deseq,
    DATASET_NEW
  ),
  fill = TRUE
)

fwrite(
  deg_comparison,
  deg_comparison_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 12. CONSTRUIR TABELA DOS GENES TESTÁVEIS EM AMBOS OS DATASETS
# ==============================================================================

old_common <- old_deseq[
  GeneID %in%
    testable_intersection,
  .(
    GeneID,
    old_log2FC = log2FoldChange,
    old_pvalue = pvalue,
    old_padj = padj,
    old_primary = primary,
    old_article_like = article_like,
    old_primary_direction = primary_direction,
    old_article_direction = article_direction
  )
]

new_common <- new_deseq[
  GeneID %in%
    testable_intersection,
  .(
    GeneID,
    Symbol,
    new_log2FC = log2FoldChange,
    new_pvalue = pvalue,
    new_padj = padj,
    new_primary = primary,
    new_article_like = article_like,
    new_primary_direction = primary_direction,
    new_article_direction = article_direction
  )
]

common_comparison <- merge(
  old_common,
  new_common,
  by = "GeneID",
  all = FALSE,
  sort = FALSE
)

if (nrow(
  common_comparison
) != length(
  testable_intersection
)) {
  stop(
    "Número inesperado de genes na comparação AD antigo vs AD novo."
  )
}

common_comparison[
  ,
  same_effect_direction := (
    !is.na(
      old_log2FC
    ) &
      !is.na(
        new_log2FC
      ) &
      sign(
        old_log2FC
      ) ==
      sign(
        new_log2FC
      )
  )
]

fwrite(
  common_comparison,
  common_genes_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 13. CORRELAÇÃO DE log2FC ENTRE OS DOIS DATASETS DE AD
# ==============================================================================

valid_effects <- common_comparison[
  !is.na(
    old_log2FC
  ) &
    !is.na(
      new_log2FC
    ) &
    is.finite(
      old_log2FC
    ) &
    is.finite(
      new_log2FC
    )
]

rho_old_new <- if (
  nrow(
    valid_effects
  ) >= 3
) {

  suppressWarnings(
    cor(
      valid_effects$old_log2FC,
      valid_effects$new_log2FC,
      method = "spearman"
    )
  )

} else {

  NA_real_
}


# ==============================================================================
# 14. SOBREPOSIÇÃO DOS DEGs
# ==============================================================================

old_primary_ids <- old_deseq[
  primary == TRUE,
  GeneID
]

new_primary_ids <- new_deseq[
  primary == TRUE,
  GeneID
]

old_article_ids <- old_deseq[
  article_like == TRUE,
  GeneID
]

new_article_ids <- new_deseq[
  article_like == TRUE,
  GeneID
]


primary_overlap <- intersect(
  old_primary_ids,
  new_primary_ids
)

article_overlap <- intersect(
  old_article_ids,
  new_article_ids
)


article_overlap_dt <- common_comparison[
  GeneID %in%
    article_overlap
]

n_article_same_direction <- sum(
  article_overlap_dt$old_article_direction ==
    article_overlap_dt$new_article_direction,
  na.rm = TRUE
)

n_article_opposite_direction <- (
  nrow(
    article_overlap_dt
  ) -
    n_article_same_direction
)


deg_overlap <- data.table(
  criterion = c(
    "primary_FDR_0.05_absLFC_0.58",
    "article_like_p_0.05_absLFC_1"
  ),

  GSE53697 = c(
    length(
      old_primary_ids
    ),
    length(
      old_article_ids
    )
  ),

  GSE125583 = c(
    length(
      new_primary_ids
    ),
    length(
      new_article_ids
    )
  ),

  overlap = c(
    length(
      primary_overlap
    ),
    length(
      article_overlap
    )
  )
)

fwrite(
  deg_overlap,
  deg_overlap_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 15. COMPATIBILIDADE DO NOVO AD COM HD E PD
# ==============================================================================

candidate_hd_intersection <- intersect(
  new_testable_genes,
  hd_testable_genes
)

candidate_pd_intersection <- intersect(
  new_testable_genes,
  pd_testable_genes
)

hd_pd_intersection <- intersect(
  hd_testable_genes,
  pd_testable_genes
)

candidate_hd_pd_intersection <- Reduce(
  intersect,
  list(
    new_testable_genes,
    hd_testable_genes,
    pd_testable_genes
  )
)

old_hd_pd_intersection <- Reduce(
  intersect,
  list(
    old_testable_genes,
    hd_testable_genes,
    pd_testable_genes
  )
)


cross_disease_universe <- data.table(
  metric = c(
    "AD_new_testable",
    "HD_testable",
    "PD_testable",
    "AD_new_intersection_HD",
    "AD_new_intersection_PD",
    "HD_intersection_PD",
    "AD_new_intersection_HD_intersection_PD",
    "AD_old_intersection_HD_intersection_PD"
  ),

  n_genes = c(
    length(
      new_testable_genes
    ),
    length(
      hd_testable_genes
    ),
    length(
      pd_testable_genes
    ),
    length(
      candidate_hd_intersection
    ),
    length(
      candidate_pd_intersection
    ),
    length(
      hd_pd_intersection
    ),
    length(
      candidate_hd_pd_intersection
    ),
    length(
      old_hd_pd_intersection
    )
  )
)

fwrite(
  cross_disease_universe,
  cross_disease_universe_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 16. CARACTERÍSTICAS GERAIS DOS DOIS DATASETS
# ==============================================================================

characteristics <- data.table(
  dataset = c(
    DATASET_OLD,
    DATASET_NEW
  ),

  role = c(
    "original_AD",
    "candidate_AD"
  ),

  tissue = c(
    TISSUE_OLD,
    TISSUE_NEW
  ),

  n_samples = c(
    n_old_samples,
    n_new_samples
  ),

  n_AD = c(
    n_old_ad,
    n_new_ad
  ),

  n_Control = c(
    n_old_control,
    n_new_control
  ),

  n_raw_genes = c(
    length(
      old_raw_genes
    ),
    length(
      new_raw_genes
    )
  ),

  n_testable_genes = c(
    length(
      old_testable_genes
    ),
    length(
      new_testable_genes
    )
  )
)

fwrite(
  characteristics,
  characteristics_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 17. GRÁFICO — TAMANHO AMOSTRAL
# ==============================================================================

plot_samples_dt <- melt(
  characteristics[
    ,
    .(
      dataset,
      AD = n_AD,
      Control = n_Control
    )
  ],
  id.vars = "dataset",
  variable.name = "group",
  value.name = "n"
)

p_samples <- ggplot(
  plot_samples_dt,
  aes(
    x = dataset,
    y = n,
    fill = group
  )
) +
  geom_col(
    position = "stack",
    width = 0.65
  ) +
  geom_text(
    aes(
      label = n
    ),
    position = position_stack(
      vjust = 0.5
    ),
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "#7A7A7A",
      "AD" = "#E67E22"
    )
  ) +
  labs(
    title = "Tamanho amostral dos datasets de AD",
    subtitle = "GSE53697 original vs GSE125583 candidato",
    x = NULL,
    y = "Número de amostras",
    fill = "Grupo"
  ) +
  comparison_theme()

ggsave(
  filename = file.path(
    plot_dir,
    "AD_dataset_sample_size.png"
  ),
  plot = p_samples,
  width = 7,
  height = 5,
  dpi = 300
)


# ==============================================================================
# 18. GRÁFICO — GENES TESTÁVEIS
# ==============================================================================
#
# Mantemos uma convenção visual consistente:
#
#   GSE53697  -> turquesa
#   GSE125583 -> laranja
#
# Isso permite que o leitor reconheça imediatamente cada dataset nos gráficos
# subsequentes.

p_testable <- ggplot(
  characteristics,
  aes(
    x = dataset,
    y = n_testable_genes,
    fill = dataset
  )
) +
  geom_col(
    width = 0.65
  ) +
  geom_text(
    aes(
      label = format_int(
        n_testable_genes
      )
    ),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(
    values = DATASET_COLORS,
    guide = "none"
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      big.mark = ".",
      decimal.mark = ","
    ),
    expand = expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  labs(
    title = "Universo testável dos datasets de AD",
    subtitle = "Genes retidos após o pré-filtro para expressão diferencial",
    x = NULL,
    y = "Genes testáveis"
  ) +
  comparison_theme()

ggsave(
  filename = file.path(
    plot_dir,
    "AD_dataset_testable_genes.png"
  ),
  plot = p_testable,
  width = 7,
  height = 5,
  dpi = 300
)


# ==============================================================================
# 19. GRÁFICO — CONTAGEM DE DEGs
# ==============================================================================
#
# Os três critérios são mostrados lado a lado:
#
#   Primary:
#     FDR < 0,05 e |log2FC| > 0,58
#
#   Sensitivity:
#     FDR < 0,05 e |log2FC| > 1
#
#   Article_like:
#     p < 0,05 e |log2FC| > 1
#
# A mesma paleta é usada para os dois datasets em todos os gráficos.

deg_plot_dt <- melt(
  deg_comparison[
    ,
    .(
      dataset,
      Primary = n_primary,
      Sensitivity = n_sensitivity,
      Article_like = n_article_like
    )
  ],
  id.vars = "dataset",
  variable.name = "criterion",
  value.name = "n_genes"
)

deg_plot_dt[
  ,
  criterion := factor(
    criterion,
    levels = c(
      "Primary",
      "Sensitivity",
      "Article_like"
    )
  )
]

p_deg <- ggplot(
  deg_plot_dt,
  aes(
    x = criterion,
    y = n_genes,
    fill = dataset
  )
) +
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7
  ) +
  geom_text(
    aes(
      label = format_int(
        n_genes
      )
    ),
    position = position_dodge(
      width = 0.8
    ),
    vjust = -0.35,
    size = 3.6
  ) +
  scale_fill_manual(
    values = DATASET_COLORS
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      big.mark = ".",
      decimal.mark = ","
    ),
    expand = expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  labs(
    title = "Expressão diferencial nos dois datasets de AD",
    subtitle = "Mesmos thresholds reaplicados aos dois resultados DESeq2",
    x = NULL,
    y = "Número de genes",
    fill = "Dataset"
  ) +
  comparison_theme()

ggsave(
  filename = file.path(
    plot_dir,
    "AD_dataset_DEG_counts.png"
  ),
  plot = p_deg,
  width = 8,
  height = 5.5,
  dpi = 300
)


# ==============================================================================
# 20. GRÁFICO — CORRELAÇÃO DOS EFEITOS
# ==============================================================================
#
# Cada ponto representa um gene testável em ambos os datasets.
#
# Para evitar um scatter visualmente dominado por milhares de pontos cinza,
# classificamos os genes segundo o critério nominal comparativo:
#
#   p < 0,05 e |log2FC| > 1
#
# Classes:
#
#   Nenhum     -> não satisfaz o critério em nenhum dataset
#   GSE53697   -> satisfaz apenas no dataset original
#   GSE125583  -> satisfaz apenas no dataset candidato
#   Ambos      -> satisfaz nos dois datasets
#
# Essa escolha é deliberada:
#
#   - o GSE53697 não possui DEGs pelo critério principal;
#   - o critério nominal comparativo permite visualizar os genes com maior
#     evidência/effect size em ambos os estudos;
#   - ele NÃO substitui o critério principal baseado em FDR.
#
# A linha diagonal y = x é apenas uma referência visual de perfeita concordância,
# e NÃO representa uma regressão.
# ==============================================================================

valid_effects[
  ,
  comparison_class := fcase(
    old_article_like &
      new_article_like,
    "Ambos",

    old_article_like,
    "GSE53697",

    new_article_like,
    "GSE125583",

    default = "Nenhum"
  )
]

valid_effects[
  ,
  comparison_class := factor(
    comparison_class,
    levels = c(
      "Nenhum",
      "GSE53697",
      "GSE125583",
      "Ambos"
    )
  )
]

effect_class_counts <- valid_effects[
  ,
  .N,
  by = comparison_class
]

p_effect <- ggplot(
  valid_effects,
  aes(
    x = old_log2FC,
    y = new_log2FC,
    color = comparison_class
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.15
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "grey35",
    linewidth = 0.65
  ) +
  scale_color_manual(
    values = EFFECT_CLASS_COLORS,
    drop = FALSE
  ) +
  labs(
    title = "Efeito de AD — GSE53697 vs GSE125583",
    subtitle = paste0(
      "Genes testáveis em ambos; Spearman rho = ",
      format_num(
        rho_old_new,
        digits = 3
      )
    ),
    x = "log2FC — GSE53697",
    y = "log2FC — GSE125583",
    color = paste0(
      "Critério nominal\np < 0,05 e |log2FC| > ",
      SENSITIVITY_LFC
    )
  ) +
  comparison_theme()

ggsave(
  filename = file.path(
    plot_dir,
    "AD_old_vs_new_log2FC.png"
  ),
  plot = p_effect,
  width = 8.5,
  height = 7,
  dpi = 300
)


# ==============================================================================
# 21. RELATÓRIO TEXTUAL
# ==============================================================================

old_summary <- deg_comparison[
  dataset == DATASET_OLD
]

new_summary <- deg_comparison[
  dataset == DATASET_NEW
]

sink(
  audit_file
)

cat("============================================================\n")
cat("07 — COMPARAÇÃO DOS DATASETS DE AD\n")
cat("GSE53697 vs GSE125583\n")
cat("============================================================\n\n")


cat("------------------------------------------------------------\n")
cat("1. TAMANHO AMOSTRAL E TECIDO\n")
cat("------------------------------------------------------------\n")

cat(
  DATASET_OLD,
  ": ",
  format_int(
    n_old_samples
  ),
  " amostras (",
  format_int(
    n_old_ad
  ),
  " AD / ",
  format_int(
    n_old_control
  ),
  " Control)\n",
  sep = ""
)

cat(
  "Tecido: ",
  TISSUE_OLD,
  "\n\n",
  sep = ""
)

cat(
  DATASET_NEW,
  ": ",
  format_int(
    n_new_samples
  ),
  " amostras (",
  format_int(
    n_new_ad
  ),
  " AD / ",
  format_int(
    n_new_control
  ),
  " Control)\n",
  sep = ""
)

cat(
  "Tecido: ",
  TISSUE_NEW,
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("2. UNIVERSO GÊNICO BRUTO\n")
cat("------------------------------------------------------------\n")

cat(
  DATASET_OLD,
  ": ",
  format_int(
    length(
      old_raw_genes
    )
  ),
  "\n",
  sep = ""
)

cat(
  DATASET_NEW,
  ": ",
  format_int(
    length(
      new_raw_genes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Interseção bruta: ",
  format_int(
    length(
      raw_intersection
    )
  ),
  "\n",
  sep = ""
)

cat(
  "União bruta: ",
  format_int(
    length(
      raw_union
    )
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("3. UNIVERSO TESTÁVEL\n")
cat("------------------------------------------------------------\n")

cat(
  DATASET_OLD,
  ": ",
  format_int(
    length(
      old_testable_genes
    )
  ),
  "\n",
  sep = ""
)

cat(
  DATASET_NEW,
  ": ",
  format_int(
    length(
      new_testable_genes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Interseção testável: ",
  format_int(
    length(
      testable_intersection
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Somente ",
  DATASET_OLD,
  ": ",
  format_int(
    length(
      old_testable_only
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Somente ",
  DATASET_NEW,
  ": ",
  format_int(
    length(
      new_testable_only
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Jaccard do universo testável: ",
  format_num(
    testable_jaccard,
    digits = 4
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("4. EXPRESSÃO DIFERENCIAL — GSE53697\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes testados: ",
  format_int(
    old_summary$n_tested
  ),
  "\n",
  sep = ""
)

cat(
  "FDR < 0,05: ",
  format_int(
    old_summary$n_FDR_lt_005
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs principais: ",
  format_int(
    old_summary$n_primary
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs sensibilidade: ",
  format_int(
    old_summary$n_sensitivity
  ),
  "\n",
  sep = ""
)

cat(
  "Article-like: ",
  format_int(
    old_summary$n_article_like
  ),
  "\n",
  sep = ""
)

cat(
  "Menor p-value: ",
  format_p(
    old_summary$min_pvalue
  ),
  "\n",
  sep = ""
)

cat(
  "Menor padj: ",
  format_p(
    old_summary$min_padj
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("5. EXPRESSÃO DIFERENCIAL — GSE125583\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes testados: ",
  format_int(
    new_summary$n_tested
  ),
  "\n",
  sep = ""
)

cat(
  "FDR < 0,05: ",
  format_int(
    new_summary$n_FDR_lt_005
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs principais: ",
  format_int(
    new_summary$n_primary
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs principais Up: ",
  format_int(
    new_summary$n_primary_up
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs principais Down: ",
  format_int(
    new_summary$n_primary_down
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs sensibilidade: ",
  format_int(
    new_summary$n_sensitivity
  ),
  "\n",
  sep = ""
)

cat(
  "Article-like: ",
  format_int(
    new_summary$n_article_like
  ),
  "\n",
  sep = ""
)

cat(
  "Menor p-value: ",
  format_p(
    new_summary$min_pvalue
  ),
  "\n",
  sep = ""
)

cat(
  "Menor padj: ",
  format_p(
    new_summary$min_padj
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("6. CONCORDÂNCIA GENE-A-GENE ENTRE OS DOIS DATASETS\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes testáveis em ambos: ",
  format_int(
    nrow(
      common_comparison
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes com log2FC comparável: ",
  format_int(
    nrow(
      valid_effects
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Spearman rho dos log2FC: ",
  format_num(
    rho_old_new,
    digits = 4
  ),
  "\n",
  sep = ""
)

cat(
  "Genes article-like apenas no GSE53697: ",
  format_int(
    sum(
      valid_effects$comparison_class == "GSE53697",
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes article-like apenas no GSE125583: ",
  format_int(
    sum(
      valid_effects$comparison_class == "GSE125583",
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes article-like em ambos: ",
  format_int(
    sum(
      valid_effects$comparison_class == "Ambos",
      na.rm = TRUE
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "NOTA: esses valores são restritos aos genes testáveis nos dois datasets. ",
  "Por isso, diferem dos totais article-like das análises individuais ",
  "(52 no GSE53697 e 155 no GSE125583).\n\n"
)

cat(
  "IMPORTANTE: a correlação e a sobreposição gene-a-gene NÃO devem ser ",
  "interpretadas como teste de replicação direta, pois os estudos utilizam ",
  "regiões cerebrais, tamanhos amostrais e coortes diferentes.\n\n"
)


cat("------------------------------------------------------------\n")
cat("7. SOBREPOSIÇÃO DE DEGs\n")
cat("------------------------------------------------------------\n")

cat(
  "Critério principal — ",
  DATASET_OLD,
  ": ",
  format_int(
    length(
      old_primary_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Critério principal — ",
  DATASET_NEW,
  ": ",
  format_int(
    length(
      new_primary_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Critério principal — interseção: ",
  format_int(
    length(
      primary_overlap
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Article-like — ",
  DATASET_OLD,
  ": ",
  format_int(
    length(
      old_article_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Article-like — ",
  DATASET_NEW,
  ": ",
  format_int(
    length(
      new_article_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Article-like — interseção: ",
  format_int(
    length(
      article_overlap
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Article-like compartilhados com mesma direção: ",
  format_int(
    n_article_same_direction
  ),
  "\n",
  sep = ""
)

cat(
  "Article-like compartilhados com direção oposta: ",
  format_int(
    n_article_opposite_direction
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("8. COMPATIBILIDADE DO GSE125583 COM HD E PD\n")
cat("------------------------------------------------------------\n")

cat(
  "AD novo testável: ",
  format_int(
    length(
      new_testable_genes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "HD testável: ",
  format_int(
    length(
      hd_testable_genes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "PD testável: ",
  format_int(
    length(
      pd_testable_genes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "AD novo ∩ HD: ",
  format_int(
    length(
      candidate_hd_intersection
    )
  ),
  "\n",
  sep = ""
)

cat(
  "AD novo ∩ PD: ",
  format_int(
    length(
      candidate_pd_intersection
    )
  ),
  "\n",
  sep = ""
)

cat(
  "HD ∩ PD: ",
  format_int(
    length(
      hd_pd_intersection
    )
  ),
  "\n",
  sep = ""
)

cat(
  "AD novo ∩ HD ∩ PD: ",
  format_int(
    length(
      candidate_hd_pd_intersection
    )
  ),
  "\n",
  sep = ""
)

cat(
  "AD antigo ∩ HD ∩ PD: ",
  format_int(
    length(
      old_hd_pd_intersection
    )
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("9. INTERPRETAÇÃO PARA SELEÇÃO DO DATASET\n")
cat("------------------------------------------------------------\n")

cat(
  "- O GSE53697 é preservado integralmente como resultado da análise original; nenhum resultado é apagado ou substituído retroativamente.\n"
)

cat(
  "- O GSE125583 foi avaliado como um dataset candidato independente, utilizando matriz raw_counts_GRCh38.p13_NCBI e pipeline documentado nas etapas 00–06 da Fase 01b.\n"
)

cat(
  "- A escolha do GSE125583 não deve ser justificada apenas pelo maior número de DEGs.\n"
)

cat(
  "- A decisão deve considerar em conjunto tamanho amostral, universo testável, robustez após FDR, estabilidade após ajuste por idade/sexo, coerência funcional e compatibilidade com HD/PD.\n"
)

cat(
  "- A comparação gene-a-gene entre GSE53697 e GSE125583 não constitui replicação direta, pois os estudos utilizam regiões cerebrais distintas.\n"
)

cat(
  "- Baixa sobreposição de DEGs entre os dois estudos, se observada, deve ser interpretada à luz dessas diferenças biológicas e técnicas.\n"
)

cat(
  "- Para a análise transdiagnóstica, o GSE125583 apresenta maior compatibilidade técnica e gênica com os datasets de HD e PD, incluindo maior universo testável compartilhado.\n"
)

cat(
  "- Entretanto, essa vantagem NÃO se estende à região cerebral: GSE64810 (HD) e GSE68719 (PD) utilizam córtex pré-frontal BA9, enquanto GSE125583 utiliza giro fusiforme.\n"
)

cat(
  "- Portanto, a adoção do GSE125583 melhora a comparabilidade técnica e estatística, mas introduz maior heterogeneidade anatômica na análise transdiagnóstica.\n"
)

cat(
  "- A análise transdiagnóstica deverá comparar efeitos relativos dentro de cada estudo (caso vs controle), e não níveis absolutos de expressão entre regiões cerebrais distintas.\n"
)

cat(
  "- A diferença de região cerebral permanece uma limitação biológica importante e deve ser explicitamente considerada na interpretação de expressão, composição celular e deconvolução.\n"
)

cat(
  "- O GSE125583 deve ser descrito como dataset candidato tecnicamente e estatisticamente mais adequado para a etapa transdiagnóstica, mas com a importante limitação de utilizar uma região cerebral distinta daquela empregada nos datasets de HD e PD.\n"
)

cat(
  "- A decisão final sobre substituição do dataset deve ser registrada no relatório da Fase 01b e mantida separada da análise original.\n"
)

cat("\n")
cat("============================================================\n")
cat("07_compare_AD_datasets.R CONCLUÍDO\n")
cat("============================================================\n")

sink()


# ==============================================================================
# 22. RESUMO NO TERMINAL
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("07 — COMPARAÇÃO AD CONCLUÍDA\n")
cat("============================================================\n")

cat(
  DATASET_OLD,
  " — amostras/testáveis/DEGs principais: ",
  format_int(
    n_old_samples
  ),
  "/",
  format_int(
    length(
      old_testable_genes
    )
  ),
  "/",
  format_int(
    old_summary$n_primary
  ),
  "\n",
  sep = ""
)

cat(
  DATASET_NEW,
  " — amostras/testáveis/DEGs principais: ",
  format_int(
    n_new_samples
  ),
  "/",
  format_int(
    length(
      new_testable_genes
    )
  ),
  "/",
  format_int(
    new_summary$n_primary
  ),
  "\n",
  sep = ""
)

cat(
  "AD antigo ∩ AD novo — genes testáveis: ",
  format_int(
    length(
      testable_intersection
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Spearman log2FC antigo vs novo: ",
  format_num(
    rho_old_new,
    digits = 3
  ),
  "\n",
  sep = ""
)

cat(
  "AD novo ∩ HD ∩ PD — universo testável: ",
  format_int(
    length(
      candidate_hd_pd_intersection
    )
  ),
  "\n",
  sep = ""
)

cat(
  "AD antigo ∩ HD ∩ PD — universo testável: ",
  format_int(
    length(
      old_hd_pd_intersection
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Relatório: ",
  audit_file,
  "\n",
  sep = ""
)

cat("============================================================\n")
