#!/usr/bin/env Rscript

# ==============================================================================
# 05_deseq2_GSE125583.R
#
# FASE 01b — EXPRESSÃO DIFERENCIAL DO DATASET CANDIDATO DE AD
# GEO: GSE125583
#
# OBJETIVO GERAL
# ------------------------------------------------------------------------------
# Executar a análise de expressão diferencial do GSE125583 com DESeq2 após o
# pré-filtro definido na etapa 04 e avaliar a robustez do sinal de AD.
#
# Esta etapa possui DOIS modelos:
#
#   MODELO PRINCIPAL
#   ---------------------------------------------------------------------------
#   design = ~ group
#
#   - utiliza todas as 289 amostras;
#   - Control é a referência;
#   - o contraste é AD vs Control;
#   - mantém a mesma estrutura usada na análise principal dos datasets de HD e
#     PD, favorecendo a comparabilidade transdiagnóstica.
#
#
#   MODELO AJUSTADO / ANÁLISE DE SENSIBILIDADE
#   ---------------------------------------------------------------------------
#   design = ~ age_centered + sex + group
#
#   - utiliza apenas amostras com idade e sexo conhecidos;
#   - serve para verificar se o efeito de AD permanece estável após ajuste por
#     idade e sexo;
#   - NÃO substitui o modelo principal;
#   - Braak NÃO entra no modelo, pois está fortemente relacionado à própria
#     neuropatologia/diagnóstico e poderia remover parte do sinal de interesse.
#
#
# CRITÉRIOS DE CLASSIFICAÇÃO
# ------------------------------------------------------------------------------
# O DESeq2 estima log2FoldChange, p-value e padj normalmente.
# NÃO usamos lfcThreshold dentro de results().
#
# A classificação é feita DEPOIS da estimação:
#
# 1. Critério principal:
#
#      padj < 0,05  E  |log2FoldChange| > 0,58
#
# 2. Análise de sensibilidade:
#
#      padj < 0,05  E  |log2FoldChange| > 1
#
# 3. Critério nominal comparativo ("article-like"):
#
#      pvalue < 0,05  E  |log2FoldChange| > 1
#
# O terceiro critério é exclusivamente comparativo/exploratório e não substitui
# o critério principal baseado em FDR.
#
#
# INDEPENDENT FILTERING
# ------------------------------------------------------------------------------
# results(..., independentFiltering = TRUE, alpha = 0.05)
#
# Isso permite que o DESeq2 aplique o independent filtering padrão ao cálculo
# de padj. Genes com p-value válido podem, em alguns casos, ficar com padj = NA.
#
#
# VST E HEATMAP
# ------------------------------------------------------------------------------
# VST é utilizado apenas para visualização do heatmap.
# Os testes de expressão diferencial usam raw counts.
#
#
# O QUE ESTE SCRIPT FAZ
# ------------------------------------------------------------------------------
# 1. Lê counts filtrados da etapa 04 e metadata da etapa 01.
# 2. Valida novamente counts, GeneID e correspondência counts <-> metadata.
# 3. Executa DESeq2 principal: ~ group.
# 4. Extrai resultados AD vs Control.
# 5. Classifica genes pelos três critérios acima.
# 6. Integra a anotação do gene master da etapa 03.
# 7. Salva:
#      - resultados completos;
#      - DEGs do critério principal;
#      - DEGs da análise de sensibilidade;
#      - genes do critério nominal comparativo.
# 8. Gera:
#      - MA plot;
#      - volcano plot;
#      - heatmap dos até 50 principais DEGs.
# 9. Executa o modelo ajustado ~ age + sex + group.
# 10. Compara log2FC e classificação entre modelo principal e ajustado.
# 11. Gera um scatter de log2FC principal vs ajustado.
# 12. Gera resumo TSV e relatório textual.
#
#
# O QUE ESTE SCRIPT NÃO FAZ
# ------------------------------------------------------------------------------
# - NÃO remove amostras com base em PCA.
# - NÃO usa Braak como covariável.
# - NÃO altera o pré-filtro da etapa 04.
# - NÃO executa enriquecimento funcional.
# - NÃO integra ainda AD com HD e PD.
#
#
# ENTRADAS
# ------------------------------------------------------------------------------
# results/candidate_GSE125583/qc/counts_GSE125583_DE_filtered.tsv
#
# results/candidate_GSE125583/metadata/metadata_GSE125583.tsv
#
# results/candidate_GSE125583/harmonization/
#   gene_master_GSE125583_HD_PD.tsv
#
#
# SAÍDAS PRINCIPAIS
# ------------------------------------------------------------------------------
# results/candidate_GSE125583/deseq2/
# ├── GSE125583_deseq2_results.tsv
# ├── GSE125583_DEG_primary.tsv
# ├── GSE125583_DEG_sensitivity.tsv
# ├── GSE125583_DEG_article_like.tsv
# ├── GSE125583_deseq2_adjusted_results.tsv
# ├── GSE125583_main_vs_adjusted.tsv
# ├── deseq2_summary_GSE125583.tsv
# └── 05_deseq2_GSE125583.txt
#
# results/candidate_GSE125583/plots/deseq2/
# ├── MA_GSE125583.png
# ├── volcano_GSE125583.png
# ├── heatmap_GSE125583_top50_primary.png      [se houver DEGs principais]
# └── log2FC_main_vs_adjusted_GSE125583.png
# ==============================================================================


# ==============================================================================
# 1. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table",
  "DESeq2",
  "ggplot2",
  "pheatmap",
  "ggrepel"
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
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
})


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

counts_file <- file.path(
  "results",
  "candidate_GSE125583",
  "qc",
  "counts_GSE125583_DE_filtered.tsv"
)

metadata_file <- file.path(
  "results",
  "candidate_GSE125583",
  "metadata",
  "metadata_GSE125583.tsv"
)

gene_master_file <- file.path(
  "results",
  "candidate_GSE125583",
  "harmonization",
  "gene_master_GSE125583_HD_PD.tsv"
)

output_dir <- file.path(
  "results",
  "candidate_GSE125583",
  "deseq2"
)

plot_dir <- file.path(
  "results",
  "candidate_GSE125583",
  "plots",
  "deseq2"
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

main_results_file <- file.path(
  output_dir,
  "GSE125583_deseq2_results.tsv"
)

primary_deg_file <- file.path(
  output_dir,
  "GSE125583_DEG_primary.tsv"
)

sensitivity_deg_file <- file.path(
  output_dir,
  "GSE125583_DEG_sensitivity.tsv"
)

article_like_file <- file.path(
  output_dir,
  "GSE125583_DEG_article_like.tsv"
)

adjusted_results_file <- file.path(
  output_dir,
  "GSE125583_deseq2_adjusted_results.tsv"
)

model_comparison_file <- file.path(
  output_dir,
  "GSE125583_main_vs_adjusted.tsv"
)

summary_file <- file.path(
  output_dir,
  "deseq2_summary_GSE125583.tsv"
)

audit_file <- file.path(
  output_dir,
  "05_deseq2_GSE125583.txt"
)

ma_plot_file <- file.path(
  plot_dir,
  "MA_GSE125583.png"
)

volcano_plot_file <- file.path(
  plot_dir,
  "volcano_GSE125583.png"
)

heatmap_plot_file <- file.path(
  plot_dir,
  "heatmap_GSE125583_top50_primary.png"
)

model_comparison_plot_file <- file.path(
  plot_dir,
  "log2FC_main_vs_adjusted_GSE125583.png"
)


# ==============================================================================
# 3. CONSTANTES ANALÍTICAS
# ==============================================================================

EXPECTED_N_SAMPLES <- 289L
EXPECTED_N_AD <- 219L
EXPECTED_N_CONTROL <- 70L

ALPHA <- 0.05
PRIMARY_LFC <- 0.58
SENSITIVITY_LFC <- 1.00

GROUP_LEVELS <- c(
  "Control",
  "AD"
)

PRIMARY_COLORS <- c(
  "Down" = "#2166AC",
  "Not_DEG" = "#BDBDBD",
  "Up" = "#B2182B"
)


# ==============================================================================
# 4. FUNÇÕES AUXILIARES
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 Validar arquivo
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# 4.2 Validar GeneID
# ------------------------------------------------------------------------------

validate_gene_id <- function(
  counts,
  dataset_name
) {

  if (!"GeneID" %in% names(counts)) {
    stop(
      dataset_name,
      ": coluna GeneID não encontrada."
    )
  }

  if (names(counts)[1] != "GeneID") {
    stop(
      dataset_name,
      ": GeneID deve ser a primeira coluna."
    )
  }

  gene_id <- as.character(
    counts$GeneID
  )

  missing_gene_id <- (
    is.na(gene_id) |
      trimws(gene_id) == ""
  )

  if (any(missing_gene_id)) {
    stop(
      dataset_name,
      ": existem GeneIDs ausentes."
    )
  }

  if (anyDuplicated(gene_id)) {
    stop(
      dataset_name,
      ": existem GeneIDs duplicados."
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 4.3 Validar metadata
# ------------------------------------------------------------------------------

validate_metadata <- function(metadata) {

  required_columns <- c(
    "sample_id",
    "group",
    "age",
    "sex",
    "braak_score"
  )

  missing_columns <- setdiff(
    required_columns,
    names(metadata)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Metadata sem colunas obrigatórias: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  if (any(is.na(metadata$sample_id))) {
    stop(
      "Existem sample IDs ausentes."
    )
  }

  if (anyDuplicated(metadata$sample_id)) {
    stop(
      "Existem sample IDs duplicados."
    )
  }

  if (any(is.na(metadata$group))) {
    stop(
      "Existem grupos ausentes."
    )
  }

  unexpected_groups <- setdiff(
    unique(
      metadata$group
    ),
    GROUP_LEVELS
  )

  if (length(unexpected_groups) > 0) {
    stop(
      "Grupos inesperados: ",
      paste(
        unexpected_groups,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 4.4 Validar counts <-> metadata
# ------------------------------------------------------------------------------

validate_sample_mapping <- function(
  counts,
  metadata,
  dataset_name
) {

  count_samples <- setdiff(
    names(counts),
    "GeneID"
  )

  metadata_samples <- metadata$sample_id

  only_counts <- setdiff(
    count_samples,
    metadata_samples
  )

  only_metadata <- setdiff(
    metadata_samples,
    count_samples
  )

  if (length(only_counts) > 0) {
    stop(
      dataset_name,
      ": amostras somente nos counts: ",
      paste(
        only_counts,
        collapse = ", "
      )
    )
  }

  if (length(only_metadata) > 0) {
    stop(
      dataset_name,
      ": amostras somente no metadata: ",
      paste(
        only_metadata,
        collapse = ", "
      )
    )
  }

  if (length(count_samples) != nrow(metadata)) {
    stop(
      dataset_name,
      ": número de amostras diferente entre counts e metadata."
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 4.5 Formatação pt-BR
# ------------------------------------------------------------------------------

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
# 4.6 Tema dos gráficos
# ------------------------------------------------------------------------------

plot_theme <- function() {

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


# ------------------------------------------------------------------------------
# 4.7 Classificar resultados do DESeq2
# ------------------------------------------------------------------------------

classify_results <- function(res_dt) {

  res_dt[
    ,
    fdr_significant := (
      !is.na(padj) &
        padj < ALPHA
    )
  ]

  res_dt[
    ,
    is_primary := (
      !is.na(padj) &
        padj < ALPHA &
        !is.na(log2FoldChange) &
        abs(log2FoldChange) > PRIMARY_LFC
    )
  ]

  res_dt[
    ,
    is_sensitivity := (
      !is.na(padj) &
        padj < ALPHA &
        !is.na(log2FoldChange) &
        abs(log2FoldChange) > SENSITIVITY_LFC
    )
  ]

  res_dt[
    ,
    is_article_like := (
      !is.na(pvalue) &
        pvalue < ALPHA &
        !is.na(log2FoldChange) &
        abs(log2FoldChange) > SENSITIVITY_LFC
    )
  ]

  res_dt[
    ,
    primary_direction := fifelse(
      is_primary &
        log2FoldChange > 0,
      "Up",
      fifelse(
        is_primary &
          log2FoldChange < 0,
        "Down",
        "Not_DEG"
      )
    )
  ]

  res_dt[
    ,
    sensitivity_direction := fifelse(
      is_sensitivity &
        log2FoldChange > 0,
      "Up",
      fifelse(
        is_sensitivity &
          log2FoldChange < 0,
        "Down",
        "Not_DEG"
      )
    )
  ]

  res_dt[
    ,
    article_like_direction := fifelse(
      is_article_like &
        log2FoldChange > 0,
      "Up",
      fifelse(
        is_article_like &
          log2FoldChange < 0,
        "Down",
        "Not_DEG"
      )
    )
  ]

  res_dt
}


# ------------------------------------------------------------------------------
# 4.8 Converter DESeqResults em data.table
# ------------------------------------------------------------------------------

results_to_dt <- function(
  res,
  model_name
) {

  dt <- data.table(
    GeneID = rownames(
      res
    ),
    baseMean = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    lfcSE = res$lfcSE,
    stat = res$stat,
    pvalue = res$pvalue,
    padj = res$padj
  )

  dt[
    ,
    model := model_name
  ]

  classify_results(
    dt
  )
}


# ------------------------------------------------------------------------------
# 4.9 Integrar anotação
# ------------------------------------------------------------------------------

annotate_results <- function(
  results_dt,
  annotation_dt
) {

  annotation_columns <- c(
    "GeneID",
    "Symbol",
    "Description",
    "GeneType",
    "EnsemblGeneID"
  )

  missing_cols <- setdiff(
    annotation_columns,
    names(
      annotation_dt
    )
  )

  if (length(missing_cols) > 0) {
    stop(
      "Gene master sem colunas de anotação esperadas: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }

  original_order <- results_dt$GeneID

  out <- merge(
    results_dt,
    annotation_dt[
      ,
      ..annotation_columns
    ],
    by = "GeneID",
    all.x = TRUE,
    sort = FALSE
  )

  out <- out[
    match(
      original_order,
      GeneID
    )
  ]

  if (!identical(
    out$GeneID,
    original_order
  )) {
    stop(
      "Falha ao restaurar a ordem após anotação."
    )
  }

  setcolorder(
    out,
    c(
      "GeneID",
      "Symbol",
      "Description",
      "GeneType",
      "EnsemblGeneID",
      setdiff(
        names(out),
        c(
          "GeneID",
          "Symbol",
          "Description",
          "GeneType",
          "EnsemblGeneID"
        )
      )
    )
  )

  out
}


# ------------------------------------------------------------------------------
# 4.10 Resumo de um conjunto de resultados
# ------------------------------------------------------------------------------

summarize_results <- function(
  dt,
  model_label,
  n_samples
) {

  data.table(
    model = model_label,
    n_samples = n_samples,
    n_genes = nrow(dt),

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

    n_DEG_primary = sum(
      dt$is_primary,
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

    n_DEG_sensitivity = sum(
      dt$is_sensitivity,
      na.rm = TRUE
    ),

    n_sensitivity_up = sum(
      dt$sensitivity_direction == "Up",
      na.rm = TRUE
    ),

    n_sensitivity_down = sum(
      dt$sensitivity_direction == "Down",
      na.rm = TRUE
    ),

    n_article_like = sum(
      dt$is_article_like,
      na.rm = TRUE
    ),

    n_article_like_up = sum(
      dt$article_like_direction == "Up",
      na.rm = TRUE
    ),

    n_article_like_down = sum(
      dt$article_like_direction == "Down",
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


# ==============================================================================
# 5. VALIDAR ARQUIVOS
# ==============================================================================

input_files <- c(
  counts_file,
  metadata_file,
  gene_master_file
)

invisible(
  lapply(
    input_files,
    validate_file
  )
)


# ==============================================================================
# 6. LER METADATA
# ==============================================================================

metadata <- fread(
  metadata_file,
  na.strings = c(
    "NA",
    ""
  ),
  showProgress = FALSE
)

validate_metadata(
  metadata
)

if (nrow(metadata) != EXPECTED_N_SAMPLES) {
  stop(
    "Número inesperado de amostras no metadata. Esperado: ",
    EXPECTED_N_SAMPLES,
    "; observado: ",
    nrow(metadata),
    "."
  )
}

n_ad <- metadata[
  group == "AD",
  .N
]

n_control <- metadata[
  group == "Control",
  .N
]

if (n_ad != EXPECTED_N_AD) {
  stop(
    "Número inesperado de AD. Esperado: ",
    EXPECTED_N_AD,
    "; observado: ",
    n_ad,
    "."
  )
}

if (n_control != EXPECTED_N_CONTROL) {
  stop(
    "Número inesperado de controles. Esperado: ",
    EXPECTED_N_CONTROL,
    "; observado: ",
    n_control,
    "."
  )
}


# ==============================================================================
# 7. LER COUNTS FILTRADOS DA ETAPA 04
# ==============================================================================

counts <- fread(
  counts_file,
  showProgress = TRUE
)

counts[
  ,
  GeneID := as.character(
    GeneID
  )
]

validate_gene_id(
  counts,
  "GSE125583 filtered"
)

validate_sample_mapping(
  counts,
  metadata,
  "GSE125583 filtered"
)

# Reordenar explicitamente na ordem do metadata.
counts <- counts[
  ,
  c(
    "GeneID",
    metadata$sample_id
  ),
  with = FALSE
]

sample_ids <- metadata$sample_id

count_matrix <- as.matrix(
  counts[
    ,
    ..sample_ids
  ]
)

rownames(
  count_matrix
) <- counts$GeneID

if (anyNA(count_matrix)) {
  stop(
    "Matriz filtrada contém valores ausentes."
  )
}

if (any(!is.finite(count_matrix))) {
  stop(
    "Matriz filtrada contém valores não finitos."
  )
}

if (any(count_matrix < 0)) {
  stop(
    "Matriz filtrada contém counts negativos."
  )
}

if (any(
  abs(
    count_matrix -
      round(
        count_matrix
      )
  ) > 1e-8
)) {
  stop(
    "Matriz filtrada contém counts fracionários."
  )
}

storage.mode(
  count_matrix
) <- "integer"


# ==============================================================================
# 8. VALIDAR QUE O PRÉ-FILTRO DA ETAPA 04 FOI RESPEITADO
# ==============================================================================

n_samples_ge_10 <- rowSums(
  count_matrix >= 10
)

if (any(
  n_samples_ge_10 < 70
)) {
  stop(
    "A matriz filtrada contém genes que não satisfazem ",
    "count >= 10 em pelo menos 70 amostras."
  )
}

n_genes_tested <- nrow(
  count_matrix
)


# ==============================================================================
# 9. LER GENE MASTER / ANOTAÇÃO
# ==============================================================================

gene_master <- fread(
  gene_master_file,
  na.strings = c(
    "NA",
    ""
  ),
  showProgress = FALSE
)

gene_master[
  ,
  GeneID := as.character(
    GeneID
  )
]

if (anyDuplicated(
  gene_master$GeneID
)) {
  stop(
    "Gene master contém GeneIDs duplicados."
  )
}

genes_without_master <- setdiff(
  rownames(
    count_matrix
  ),
  gene_master$GeneID
)

if (length(genes_without_master) > 0) {
  stop(
    "Existem genes da matriz filtrada ausentes no gene master: ",
    length(
      genes_without_master
    )
  )
}


# ==============================================================================
# 10. PREPARAR colData DO MODELO PRINCIPAL
# ==============================================================================

col_data_main <- as.data.frame(
  metadata
)

rownames(
  col_data_main
) <- col_data_main$sample_id

col_data_main$group <- factor(
  col_data_main$group,
  levels = GROUP_LEVELS
)

if (!identical(
  rownames(
    col_data_main
  ),
  colnames(
    count_matrix
  )
)) {
  stop(
    "A ordem de colData e count matrix não coincide."
  )
}


# ==============================================================================
# 11. MODELO PRINCIPAL — design ~ group
# ==============================================================================

dds_main <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data_main,
  design = ~ group
)

# DESeq() executa:
#   - estimação de size factors;
#   - estimação de dispersões;
#   - ajuste do modelo binomial negativo;
#   - teste de Wald.
#
# Mantemos os defaults do DESeq2 para comparabilidade com a análise anterior.
dds_main <- DESeq(
  dds_main,
  quiet = TRUE
)

main_results_raw <- results(
  dds_main,
  contrast = c(
    "group",
    "AD",
    "Control"
  ),
  alpha = ALPHA,
  independentFiltering = TRUE
)

main_results <- results_to_dt(
  main_results_raw,
  model_name = "main_group"
)

main_results <- annotate_results(
  main_results,
  gene_master
)


# ==============================================================================
# 12. EXTRAIR LISTAS DO MODELO PRINCIPAL
# ==============================================================================

primary_deg <- main_results[
  is_primary == TRUE
][
  order(
    padj,
    -abs(
      log2FoldChange
    )
  )
]

sensitivity_deg <- main_results[
  is_sensitivity == TRUE
][
  order(
    padj,
    -abs(
      log2FoldChange
    )
  )
]

article_like_deg <- main_results[
  is_article_like == TRUE
][
  order(
    pvalue,
    -abs(
      log2FoldChange
    )
  )
]


# ==============================================================================
# 13. SALVAR RESULTADOS DO MODELO PRINCIPAL
# ==============================================================================

fwrite(
  main_results,
  main_results_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  primary_deg,
  primary_deg_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  sensitivity_deg,
  sensitivity_deg_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  article_like_deg,
  article_like_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 14. SIZE FACTORS DO MODELO PRINCIPAL
# ==============================================================================

main_size_factors <- sizeFactors(
  dds_main
)


# ==============================================================================
# 15. MA PLOT
# ==============================================================================

ma_data <- copy(
  main_results
)

ma_data[
  ,
  baseMean_plot := pmax(
    baseMean,
    1e-6
  )
]

p_ma <- ggplot(
  ma_data[
    !is.na(
      log2FoldChange
    )
  ],
  aes(
    x = baseMean_plot,
    y = log2FoldChange,
    color = primary_direction
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.2
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = c(
      -PRIMARY_LFC,
      PRIMARY_LFC
    ),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_x_log10() +
  scale_color_manual(
    values = PRIMARY_COLORS,
    breaks = c(
      "Up",
      "Down",
      "Not_DEG"
    )
  ) +
  labs(
    title = "MA plot — GSE125583",
    subtitle = paste0(
      "AD vs Control; principal: FDR < 0,05 e |log2FC| > ",
      PRIMARY_LFC
    ),
    x = "Mean normalized counts (baseMean, log10)",
    y = "log2FoldChange",
    color = "Classificação"
  ) +
  plot_theme()

ggsave(
  filename = ma_plot_file,
  plot = p_ma,
  width = 8,
  height = 6,
  dpi = 300
)


# ==============================================================================
# 16. VOLCANO PLOT
# ==============================================================================
#
# O volcano plot combina:
#
#   eixo x -> magnitude do efeito (log2FoldChange)
#   eixo y -> significância estatística (-log10 do padj)
#
# Os genes são classificados segundo o critério principal:
#
#   FDR < 0,05 e |log2FC| > 0,58
#
# Para manter o padrão visual utilizado nos volcano plots de HD e PD,
# rotulamos os 10 DEGs principais com menor padj.
# ==============================================================================


# ------------------------------------------------------------------------------
# 16.1 Preparar dados do volcano
# ------------------------------------------------------------------------------

volcano_data <- copy(
  main_results
)


# Para visualização, genes sem padj ou log2FC não entram no volcano.
volcano_data <- volcano_data[
  !is.na(
    padj
  ) &
    !is.na(
      log2FoldChange
    )
]


# ------------------------------------------------------------------------------
# 16.2 Evitar valores infinitos em -log10(padj)
# ------------------------------------------------------------------------------
#
# Se algum padj for numericamente igual a zero, -log10(0) produziria Inf.
# Substituímos apenas para fins gráficos pelo menor padj positivo observado.

min_positive_padj <- min(
  volcano_data$padj[
    volcano_data$padj > 0
  ],
  na.rm = TRUE
)

if (
  !is.finite(
    min_positive_padj
  )
) {
  min_positive_padj <- .Machine$double.xmin
}


volcano_data[
  ,
  padj_plot := pmax(
    padj,
    min_positive_padj
  )
]


volcano_data[
  ,
  minus_log10_padj := -log10(
    padj_plot
  )
]


# ------------------------------------------------------------------------------
# 16.3 Selecionar genes para rotular
# ------------------------------------------------------------------------------
#
# Selecionamos somente genes classificados como DEG pelo critério principal.
#
# Critério de ordenação:
#   1. menor padj;
#   2. maior |log2FoldChange| em caso de empate.
#
# Apenas genes com Symbol disponível são considerados.
#
# A seleção é objetiva e reproduzível, evitando escolher genes manualmente.

volcano_labels <- volcano_data[
  primary_direction != "Not_DEG" &
    !is.na(Symbol) &
    trimws(Symbol) != ""
]


# order() permite utilizar expressões como abs(log2FoldChange),
# ao contrário de setorder(), que trabalha diretamente com nomes de colunas.
volcano_labels <- volcano_labels[
  order(
    padj,
    -abs(log2FoldChange)
  )
]


# Manter apenas os 10 DEGs principais mais significativos.
volcano_labels <- head(
  volcano_labels,
  10L
)

# ------------------------------------------------------------------------------
# 16.4 Construir volcano plot
# ------------------------------------------------------------------------------

p_volcano <- ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = minus_log10_padj,
    color = primary_direction
  )
) +

  geom_point(
    alpha = 0.55,
    size = 1.3
  ) +

  # Limiares de magnitude do efeito
  geom_vline(
    xintercept = c(
      -PRIMARY_LFC,
      PRIMARY_LFC
    ),
    linetype = "dashed",
    linewidth = 0.5
  ) +

  # Limiar de FDR = 0,05
  geom_hline(
    yintercept = -log10(
      ALPHA
    ),
    linetype = "dashed",
    linewidth = 0.5
  ) +

  scale_color_manual(
    values = PRIMARY_COLORS,
    breaks = c(
      "Down",
      "Not_DEG",
      "Up"
    )
  ) +

  # --------------------------------------------------------------------------
  # Rótulos dos top 10 DEGs
  # --------------------------------------------------------------------------
  #
  # geom_text_repel reposiciona os nomes automaticamente para reduzir
  # sobreposição entre os rótulos.

  ggrepel::geom_text_repel(
    data = volcano_labels,
    aes(
      label = Symbol
    ),
    size = 3.8,
    box.padding = 0.45,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123,
    show.legend = FALSE
  ) +

  labs(
    title = "Volcano plot — GSE125583",

    subtitle = paste0(
      "AD vs Control; principal: FDR < 0,05 e |log2FC| > ",
      PRIMARY_LFC
    ),

    x = "log2FoldChange",

    # Mesmo significado de FDR/padj, mas mais próximo da nomenclatura
    # usada nos gráficos anteriores de HD e PD.
    y = expression(
      -log[10]("adjusted p-value")
    ),

    color = "Classificação"
  ) +

  plot_theme()


# ------------------------------------------------------------------------------
# 16.5 Salvar figura
# ------------------------------------------------------------------------------

ggsave(
  filename = volcano_plot_file,
  plot = p_volcano,
  width = 9,
  height = 7,
  dpi = 300
)


# ==============================================================================
# 17. HEATMAP DOS ATÉ 50 DEGs PRINCIPAIS
# ==============================================================================
#
# OBJETIVO
# ------------------------------------------------------------------------------
# Gerar o heatmap do GSE125583 com o MESMO PADRÃO VISUAL utilizado nos
# heatmaps de HD e PD:
#
#   - pheatmap;
#   - top 50 DEGs principais;
#   - matriz VST;
#   - scale = "row" diretamente dentro do pheatmap;
#   - clustering hierárquico padrão de genes e amostras;
#   - paleta padrão do pheatmap;
#   - símbolos gênicos como rótulos das linhas;
#   - anotação superior por grupo;
#   - proporção gráfica semelhante à utilizada em HD/PD.
#
# IMPORTANTE
# ------------------------------------------------------------------------------
# O heatmap continua sendo exclusivamente exploratório.
# O DESeq2 foi ajustado sobre RAW COUNTS. VST e scale = "row" são usados
# somente para visualização.
# ==============================================================================

heatmap_created <- FALSE

n_heatmap_genes <- min(
  50L,
  nrow(
    primary_deg
  )
)

if (n_heatmap_genes > 1) {

  # ---------------------------------------------------------------------------
  # 17.1 Selecionar os top 50 DEGs principais
  # ---------------------------------------------------------------------------
  #
  # primary_deg já está ordenado por:
  #   1) menor padj;
  #   2) maior |log2FoldChange| em caso de necessidade de desempate.

  top_genes <- primary_deg[
    seq_len(
      n_heatmap_genes
    ),
    GeneID
  ]


  # ---------------------------------------------------------------------------
  # 17.2 VST para visualização
  # ---------------------------------------------------------------------------

  vsd_main <- tryCatch(
    {
      vst(
        dds_main,
        blind = FALSE
      )
    },
    error = function(e) {

      message(
        "vst() não pôde ser utilizado diretamente; ",
        "usando varianceStabilizingTransformation()."
      )

      varianceStabilizingTransformation(
        dds_main,
        blind = FALSE
      )
    }
  )

  vst_matrix <- assay(
    vsd_main
  )

  missing_top_genes <- setdiff(
    top_genes,
    rownames(
      vst_matrix
    )
  )

  if (length(missing_top_genes) > 0) {
    stop(
      "Genes selecionados para heatmap ausentes na matriz VST: ",
      paste(
        missing_top_genes,
        collapse = ", "
      )
    )
  }

  mat_heatmap <- vst_matrix[
    top_genes,
    ,
    drop = FALSE
  ]


  # ---------------------------------------------------------------------------
  # 17.3 Usar símbolos gênicos como labels
  # ---------------------------------------------------------------------------
  #
  # Os heatmaps de HD/PD exibem essencialmente o símbolo do gene.
  # GeneID só é usado como fallback se Symbol estiver ausente.

  heatmap_labels <- primary_deg[
    match(
      top_genes,
      GeneID
    ),
    .(
      GeneID,
      Symbol
    )
  ]

  if (!identical(
    heatmap_labels$GeneID,
    top_genes
  )) {
    stop(
      "Falha ao alinhar os labels dos genes do heatmap."
    )
  }

  gene_display <- fifelse(
    !is.na(
      heatmap_labels$Symbol
    ) &
      trimws(
        heatmap_labels$Symbol
      ) != "",
    heatmap_labels$Symbol,
    heatmap_labels$GeneID
  )

  # make.unique evita nomes de linhas duplicados caso símbolos se repitam.
  rownames(
    mat_heatmap
  ) <- make.unique(
    gene_display
  )


  # ---------------------------------------------------------------------------
  # 17.4 Anotação das colunas
  # ---------------------------------------------------------------------------

  annotation_col <- data.frame(
    Group = factor(
      metadata$group,
      levels = GROUP_LEVELS
    )
  )

  rownames(
    annotation_col
  ) <- metadata$sample_id

  annotation_col <- annotation_col[
    colnames(
      mat_heatmap
    ),
    ,
    drop = FALSE
  ]

  if (!identical(
    rownames(
      annotation_col
    ),
    colnames(
      mat_heatmap
    )
  )) {
    stop(
      "Falha ao alinhar annotation_col e matriz do heatmap."
    )
  }

  # Mantemos as cores de grupo adotadas no projeto:
  # Control = cinza; AD = laranja.
  annotation_colors <- list(
    Group = c(
      "Control" = "#7A7A7A",
      "AD" = "#E67E22"
    )
  )


  # ---------------------------------------------------------------------------
  # 17.5 Gerar heatmap no padrão HD/PD
  # ---------------------------------------------------------------------------
  #
  # Diferenças deliberadamente removidas em relação à versão anterior:
  #
  #   - não pré-calculamos z-score manualmente;
  #   - não definimos paleta vermelho/branco/azul customizada;
  #   - não usamos "SYMBOL [GeneID]" como label;
  #   - não usamos figura larga de 13 x 10.
  #
  # O pheatmap recebe a matriz VST e executa scale = "row", reproduzindo
  # o estilo de geração usado nos heatmaps anteriores.
  #
  # cluster_rows e cluster_cols usam os defaults do pheatmap:
  #   distância euclidiana + complete linkage.
  #
  # Com 289 amostras, as colunas do GSE125583 serão inevitavelmente mais
  # estreitas que em HD/PD (que têm bem menos amostras). Isso é uma diferença
  # real de tamanho amostral e não deve ser "corrigido" removendo indivíduos.

  pheatmap::pheatmap(
    mat_heatmap,

    scale = "row",

    annotation_col = annotation_col,
    annotation_colors = annotation_colors,

    cluster_rows = TRUE,
    cluster_cols = TRUE,

    show_colnames = FALSE,
    show_rownames = TRUE,

    border_color = NA,

    fontsize = 10,
    fontsize_row = 8,

    main = "Top 50 DEGs - AD",

    filename = heatmap_plot_file,

    width = 8,
    height = 9
  )

  heatmap_created <- file.exists(
    heatmap_plot_file
  )
}


# ==============================================================================
# 18. PREPARAR MODELO AJUSTADO ~ age + sex + group
# ==============================================================================

adjusted_complete <- (
  !is.na(
    metadata$age
  ) &
    !is.na(
      metadata$sex
    )
)

metadata_adjusted <- metadata[
  adjusted_complete
]

excluded_adjusted <- metadata[
  !adjusted_complete,
  .(
    sample_id,
    title,
    case_id,
    group,
    age,
    sex,
    braak_score
  )
]

n_adjusted_samples <- nrow(
  metadata_adjusted
)

n_adjusted_ad <- metadata_adjusted[
  group == "AD",
  .N
]

n_adjusted_control <- metadata_adjusted[
  group == "Control",
  .N
]

# Criar idade centralizada.
# Isso não muda a comparação AD vs Control; apenas melhora a parametrização
# numérica e a interpretação do intercepto.
metadata_adjusted[
  ,
  age_centered := age - mean(
    age
  )
]

metadata_adjusted[
  ,
  group := factor(
    group,
    levels = GROUP_LEVELS
  )
]

metadata_adjusted[
  ,
  sex := factor(
    sex,
    levels = c(
      "F",
      "M"
    )
  )
]


# ==============================================================================
# 19. VALIDAR MATRIZ DE DESENHO DO MODELO AJUSTADO
# ==============================================================================

design_matrix_adjusted <- model.matrix(
  ~ age_centered + sex + group,
  data = metadata_adjusted
)

if (
  qr(
    design_matrix_adjusted
  )$rank <
    ncol(
      design_matrix_adjusted
    )
) {
  stop(
    "O modelo ajustado ~ age_centered + sex + group não possui posto completo."
  )
}


# ==============================================================================
# 20. EXECUTAR MODELO AJUSTADO
# ==============================================================================

adjusted_sample_ids <- metadata_adjusted$sample_id

adjusted_count_matrix <- count_matrix[
  ,
  adjusted_sample_ids,
  drop = FALSE
]

col_data_adjusted <- as.data.frame(
  metadata_adjusted
)

rownames(
  col_data_adjusted
) <- col_data_adjusted$sample_id

if (!identical(
  colnames(
    adjusted_count_matrix
  ),
  rownames(
    col_data_adjusted
  )
)) {
  stop(
    "Counts e colData do modelo ajustado não estão alinhados."
  )
}

dds_adjusted <- DESeqDataSetFromMatrix(
  countData = adjusted_count_matrix,
  colData = col_data_adjusted,
  design = ~ age_centered + sex + group
)

dds_adjusted <- DESeq(
  dds_adjusted,
  quiet = TRUE
)

adjusted_results_raw <- results(
  dds_adjusted,
  contrast = c(
    "group",
    "AD",
    "Control"
  ),
  alpha = ALPHA,
  independentFiltering = TRUE
)

adjusted_results <- results_to_dt(
  adjusted_results_raw,
  model_name = "adjusted_age_sex_group"
)

adjusted_results <- annotate_results(
  adjusted_results,
  gene_master
)

fwrite(
  adjusted_results,
  adjusted_results_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 21. COMPARAR MODELO PRINCIPAL E AJUSTADO
# ==============================================================================

main_for_comparison <- main_results[
  ,
  .(
    GeneID,
    Symbol,
    main_log2FoldChange = log2FoldChange,
    main_pvalue = pvalue,
    main_padj = padj,
    main_primary = is_primary,
    main_sensitivity = is_sensitivity
  )
]

adjusted_for_comparison <- adjusted_results[
  ,
  .(
    GeneID,
    adjusted_log2FoldChange = log2FoldChange,
    adjusted_pvalue = pvalue,
    adjusted_padj = padj,
    adjusted_primary = is_primary,
    adjusted_sensitivity = is_sensitivity
  )
]

model_comparison <- merge(
  main_for_comparison,
  adjusted_for_comparison,
  by = "GeneID",
  all = FALSE,
  sort = FALSE
)

model_comparison[
  ,
  same_direction := (
    !is.na(
      main_log2FoldChange
    ) &
      !is.na(
        adjusted_log2FoldChange
      ) &
      sign(
        main_log2FoldChange
      ) ==
      sign(
        adjusted_log2FoldChange
      )
  )
]

model_comparison[
  ,
  primary_status := fcase(
    main_primary & adjusted_primary, "Both",
    main_primary & !adjusted_primary, "Main_only",
    !main_primary & adjusted_primary, "Adjusted_only",
    default = "Neither"
  )
]

fwrite(
  model_comparison,
  model_comparison_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 22. MÉTRICAS DE ESTABILIDADE ENTRE MODELOS
# ==============================================================================

valid_lfc_comparison <- model_comparison[
  !is.na(
    main_log2FoldChange
  ) &
    !is.na(
      adjusted_log2FoldChange
    )
]

rho_main_adjusted <- suppressWarnings(
  cor(
    valid_lfc_comparison$main_log2FoldChange,
    valid_lfc_comparison$adjusted_log2FoldChange,
    method = "spearman"
  )
)

n_primary_main <- sum(
  model_comparison$main_primary,
  na.rm = TRUE
)

n_primary_adjusted <- sum(
  model_comparison$adjusted_primary,
  na.rm = TRUE
)

n_primary_both <- sum(
  model_comparison$main_primary &
    model_comparison$adjusted_primary,
  na.rm = TRUE
)

n_primary_main_only <- sum(
  model_comparison$main_primary &
    !model_comparison$adjusted_primary,
  na.rm = TRUE
)

n_primary_adjusted_only <- sum(
  !model_comparison$main_primary &
    model_comparison$adjusted_primary,
  na.rm = TRUE
)

shared_primary <- model_comparison[
  main_primary &
    adjusted_primary
]

n_shared_primary_same_direction <- sum(
  shared_primary$same_direction,
  na.rm = TRUE
)

n_shared_primary_opposite_direction <- (
  nrow(
    shared_primary
  ) -
    n_shared_primary_same_direction
)


# ==============================================================================
# 23. PLOT — log2FC PRINCIPAL vs AJUSTADO
# ==============================================================================

comparison_plot_data <- copy(
  valid_lfc_comparison
)

comparison_plot_data[
  ,
  plot_class := factor(
    primary_status,
    levels = c(
      "Neither",
      "Main_only",
      "Adjusted_only",
      "Both"
    )
  )
]

comparison_colors <- c(
  "Neither" = "#BDBDBD",
  "Main_only" = "#E69F00",
  "Adjusted_only" = "#56B4E9",
  "Both" = "#009E73"
)

p_model_comparison <- ggplot(
  comparison_plot_data,
  aes(
    x = main_log2FoldChange,
    y = adjusted_log2FoldChange,
    color = plot_class
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.2
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
    color = "grey40",
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = comparison_colors
  ) +
  labs(
    title = "Estabilidade do efeito de AD — GSE125583",
    subtitle = paste0(
      "Modelo principal vs ajustado; Spearman rho = ",
      format_num(
        rho_main_adjusted,
        digits = 3
      )
    ),
    x = "log2FC — modelo principal (~ group)",
    y = "log2FC — ajustado (~ age + sex + group)",
    color = "DEG principal"
  ) +
  plot_theme()

ggsave(
  filename = model_comparison_plot_file,
  plot = p_model_comparison,
  width = 8,
  height = 7,
  dpi = 300
)


# ==============================================================================
# 24. CRIAR RESUMO DOS DOIS MODELOS
# ==============================================================================

main_summary <- summarize_results(
  main_results,
  model_label = "main_group",
  n_samples = EXPECTED_N_SAMPLES
)

adjusted_summary <- summarize_results(
  adjusted_results,
  model_label = "adjusted_age_sex_group",
  n_samples = n_adjusted_samples
)

deseq2_summary <- rbind(
  main_summary,
  adjusted_summary,
  fill = TRUE
)

fwrite(
  deseq2_summary,
  summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 25. RELATÓRIO TEXTUAL
# ==============================================================================

sink(
  audit_file
)

cat("============================================================\n")
cat("05 — DESeq2 GSE125583\n")
cat("============================================================\n\n")


cat("------------------------------------------------------------\n")
cat("1. ENTRADA\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes após pré-filtro: ",
  format_int(
    n_genes_tested
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras totais: ",
  format_int(
    EXPECTED_N_SAMPLES
  ),
  "\n",
  sep = ""
)

cat(
  "AD: ",
  format_int(
    n_ad
  ),
  "\n",
  sep = ""
)

cat(
  "Control: ",
  format_int(
    n_control
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("2. MODELO PRINCIPAL\n")
cat("------------------------------------------------------------\n")

cat(
  "Design: ~ group\n"
)

cat(
  "Referência: Control\n"
)

cat(
  "Contraste: AD vs Control\n"
)

cat(
  "Independent filtering: TRUE\n"
)

cat(
  "alpha: ",
  ALPHA,
  "\n",
  sep = ""
)

cat(
  "lfcThreshold dentro de results(): NÃO utilizado\n\n"
)


cat("------------------------------------------------------------\n")
cat("3. SIZE FACTORS — MODELO PRINCIPAL\n")
cat("------------------------------------------------------------\n")

cat(
  "Mínimo: ",
  format_num(
    min(
      main_size_factors
    ),
    digits = 3
  ),
  "\n",
  sep = ""
)

cat(
  "Mediana: ",
  format_num(
    median(
      main_size_factors
    ),
    digits = 3
  ),
  "\n",
  sep = ""
)

cat(
  "Máximo: ",
  format_num(
    max(
      main_size_factors
    ),
    digits = 3
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("4. RESULTADOS — MODELO PRINCIPAL\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes na tabela de resultados: ",
  format_int(
    nrow(
      main_results
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes com p-value não NA: ",
  format_int(
    sum(
      !is.na(
        main_results$pvalue
      )
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes com padj não NA: ",
  format_int(
    sum(
      !is.na(
        main_results$padj
      )
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Menor p-value: ",
  format_p(
    min(
      main_results$pvalue,
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Menor padj: ",
  format_p(
    min(
      main_results$padj,
      na.rm = TRUE
    )
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("5. CLASSIFICAÇÕES — MODELO PRINCIPAL\n")
cat("------------------------------------------------------------\n")

cat(
  "FDR < 0,05 (independentemente de log2FC): ",
  format_int(
    sum(
      main_results$fdr_significant,
      na.rm = TRUE
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Critério principal: FDR < 0,05 e |log2FC| > ",
  PRIMARY_LFC,
  "\n",
  sep = ""
)

cat(
  "  Total: ",
  format_int(
    nrow(
      primary_deg
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  Up: ",
  format_int(
    sum(
      primary_deg$primary_direction == "Up",
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  Down: ",
  format_int(
    sum(
      primary_deg$primary_direction == "Down",
      na.rm = TRUE
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Sensibilidade: FDR < 0,05 e |log2FC| > ",
  SENSITIVITY_LFC,
  "\n",
  sep = ""
)

cat(
  "  Total: ",
  format_int(
    nrow(
      sensitivity_deg
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  Up: ",
  format_int(
    sum(
      sensitivity_deg$sensitivity_direction == "Up",
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  Down: ",
  format_int(
    sum(
      sensitivity_deg$sensitivity_direction == "Down",
      na.rm = TRUE
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Critério nominal comparativo: p < 0,05 e |log2FC| > ",
  SENSITIVITY_LFC,
  "\n",
  sep = ""
)

cat(
  "  Total: ",
  format_int(
    nrow(
      article_like_deg
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  Up: ",
  format_int(
    sum(
      article_like_deg$article_like_direction == "Up",
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  Down: ",
  format_int(
    sum(
      article_like_deg$article_like_direction == "Down",
      na.rm = TRUE
    )
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("6. HEATMAP\n")
cat("------------------------------------------------------------\n")

cat(
  "Heatmap gerado: ",
  heatmap_created,
  "\n",
  sep = ""
)

cat(
  "Genes no heatmap: ",
  format_int(
    n_heatmap_genes
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("7. MODELO AJUSTADO — SENSIBILIDADE A COVARIÁVEIS\n")
cat("------------------------------------------------------------\n")

cat(
  "Design: ~ age_centered + sex + group\n"
)

cat(
  "Amostras utilizadas: ",
  format_int(
    n_adjusted_samples
  ),
  "\n",
  sep = ""
)

cat(
  "AD: ",
  format_int(
    n_adjusted_ad
  ),
  "\n",
  sep = ""
)

cat(
  "Control: ",
  format_int(
    n_adjusted_control
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras excluídas por age/sex ausente: ",
  format_int(
    nrow(
      excluded_adjusted
    )
  ),
  "\n",
  sep = ""
)

if (nrow(
  excluded_adjusted
) > 0) {

  cat("\nAmostras excluídas do modelo ajustado:\n")
  print(
    excluded_adjusted
  )
}

cat("\n")

cat(
  "DEGs principais ajustados: ",
  format_int(
    sum(
      adjusted_results$is_primary,
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs sensibilidade ajustados: ",
  format_int(
    sum(
      adjusted_results$is_sensitivity,
      na.rm = TRUE
    )
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("8. ESTABILIDADE PRINCIPAL vs AJUSTADO\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes comparados: ",
  format_int(
    nrow(
      valid_lfc_comparison
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Spearman rho dos log2FC: ",
  format_num(
    rho_main_adjusted,
    digits = 4
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs principais — modelo principal: ",
  format_int(
    n_primary_main
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs principais — modelo ajustado: ",
  format_int(
    n_primary_adjusted
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs principais em ambos: ",
  format_int(
    n_primary_both
  ),
  "\n",
  sep = ""
)

cat(
  "Somente modelo principal: ",
  format_int(
    n_primary_main_only
  ),
  "\n",
  sep = ""
)

cat(
  "Somente modelo ajustado: ",
  format_int(
    n_primary_adjusted_only
  ),
  "\n",
  sep = ""
)

cat(
  "Compartilhados com mesma direção: ",
  format_int(
    n_shared_primary_same_direction
  ),
  "\n",
  sep = ""
)

cat(
  "Compartilhados com direção oposta: ",
  format_int(
    n_shared_primary_opposite_direction
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("9. INTERPRETAÇÃO METODOLÓGICA\n")
cat("------------------------------------------------------------\n")

cat(
  "- O modelo principal utiliza todas as amostras e design ~ group para manter comparabilidade com HD e PD.\n"
)

cat(
  "- Control é a referência; log2FoldChange positivo indica maior expressão em AD.\n"
)

cat(
  "- O critério principal é aplicado após a estimação do DESeq2, sem lfcThreshold no teste.\n"
)

cat(
  "- O critério nominal p < 0,05 + |log2FC| > 1 é apenas comparativo e não substitui FDR.\n"
)

cat(
  "- O modelo ajustado testa a estabilidade do efeito de AD após ajuste por idade e sexo.\n"
)

cat(
  "- Braak não foi incluído como covariável porque representa uma característica neuropatológica fortemente associada ao processo da doença.\n"
)

cat(
  "- O heatmap, quando disponível, utiliza VST apenas para visualização; os testes usam raw counts.\n"
)

cat(
  "- Nenhuma integração funcional ou transdiagnóstica é realizada nesta etapa.\n"
)

cat("\n")
cat("============================================================\n")
cat("05_deseq2_GSE125583.R CONCLUÍDO\n")
cat("============================================================\n")

sink()


# ==============================================================================
# 26. RESUMO NO TERMINAL
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("05 — DESeq2 GSE125583 CONCLUÍDO\n")
cat("============================================================\n")

cat(
  "Genes testados: ",
  format_int(
    n_genes_tested
  ),
  "\n",
  sep = ""
)

cat(
  "FDR < 0,05: ",
  format_int(
    sum(
      main_results$fdr_significant,
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs principais: ",
  format_int(
    nrow(
      primary_deg
    )
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs sensibilidade: ",
  format_int(
    nrow(
      sensitivity_deg
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes article-like: ",
  format_int(
    nrow(
      article_like_deg
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras modelo ajustado: ",
  format_int(
    n_adjusted_samples
  ),
  "\n",
  sep = ""
)

cat(
  "Spearman main vs ajustado: ",
  format_num(
    rho_main_adjusted,
    digits = 3
  ),
  "\n",
  sep = ""
)

cat(
  "Resultados principais: ",
  main_results_file,
  "\n",
  sep = ""
)

cat(
  "Resumo: ",
  summary_file,
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
