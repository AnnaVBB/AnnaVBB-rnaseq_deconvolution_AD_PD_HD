#!/usr/bin/env Rscript

# ==============================================================================
# 04_qc_counts.R
#
# OBJETIVOS
# ------------------------------------------------------------------------------
# 1. Avaliar a qualidade das matrizes de contagens preparadas no 02.
# 2. Resumir características das bibliotecas/amostras antes da filtragem.
# 3. Avaliar descritivamente thresholds de 5, 10 e 20 contagens.
# 4. Aplicar o pré-filtro definido para EXPRESSÃO DIFERENCIAL: count >= 10 em pelo menos o tamanho do menor grupo experimental
#       AD = 8 amostras
#       HD = 20 amostras
#       PD = 29 amostras
# 5. Gerar matrizes pré-filtradas para 05_deseq2.R.
# 6. Realizar VST + PCA exclusivamente para QC exploratório.
# 7. Calcular a interseção gênica após a pré-filtragem para DE.
#
# ENTRADAS
# ------------------------------------------------------------------------------
# results/counts/counts_AD.tsv
# results/counts/counts_HD.tsv
# results/counts/counts_PD.tsv
# results/metadata/metadata_AD.tsv
# results/metadata/metadata_HD.tsv
# results/metadata/metadata_PD.tsv
# results/harmonization/gene_master.tsv
#
# SAÍDAS
# ------------------------------------------------------------------------------
# results/qc/
#   counts_AD_DE_filtered.tsv
#   counts_HD_DE_filtered.tsv
#   counts_PD_DE_filtered.tsv
#   qc_gene_filter_summary.tsv
#   qc_sample_summary.tsv
#
# results/plots/
#   AD/PCA.png
#   HD/PCA.png
#   PD/PCA.png
#
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
})

# ==============================================================================
# PARÂMETROS

# Critério principal aprovado para expressão diferencial.
min_count <- 10

# Número mínimo de amostras por dataset. Corresponde ao tamanho do menor grupo experimental.
min_samples <- c(
  AD = 8,
  HD = 20,
  PD = 29
)

# Thresholds avaliados apenas descritivamente. Não modificam o critério principal.
sensitivity_counts <- c(
  5,
  10,
  20
)


# CAMINHOS
counts_files <- c(
  AD = "results/counts/counts_AD.tsv",
  HD = "results/counts/counts_HD.tsv",
  PD = "results/counts/counts_PD.tsv"
)

metadata_files <- c(
  AD = "results/metadata/metadata_AD.tsv",
  HD = "results/metadata/metadata_HD.tsv",
  PD = "results/metadata/metadata_PD.tsv"
)

gene_master_file <- "results/harmonization/gene_master.tsv"

# Diretório para tabelas e matrizes de QC
output_dir <- "results/qc"

# Diretório geral para figuras.
plots_dir <- "results/plots"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Criar uma pasta de figuras para cada doença.
for (dataset in c("AD", "HD", "PD")) {
  dir.create(
    file.path(
      plots_dir,
      dataset
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# VERIFICAR ARQUIVOS DE ENTRADA
input_files <- c(
  counts_files,
  metadata_files,
  gene_master_file
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Arquivos necessários não encontrados:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}

# ==============================================================================
# FUNÇÃO DE QC E PRÉ-FILTRAGEM

run_qc <- function(
  dataset_name,
  disease_level,
  counts_file,
  metadata_file
) {

  cat("\n")
  cat("============================================\n")
  cat("QC -", dataset_name, "\n")
  cat("============================================\n")

  # ============================================================================
  # Ler matriz de contagens
  # ============================================================================
  counts_dt <- fread(
    counts_file,
    check.names = FALSE
  )

  if (!"GeneID" %in% names(counts_dt)) {
    stop(
      dataset_name,
      ": coluna GeneID não encontrada."
    )
  }

  gene_ids <- as.character(
    counts_dt$GeneID
  )

  if (anyNA(gene_ids)) {
    stop(
      dataset_name,
      ": existem GeneIDs ausentes."
    )
  }

  if (anyDuplicated(gene_ids)) {
    stop(
      dataset_name,
      ": existem GeneIDs duplicados."
    )
  }

  sample_columns <- setdiff(
    names(counts_dt),
    "GeneID"
  )


  count_matrix <- as.matrix(
    counts_dt[
      ,
      ..sample_columns
    ]
  )

  storage.mode(count_matrix) <- "numeric"
  rownames(count_matrix) <- gene_ids


  # Validar matriz de contagens
  if (anyNA(count_matrix)) {

    stop(
      dataset_name,
      ": existem NA na matriz de contagens."
    )
  }

  if (any(count_matrix < 0)) {
    stop(
      dataset_name,
      ": existem contagens negativas."
    )
  }

  n_fractional <- sum(
    abs(
      count_matrix -
        round(count_matrix)
    ) > 1e-8
  )


  if (n_fractional > 0) {
    stop(
      dataset_name,
      ": foram encontrados ",
      n_fractional,
      " valores não inteiros."
    )
  }

  # A integralidade já foi validada. Portanto, esta conversão NÃO realiza novo arredondamento.
  storage.mode(count_matrix) <- "integer"

  # Ler metadata
    metadata <- fread(
    metadata_file,
    na.strings = c(
      "NA",
      ""
    )
  )

  required_metadata <- c(
    "sample_id",
    "group"
  )

  missing_metadata <- setdiff(
    required_metadata,
    names(metadata)
  )

  if (length(missing_metadata) > 0) {
    stop(
      dataset_name,
      ": campos ausentes no metadata: ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  }

  if (anyNA(metadata$sample_id)) {
    stop(
      dataset_name,
      ": existem sample_id ausentes."
    )
  }

  if (anyDuplicated(metadata$sample_id)) {
    stop(
      dataset_name,
      ": existem sample_id duplicados."
    )
  }

  if (anyNA(metadata$group)) {
    stop(
      dataset_name,
      ": existem grupos ausentes."
    )
  }

  # Validar counts <-> metadata
  if (
    !setequal(
      colnames(count_matrix),
      metadata$sample_id
    )
  ) {
    stop(
      dataset_name,
      ": counts e metadata possuem amostras diferentes."
    )
  }


  # Reordenar metadata explicitamente para acompanhar a matriz.
  metadata <- metadata[
    match(
      colnames(count_matrix),
      sample_id
    )
  ]

  if (
    !identical(
      metadata$sample_id,
      colnames(count_matrix)
    )
  ) {
    stop(
      dataset_name,
      ": falha ao alinhar metadata e counts."
    )
  }


  # Validar grupos experimentais
  metadata[
    ,
    group := factor(
      group,
      levels = c(
        "Control",
        disease_level
      )
    )
  ]

  if (anyNA(metadata$group)) {
    stop(
      dataset_name,
      ": grupo inesperado encontrado."
    )
  }

  observed_groups <- table(
    metadata$group
  )

  cat("\nAmostras por grupo:\n")

  print(
    observed_groups
  )

  observed_min_group <- min(
    observed_groups
  )

  if (
    observed_min_group !=
      min_samples[[dataset_name]]
  ) {
    stop(
      dataset_name,
      ": tamanho do menor grupo diferente do esperado. ",
      "Observado = ",
      observed_min_group,
      "; esperado = ",
      min_samples[[dataset_name]],
      "."
    )
  }

  
  # QC DAS AMOSTRAS ANTES DA FILTRAGEM
  library_size <- colSums(
    count_matrix
  )

  genes_nonzero <- colSums(
    count_matrix > 0
  )

  zero_fraction <- colMeans(
    count_matrix == 0
  )

  sample_summary <- data.table(
    dataset = dataset_name,
    sample_id = colnames(count_matrix),
    group = as.character(metadata$group),
    library_size = as.numeric(library_size),
    genes_nonzero = as.integer(genes_nonzero),
    zero_fraction = as.numeric(zero_fraction)
  )

  cat("\nQC das bibliotecas:\n")

  cat(
    "Library size mínima:",
    min(library_size),
    "\n"
  )

  cat(
    "Library size mediana:",
    median(library_size),
    "\n"
  )

  cat(
    "Library size máxima:",
    max(library_size),
    "\n"
  )

  cat(
    "Genes detectados/amostra - mínimo:",
    min(genes_nonzero),
    "\n"
  )

  cat(
    "Genes detectados/amostra - mediana:",
    median(genes_nonzero),
    "\n"
  )

  cat(
    "Genes detectados/amostra - máximo:",
    max(genes_nonzero),
    "\n"
  )

  cat(
    "Fração de zeros/amostra - mínima:",
    min(zero_fraction),
    "\n"
  )

  cat(
    "Fração de zeros/amostra - mediana:",
    median(zero_fraction),
    "\n"
  )

  cat(
    "Fração de zeros/amostra - máxima:",
    max(zero_fraction),
    "\n"
  )

  # GENES ALL-ZERO
  all_zero <- rowSums(
    count_matrix
  ) == 0

  n_all_zero <- sum(
    all_zero
  )

  cat(
    "\nGenes all-zero:",
    n_all_zero,
    "\n"
  )


  # ============================================================================
  # ANÁLISE DESCRITIVA DE SENSIBILIDADE
  # Mantemos o número mínimo de amostras definido pelo menor grupo e variamos apenas o número mínimo de contagens.

  sensitivity_dt <- rbindlist(
    lapply(
      sensitivity_counts,
      function(threshold) {

        keep_sensitivity <- rowSums(
          count_matrix >= threshold
        ) >= min_samples[[dataset_name]]

        data.table(
          dataset = dataset_name,
          min_count = threshold,
          min_samples = min_samples[[dataset_name]],
          genes_before = nrow(count_matrix),
          genes_retained = sum(keep_sensitivity),
          genes_removed = sum(!keep_sensitivity),
          percent_retained = (
            100 *
              sum(keep_sensitivity) /
              nrow(count_matrix)
          )
        )
      }
    )
  )

  cat("\nSensibilidade da pré-filtragem:\n")

  print(
    sensitivity_dt
  )

  # APLICAR PRÉ-FILTRO PRINCIPAL PARA EXPRESSÃO DIFERENCIAL
  keep <- rowSums(
    count_matrix >= min_count
  ) >= min_samples[[dataset_name]]

  filtered_matrix <- count_matrix[
    keep,
    ,
    drop = FALSE
  ]

  cat("\n")
  cat("--------------------------------------------\n")
  cat("PRÉ-FILTRAGEM PARA EXPRESSÃO DIFERENCIAL\n")
  cat("--------------------------------------------\n")

  cat(
    "Critério: count >=",
    min_count,
    "em pelo menos",
    min_samples[[dataset_name]],
    "amostras\n"
  )

  cat(
    "Genes antes:",
    nrow(count_matrix),
    "\n"
  )

  cat(
    "Genes mantidos:",
    nrow(filtered_matrix),
    "\n"
  )

  cat(
    "Genes removidos:",
    sum(!keep),
    "\n"
  )

  cat(
    "Percentual mantido:",
    round(
      100 *
        nrow(filtered_matrix) /
        nrow(count_matrix),
      2
    ),
    "%\n"
  )

  # VALIDAR MATRIZ FILTRADA
  if (nrow(filtered_matrix) == 0) {

    stop(
      dataset_name,
      ": nenhum gene permaneceu após a pré-filtragem."
    )
  }

  if (anyNA(filtered_matrix)) {
    stop(
      dataset_name,
      ": NA encontrado após pré-filtragem."
    )
  }


  if (any(filtered_matrix < 0)) {
    stop(
      dataset_name,
      ": valor negativo encontrado após pré-filtragem."
    )
  }


  if (
    any(
      abs(
        filtered_matrix -
          round(filtered_matrix)
      ) > 1e-8
    )
  ) {

    stop(
      dataset_name,
      ": valor não inteiro encontrado após pré-filtragem."
    )
  }

  # SALVAR MATRIZ FILTRADA PARA EXPRESSÃO DIFERENCIAL
  
  filtered_dt <- as.data.table(
    filtered_matrix,
    keep.rownames = "GeneID"
  )

  filtered_output <- file.path(
    output_dir,
    paste0(
      "counts_",
      dataset_name,
      "_DE_filtered.tsv"
    )
  )

  fwrite(
    filtered_dt,
    filtered_output,
    sep = "\t",
    quote = FALSE
  )

  # ============================================================================
  # VST PARA QC / PCA
  # IMPORTANTE: O VST NÃO será utilizado como entrada da análise diferencial. blind = TRUE é utilizado porque esta é uma transformação exploratória destinada ao QC global das amostras.
  # ============================================================================
  cat(
    "\nExecutando VST para PCA...\n"
  )

  coldata <- data.frame(
    group = metadata$group,
    row.names = metadata$sample_id
  )

  dds_qc <- DESeqDataSetFromMatrix(
    countData = filtered_matrix,
    colData = coldata,
    design = ~ group
  )

  vsd <- vst(
    dds_qc,
    blind = TRUE
  )

  # PCA
  pca_data <- plotPCA(
    vsd,
    intgroup = "group",
    returnData = TRUE
  )

  percent_var <- round(
    100 *
      attr(
        pca_data,
        "percentVar"
      )
  )

  pca_data$sample_id <- rownames(
    pca_data
  )
    group_colors <- switch(
    dataset_name,
    AD = c(
        "Control" = "#4D4D4D",
        "AD" = "#723564"
    ),
    HD = c(
        "Control" = "#4D4D4D",
        "HD" = "#004aad"
    ),
    PD = c(
        "Control" = "#4D4D4D",
        "PD" = "#ad4931"
    )
    )
 
    p <- ggplot(
    pca_data,
    aes(
        x = PC1,
        y = PC2,
        color = group,
        shape = group
    )
    ) +
    geom_point(
        size = 3
    ) +
    scale_color_manual(
        values = group_colors
    ) +
    labs(
        title = paste0(
        "PCA - ",
        dataset_name
        ),
        x = paste0(
        "PC1: ",
        percent_var[1],
        "% variance"
        ),
        y = paste0(
        "PC2: ",
        percent_var[2],
        "% variance"
        ),
        color = "Group",
        shape = "Group"
    ) +
    theme_classic()

    pca_output <- file.path(
        plots_dir,
        dataset_name,
        "PCA.png"
    )

    ggsave(
        filename = pca_output,
        plot = p,
        width = 7,
        height = 5,
        dpi = 300
    )

    cat(
        "PCA salva em:",
        pca_output,
        "\n"
    )

  # RESUMO DO FILTRO PRINCIPAL
  
  filter_summary <- data.table(
    dataset = dataset_name,
    n_samples = ncol(count_matrix),
    min_group_size = min_samples[[dataset_name]],
    min_count = min_count,
    genes_before = nrow(count_matrix),
    genes_all_zero = n_all_zero,
    genes_retained = nrow(filtered_matrix),
    genes_removed = sum(!keep),
    percent_retained = (
      100 *
        nrow(filtered_matrix) /
        nrow(count_matrix)
    )
  )

  # RETORNAR RESULTADOS
  return(
    list(
      sample_summary = sample_summary,
      sensitivity = sensitivity_dt,
      filter_summary = filter_summary,
      retained_genes = rownames(filtered_matrix)
    )
  )
}

# EXECUTAR AD
qc_AD <- run_qc(
  dataset_name = "AD",
  disease_level = "AD",
  counts_file = counts_files[["AD"]],
  metadata_file = metadata_files[["AD"]]
)

# EXECUTAR HD
qc_HD <- run_qc(
  dataset_name = "HD",
  disease_level = "HD",
  counts_file = counts_files[["HD"]],
  metadata_file = metadata_files[["HD"]]
)

# EXECUTAR PD
qc_PD <- run_qc(
  dataset_name = "PD",
  disease_level = "PD",
  counts_file = counts_files[["PD"]],
  metadata_file = metadata_files[["PD"]]
)

# CONSOLIDAR QC DAS AMOSTRAS
sample_summary_all <- rbindlist(
  list(
    qc_AD$sample_summary,
    qc_HD$sample_summary,
    qc_PD$sample_summary
  )
)

fwrite(
  sample_summary_all,
  file.path(
    output_dir,
    "qc_sample_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

# CONSOLIDAR RESULTADOS DO PRÉ-FILTRO
filter_summary_all <- rbindlist(
  list(
    qc_AD$filter_summary,
    qc_HD$filter_summary,
    qc_PD$filter_summary
  )
)

sensitivity_all <- rbindlist(
  list(
    qc_AD$sensitivity,
    qc_HD$sensitivity,
    qc_PD$sensitivity
  )
)

# Linhas referentes ao filtro efetivamente utilizado.
filter_summary_selected <- copy(
  filter_summary_all
)

filter_summary_selected[
  ,
  `:=`(
    analysis_type = "selected",
    min_samples_sensitivity = NA_integer_
  )
]

filter_summary_selected <- filter_summary_selected[
  ,
  .(
    dataset,
    analysis_type,
    n_samples,
    min_group_size,
    min_count,
    min_samples_sensitivity,
    genes_before,
    genes_all_zero,
    genes_retained,
    genes_removed,
    percent_retained
  )
]

# Linhas referentes aos cenários descritivos de sensibilidade.
sensitivity_export <- copy(
  sensitivity_all
)

setnames(
  sensitivity_export,
  old = "min_samples",
  new = "min_samples_sensitivity"
)

sensitivity_export[
  ,
  `:=`(
    analysis_type = "sensitivity",
    n_samples = NA_integer_,
    min_group_size = min_samples_sensitivity,
    genes_all_zero = NA_integer_
  )
]

sensitivity_export <- sensitivity_export[
  ,
  .(
    dataset,
    analysis_type,
    n_samples,
    min_group_size,
    min_count,
    min_samples_sensitivity,
    genes_before,
    genes_all_zero,
    genes_retained,
    genes_removed,
    percent_retained
  )
]

gene_filter_summary <- rbindlist(
  list(
    filter_summary_selected,
    sensitivity_export
  ),
  use.names = TRUE,
  fill = TRUE
)


fwrite(
  gene_filter_summary,
  file.path(
    output_dir,
    "qc_gene_filter_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# UNIVERSO COMUM APÓS PRÉ-FILTRAGEM PARA DE
common_prefiltered <- Reduce(
  intersect,
  list(
    qc_AD$retained_genes,
    qc_HD$retained_genes,
    qc_PD$retained_genes
  )
)

cat("\n")
cat("============================================\n")
cat("UNIVERSO COMUM APÓS PRÉ-FILTRAGEM PARA DE\n")
cat("============================================\n")

cat(
  "Genes mantidos em AD:",
  length(qc_AD$retained_genes),
  "\n"
)

cat(
  "Genes mantidos em HD:",
  length(qc_HD$retained_genes),
  "\n"
)

cat(
  "Genes mantidos em PD:",
  length(qc_PD$retained_genes),
  "\n"
)

cat(
  "Genes presentes após pré-filtragem nos três datasets:",
  length(common_prefiltered),
  "\n"
)


# COMPARAR COM O UNIVERSO COMUM BRUTO DO 03
gene_master <- fread(
  gene_master_file,
  colClasses = list(
    character = "gene_id"
  )
)

required_master_columns <- c(
  "gene_id",
  "in_AD",
  "in_HD",
  "in_PD"
)

missing_master_columns <- setdiff(
  required_master_columns,
  names(gene_master)
)


if (length(missing_master_columns) > 0) {
  stop(
    "gene_master.tsv não contém as colunas esperadas: ",
    paste(
      missing_master_columns,
      collapse = ", "
    )
  )
}

raw_common <- gene_master[
  in_AD == TRUE &
    in_HD == TRUE &
    in_PD == TRUE,
  gene_id
]

cat(
  "Universo comum bruto do 03:",
  length(raw_common),
  "\n"
)


if (length(raw_common) != 18941) {
  stop(
    "O universo comum bruto não corresponde aos 18.941 genes esperados."
  )
}


if (
  !all(
    common_prefiltered %in% raw_common
  )
) {
  stop(
    "Foram encontrados genes no universo comum pós-filtro ",
    "que não pertencem ao universo comum bruto."
  )
}

removed_from_common <- setdiff(
  raw_common,
  common_prefiltered
)

cat(
  "Genes comuns brutos removidos de pelo menos um dataset pelo pré-filtro:",
  length(removed_from_common),
  "\n"
)

cat(
  "Percentual do universo comum bruto preservado:",
  round(
    100 *
      length(common_prefiltered) /
      length(raw_common),
    2
  ),
  "%\n"
)


# RESUMO FINAL
cat("\n")
cat("============================================\n")
cat("RESUMO FINAL DO QC\n")
cat("============================================\n")

print(
  filter_summary_all
)

cat("\nCritério utilizado para expressão diferencial:\n")

cat(
  "count >=",
  min_count,
  "em pelo menos o tamanho do menor grupo de cada dataset.\n"
)

cat(
  "AD:",
  min_samples[["AD"]],
  "amostras\n"
)

cat(
  "HD:",
  min_samples[["HD"]],
  "amostras\n"
)

cat(
  "PD:",
  min_samples[["PD"]],
  "amostras\n"
)

cat("\n")
cat("Matrizes geradas para DESeq2:\n")

cat(
  " - results/qc/counts_AD_DE_filtered.tsv\n"
)

cat(
  " - results/qc/counts_HD_DE_filtered.tsv\n"
)

cat(
  " - results/qc/counts_PD_DE_filtered.tsv\n"
)

cat("\n")
cat("Figuras geradas:\n")

cat(
  " - results/plots/AD/PCA.png\n"
)

cat(
  " - results/plots/HD/PCA.png\n"
)

cat(
  " - results/plots/PD/PCA.png\n"
)

cat("\n")
cat("============================================\n")
cat("04_qc_counts.R CONCLUÍDO\n")
cat("============================================\n")

cat(
  "Nenhuma amostra foi removida automaticamente.\n"
)

cat(
  "VST foi utilizado somente para PCA/QC.\n"
)

cat(
  "As matrizes DE filtradas permanecem em escala de contagens inteiras.\n"
)

cat(
  "O filtro aplicado aqui é específico para expressão diferencial.\n"
)

cat(
  "A deconvolução utilizará posteriormente uma estratégia de filtragem própria.\n"
)