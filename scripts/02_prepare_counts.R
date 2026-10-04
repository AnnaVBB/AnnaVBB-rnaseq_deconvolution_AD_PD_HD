#!/usr/bin/env Rscript

# ==============================================================================
# 02_prepare_counts.R
# Objetivo: Preparar as matrizes de contagens brutas que serão utilizadas nas etapas posteriores do pipeline.
#
# Este script:
#   1. lê as matrizes brutas auditadas em 00_audit_raw_data.R e os metadados preparados em 01_prepare_metadata.R;
#   2. para AD, seleciona exclusivamente GeneID + colunas *_raw, excluindo as colunas *_rpkm da matriz preparada de AD;
#   3. audita e arredonda as contagens fracionárias de AD para o inteiro mais próximo;
#   4. preserva HD e PD sem transformação numérica;
#   4. garante que todas as matrizes tenham GeneID + amostras;
#   5. valida a correspondência e a ordem counts <-> metadata;
#   6. valida que todas as contagens finais sejam inteiras, não negativas e sem valores ausentes;
#
# O arredondamento de AD é realizado porque o arquivo original contém valores fracionários nas colunas identificadas como *_raw, enquanto o fluxo padrão do DESeq2 requer contagens inteiras não negativas.
#
# Saídas:
#   results/counts/counts_AD.tsv
#   results/counts/counts_HD.tsv
#   results/counts/counts_PD.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# Matrizes brutas
ad_file <- "data/raw/GSE53697/GSE53697_RNAseq_AD.txt.gz"

hd_file <- paste0(
  "data/raw/GSE64810/",
  "GSE64810_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)

pd_file <- paste0(
  "data/raw/GSE68719/",
  "GSE68719_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)


# Metadados preparados na etapa 01
ad_metadata_file <- "results/metadata/metadata_AD.tsv"
hd_metadata_file <- "results/metadata/metadata_HD.tsv"
pd_metadata_file <- "results/metadata/metadata_PD.tsv"


# Diretório de saída
output_dir <- "results/counts"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# VERIFICAR ARQUIVOS DE ENTRADA
input_files <- c(
  ad_file,
  hd_file,
  pd_file,
  ad_metadata_file,
  hd_metadata_file,
  pd_metadata_file
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Arquivos de entrada não encontrados:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


# ==============================================================================
# FUNÇÕES AUXILIARES
# ==============================================================================
# Verificar GeneID
validate_gene_id <- function(dt, dataset_name) {
  if (!"GeneID" %in% names(dt)) {
    stop(dataset_name, ": coluna GeneID não encontrada.")}

  gene_id <- as.character(dt$GeneID)

  missing_gene_id <- (
    is.na(gene_id) |
    trimws(gene_id) == ""
  )

  if (any(missing_gene_id)) {
    stop(dataset_name, ": existem GeneIDs ausentes.")}

  if (anyDuplicated(gene_id)) {
    stop(dataset_name, ": existem GeneIDs duplicados.")}

  invisible(TRUE)
}

# Verificar se As colunas da matriz de expressão correspondem exatamente às linhas do metadata
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

  if (anyNA(metadata_samples)) {
    stop(
      dataset_name,
      ": existem sample_id ausentes no metadata."
    )
  }

  if (anyDuplicated(metadata_samples)) {
    stop(
      dataset_name,
      ": existem sample_id duplicados no metadata."
    )
  }

  if (length(count_samples) != nrow(metadata)) {
    stop(
      dataset_name,
      ": número de amostras da matriz (",
      length(count_samples),
      ") difere do metadata (",
      nrow(metadata),
      ")."
    )
  }

  missing_in_metadata <- setdiff(
    count_samples,
    metadata_samples
  )

  missing_in_counts <- setdiff(
    metadata_samples,
    count_samples
  )

  if (length(missing_in_metadata) > 0) {
    stop(
      dataset_name,
      ": amostras da matriz ausentes no metadata: ",
      paste(
        missing_in_metadata,
        collapse = ", "
      )
    )
  }

  if (length(missing_in_counts) > 0) {
    stop(
      dataset_name,
      ": amostras do metadata ausentes na matriz: ",
      paste(
        missing_in_counts,
        collapse = ", "
      )
    )
  }

  if (!identical(
    count_samples,
    metadata_samples
  )) {

    stop(
      dataset_name,
      ": a ordem das amostras na matriz não corresponde ",
      "à ordem do metadata."
    )
  }

  invisible(TRUE)
}

# Validar matriz de contagens final
validate_final_counts <- function(
  counts,
  dataset_name
) {

  sample_columns <- setdiff(
    names(counts),
    "GeneID"
  )

  x <- as.matrix(
    counts[
      ,
      ..sample_columns
    ]
  )

  storage.mode(x) <- "numeric"

  n_na <- sum(
    is.na(x)
  )

  n_negative <- sum(
    x < 0,
    na.rm = TRUE
  )

  n_fractional <- sum(
    abs(x - round(x)) > 1e-8,
    na.rm = TRUE
  )

  if (n_na > 0) {
    stop(
      dataset_name,
      ": existem ",
      n_na,
      " valores NA na matriz final."
    )
  }

  if (n_negative > 0) {
    stop(
      dataset_name,
      ": existem ",
      n_negative,
      " contagens negativas na matriz final."
    )
  }

  if (n_fractional > 0) {
    stop(
      dataset_name,
      ": existem ",
      n_fractional,
      " valores não inteiros na matriz final."
    )
  }

  invisible(TRUE)
}

# Resumo de uma matriz preparada
summarize_counts <- function(
  counts,
  dataset_name
) {

  sample_columns <- setdiff(
    names(counts),
    "GeneID"
  )

  x <- as.matrix(
    counts[
      ,
      ..sample_columns
    ]
  )

  storage.mode(x) <- "numeric"

  cat("\n")
  cat("--------------------------------------------\n")
  cat(dataset_name, "\n")
  cat("--------------------------------------------\n")

  cat("Genes:", nrow(counts), "\n")
  cat("Amostras:", length(sample_columns), "\n")
  cat("Valores de contagem:", length(x), "\n")
  cat("NA:", sum(is.na(x)), "\n")

  cat("Negativos:",
    sum(
      x < 0,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Não inteiros:",
    sum(
      abs(x - round(x)) > 1e-8,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Zeros:",
    sum(
      x == 0,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Mínimo:",
    min(
      x,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Máximo:",
    max(
      x,
      na.rm = TRUE
    ),
    "\n"
  )
}


# ==============================================================================
# LER METADADOS
cat("\n")
cat("============================================\n")
cat("LEITURA DOS METADADOS\n")
cat("============================================\n")

ad_metadata <- fread(
  ad_metadata_file,
  na.strings = c(
    "NA",
    ""
  )
)

hd_metadata <- fread(
  hd_metadata_file,
  na.strings = c(
    "NA",
    ""
  )
)

pd_metadata <- fread(
  pd_metadata_file,
  na.strings = c(
    "NA",
    ""
  )
)

cat(
  "AD:",
  nrow(ad_metadata),
  "amostras\n"
)

cat(
  "HD:",
  nrow(hd_metadata),
  "amostras\n"
)

cat(
  "PD:",
  nrow(pd_metadata),
  "amostras\n"
)

if (nrow(ad_metadata) != 17) {
  stop("Metadata AD deveria possuir 17 amostras.")
}

if (nrow(hd_metadata) != 69) {
  stop("Metadata HD deveria possuir 69 amostras.")
}

if (nrow(pd_metadata) != 73) {
  stop("Metadata PD deveria possuir 73 amostras.")
}


# ==============================================================================
# PREPARAR AD

cat("\n")
cat("============================================\n")
cat("PREPARAÇÃO AD / GSE53697\n")
cat("============================================\n")

ad_raw <- fread(
  ad_file,
  check.names = FALSE
)

# Validar estrutura original
required_ad_annotation <- c(
  "GeneID",
  "GeneSymbol"
)

if (!all(
  required_ad_annotation %in%
    names(ad_raw)
)) {

  stop(
    "AD: GeneID e/ou GeneSymbol não encontrados."
  )
}

ad_raw_columns <- grep(
  "_raw$",
  names(ad_raw),
  value = TRUE
)

ad_rpkm_columns <- grep(
  "_rpkm$",
  names(ad_raw),
  value = TRUE
)


if (length(ad_raw_columns) != 17) {

  stop(
    "AD: esperadas 17 colunas _raw; encontradas ",
    length(ad_raw_columns),
    "."
  )
}

if (length(ad_rpkm_columns) != 17) {

  stop(
    "AD: esperadas 17 colunas _rpkm; encontradas ",
    length(ad_rpkm_columns),
    "."
  )
}


# Criar matriz contendo SOMENTE GeneID + raw counts
ad_counts <- ad_raw[
  ,
  c(
    "GeneID",
    ad_raw_columns
  ),
  with = FALSE
]

validate_gene_id(
  ad_counts,
  "AD"
)


# Auditoria ANTES do arredondamento
ad_before <- as.matrix(
  ad_counts[
    ,
    ..ad_raw_columns
  ]
)

storage.mode(ad_before) <- "numeric"


if (anyNA(ad_before)) {
  stop(
    "AD: existem valores NA antes do arredondamento."
  )
}


if (any(ad_before < 0)) {
  stop(
    "AD: existem contagens negativas antes do arredondamento."
  )
}


fractional_mask <- (
  abs(
    ad_before -
      round(ad_before)
  ) > 1e-8
)

n_fractional_before <- sum(
  fractional_mask
)

total_values_ad <- length(
  ad_before
)

fractional_percent <- (
  100 *
    n_fractional_before /
    total_values_ad
)

cat(
  "Genes:",
  nrow(ad_counts),
  "\n"
)

cat(
  "Amostras raw:",
  length(ad_raw_columns),
  "\n"
)

cat(
  "Colunas RPKM excluídas da matriz preparada:",
  length(ad_rpkm_columns),
  "\n"
)

cat(
  "Valores totais:",
  total_values_ad,
  "\n"
)

cat(
  "Valores fracionários antes do arredondamento:",
  n_fractional_before,
  "\n"
)

cat(
  sprintf(
    "Percentual fracionário: %.4f%%\n",
    fractional_percent
  )
)


# Arredondamento dos valores não inteiros no objeto ad_counts. para o inteiro mais próximo.
ad_after <- round(
  ad_before
)

# Auditoria do efeito do arredondamento
rounding_difference <- (
  ad_after -
    ad_before
)

absolute_rounding_difference <- abs(
  rounding_difference
)

n_changed <- sum(
  absolute_rounding_difference > 1e-8
)

max_abs_change <- max(
  absolute_rounding_difference
)

total_before <- sum(
  ad_before
)

total_after <- sum(
  ad_after
)

total_delta <- (
  total_after -
    total_before
)

n_fractional_after <- sum(
  abs(
    ad_after -
      round(ad_after)
  ) > 1e-8
)


cat("\nAUDITORIA DO ARREDONDAMENTO AD\n")
cat("--------------------------------------------\n")

cat(
  "Valores alterados:",
  n_changed,
  "\n"
)

cat(
  "Valores não inteiros antes:",
  n_fractional_before,
  "\n"
)

cat(
  "Valores não inteiros depois:",
  n_fractional_after,
  "\n"
)

cat(
  "Máxima alteração absoluta:",
  format(
    max_abs_change,
    digits = 12
  ),
  "\n"
)

cat(
  "Soma das contagens antes:",
  format(
    total_before,
    scientific = FALSE,
    digits = 15
  ),
  "\n"
)

cat(
  "Soma das contagens depois:",
  format(
    total_after,
    scientific = FALSE,
    digits = 15
  ),
  "\n"
)

cat(
  "Diferença total após arredondamento:",
  format(
    total_delta,
    scientific = FALSE,
    digits = 15
  ),
  "\n"
)

# Garantias específicas do arredondamento
if (n_fractional_before == 0) {
  stop(
    "AD: nenhum valor fracionário foi encontrado. ",
    "Isso difere da auditoria realizada na etapa 00."
  )
}


if (max_abs_change > 0.5 + 1e-8) {
  stop(
    "AD: alteração de arredondamento superior a 0.5 detectada."
  )
}


if (n_fractional_after != 0) {
  stop(
    "AD: ainda existem valores não inteiros após o arredondamento."
  )
}


# Substituir somente as colunas de contagem no objeto preparado
ad_counts[
  ,
  (ad_raw_columns) := as.data.table(
    ad_after
  )
]

# Liberar objeto bruto grande
rm(
  ad_raw,
  ad_before,
  ad_after,
  rounding_difference,
  absolute_rounding_difference,
  fractional_mask
)

gc()


# ==============================================================================
# PREPARAR HD
cat("\n")
cat("============================================\n")
cat("PREPARAÇÃO HD / GSE64810\n")
cat("============================================\n")

hd_counts <- fread(
  hd_file,
  check.names = FALSE
)

validate_gene_id(
  hd_counts,
  "HD"
)

hd_sample_columns <- setdiff(
  names(hd_counts),
  "GeneID"
)


if (length(hd_sample_columns) != 69) {
  stop(
    "HD: esperadas 69 amostras; encontradas ",
    length(hd_sample_columns),
    "."
  )
}


hd_matrix <- as.matrix(
  hd_counts[
    ,
    ..hd_sample_columns
  ]
)

storage.mode(hd_matrix) <- "numeric"


if (anyNA(hd_matrix)) {
  stop(
    "HD: existem valores NA."
  )
}


if (any(hd_matrix < 0)) {
  stop(
    "HD: existem contagens negativas."
  )
}


hd_fractional <- sum(
  abs(
    hd_matrix -
      round(hd_matrix)
  ) > 1e-8
)


if (hd_fractional != 0) {

  stop(
    "HD: foram encontrados ",
    hd_fractional,
    " valores não inteiros."
  )
}


cat(
  "Genes:",
  nrow(hd_counts),
  "\n"
)

cat(
  "Amostras:",
  length(hd_sample_columns),
  "\n"
)

cat(
  "Valores não inteiros:",
  hd_fractional,
  "\n"
)

cat(
  "Transformação numérica aplicada: nenhuma\n"
)


rm(
  hd_matrix
)

gc()


# PREPARAR PD
cat("\n")
cat("============================================\n")
cat("PREPARAÇÃO PD / GSE68719\n")
cat("============================================\n")

pd_counts <- fread(
  pd_file,
  check.names = FALSE
)

validate_gene_id(
  pd_counts,
  "PD"
)

pd_sample_columns <- setdiff(
  names(pd_counts),
  "GeneID"
)


if (length(pd_sample_columns) != 73) {
  stop(
    "PD: esperadas 73 amostras; encontradas ",
    length(pd_sample_columns),
    "."
  )
}

pd_matrix <- as.matrix(
  pd_counts[
    ,
    ..pd_sample_columns
  ]
)

storage.mode(pd_matrix) <- "numeric"

if (anyNA(pd_matrix)) {
  stop(
    "PD: existem valores NA."
  )
}

if (any(pd_matrix < 0)) {
  stop(
    "PD: existem contagens negativas."
  )
}

pd_fractional <- sum(
  abs(
    pd_matrix -
      round(pd_matrix)
  ) > 1e-8
)


if (pd_fractional != 0) {
  stop(
    "PD: foram encontrados ",
    pd_fractional,
    " valores não inteiros."
  )
}

cat(
  "Genes:",
  nrow(pd_counts),
  "\n"
)

cat(
  "Amostras:",
  length(pd_sample_columns),
  "\n"
)

cat(
  "Valores não inteiros:",
  pd_fractional,
  "\n"
)

cat(
  "Transformação numérica aplicada: nenhuma\n"
)

rm(
  pd_matrix
)

gc()


# ==============================================================================
# VALIDAR CORRESPONDÊNCIA COUNTS <-> METADATA
cat("\n")
cat("============================================\n")
cat("VALIDAÇÃO COUNTS <-> METADATA\n")
cat("============================================\n")

validate_sample_mapping(
  ad_counts,
  ad_metadata,
  "AD"
)

cat(
  "AD: 17 / 17 amostras correspondentes e na ordem correta\n"
)

validate_sample_mapping(
  hd_counts,
  hd_metadata,
  "HD"
)

cat(
  "HD: 69 / 69 amostras correspondentes e na ordem correta\n"
)

validate_sample_mapping(
  pd_counts,
  pd_metadata,
  "PD"
)

cat(
  "PD: 73 / 73 amostras correspondentes e na ordem correta\n"
)

# ==============================================================================
# VALIDAÇÃO FINAL DAS MATRIZES
cat("\n")
cat("============================================\n")
cat("VALIDAÇÃO FINAL DAS CONTAGENS\n")
cat("============================================\n")

validate_final_counts(
  ad_counts,
  "AD"
)

validate_final_counts(
  hd_counts,
  "HD"
)

validate_final_counts(
  pd_counts,
  "PD"
)

cat(
  "AD: matriz válida\n"
)

cat(
  "HD: matriz válida\n"
)

cat(
  "PD: matriz válida\n"
)

# RESUMO DAS MATRIZES PREPARADAS
cat("\n")
cat("============================================\n")
cat("RESUMO DAS MATRIZES PREPARADAS\n")
cat("============================================\n")

summarize_counts(
  ad_counts,
  "AD / GSE53697"
)

summarize_counts(
  hd_counts,
  "HD / GSE64810"
)

summarize_counts(
  pd_counts,
  "PD / GSE68719"
)


ad_output <- file.path(
  output_dir,
  "counts_AD.tsv"
)

hd_output <- file.path(
  output_dir,
  "counts_HD.tsv"
)

pd_output <- file.path(
  output_dir,
  "counts_PD.tsv"
)


fwrite(
  ad_counts,
  ad_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  hd_counts,
  hd_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  pd_counts,
  pd_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# VERIFICAR ARQUIVOS GRAVADOS
if (!file.exists(ad_output)) {
  stop(
    "Falha ao gerar ",
    ad_output
  )
}

if (!file.exists(hd_output)) {
  stop(
    "Falha ao gerar ",
    hd_output
  )
}

if (!file.exists(pd_output)) {
  stop(
    "Falha ao gerar ",
    pd_output
  )
}


cat("\n")
cat("============================================\n")
cat("02_prepare_counts.R CONCLUÍDO\n")
cat("============================================\n")

cat(
  "AD:",
  nrow(ad_counts),
  "genes x",
  ncol(ad_counts) - 1,
  "amostras\n"
)

cat(
  "HD:",
  nrow(hd_counts),
  "genes x",
  ncol(hd_counts) - 1,
  "amostras\n"
)

cat(
  "PD:",
  nrow(pd_counts),
  "genes x",
  ncol(pd_counts) - 1,
  "amostras\n"
)

cat("\nArredondamento AD:\n")

cat(
  " - valores fracionários antes:",
  n_fractional_before,
  "\n"
)

cat(
  " - valores fracionários depois:",
  n_fractional_after,
  "\n"
)

cat(
  " - máxima alteração absoluta:",
  format(
    max_abs_change,
    digits = 12
  ),
  "\n"
)

cat(
  " - diferença na soma total:",
  format(
    total_delta,
    scientific = FALSE,
    digits = 15
  ),
  "\n"
)


cat("\nArquivos gerados:\n")

cat(
  " -",
  ad_output,
  "\n"
)

cat(
  " -",
  hd_output,
  "\n"
)

cat(
  " -",
  pd_output,
  "\n"
)


cat("\nNenhum gene foi removido.\n")
cat("Nenhum filtro de expressão foi aplicado.\n")
cat("Nenhuma normalização foi realizada.\n")
cat("Nenhuma harmonização entre datasets foi realizada.\n")
cat("HD e PD foram preservados numericamente.\n")
cat("Os arquivos brutos originais permanecem inalterados.\n")