# ==============================================================================
# 02_prepare_counts_ROSMAP.R
# ==============================================================================
#
# Fase 01c — ROSMAP
# Preparação estrutural da matriz de counts para a análise principal
#
# Este script deve ser executado APÓS:
#
#   00_audit_ROSMAP.R
#   01_prepare_metadata_ROSMAP.R
#
# Objetivos
# ---------
# 1. Ler o metadata principal já validado:
#
#       results/candidate_ROSMAP/metadata/metadata_ROSMAP.tsv
#
# 2. Selecionar da matriz original SOMENTE as 144 amostras da análise principal:
#
#       Control = cogdx 1  → 67
#       AD      = cogdx 4  → 77
#
# 3. Preservar EXATAMENTE a ordem de sample_id definida no metadata.
#
# 4. Remover apenas as quatro linhas técnicas do featureCounts:
#
#       N_unmapped
#       N_multimapping
#       N_noFeature
#       N_ambiguous
#
# 5. Manter os identificadores gênicos ORIGINAIS e VERSIONADOS do ROSMAP.
#
#    Exemplo:
#
#       ENSG00000186092.6
#
#    A remoção da versão e o mapeamento para NCBI GeneID ocorrerão somente em:
#
#       03_harmonize_genes_ROSMAP.R
#
# 6. Auditar novamente a integridade da matriz selecionada:
#
#       - sample IDs;
#       - ordem counts ↔ metadata;
#       - feature IDs;
#       - duplicatas;
#       - NA;
#       - valores não finitos;
#       - valores negativos;
#       - valores fracionários;
#       - genes all-zero;
#       - amostras all-zero;
#       - library size;
#       - genes detectados por amostra.
#
# 7. Salvar uma matriz estruturalmente limpa, mas AINDA NÃO FILTRADA por
#    baixa expressão.
#
# IMPORTANTE
# ----------
# Este script NÃO:
#   - filtra genes por count >= 10;
#   - normaliza counts;
#   - executa DESeq2;
#   - transforma counts;
#   - arredonda counts;
#   - remove versão dos IDs Ensembl;
#   - converte Ensembl para GeneID;
#   - colapsa genes;
#   - remove genes all-zero.
#
# O filtro de baixa expressão será realizado posteriormente no 04_qc_counts,
# após a harmonização dos identificadores gênicos.
#
# ==============================================================================


# ==============================================================================
# 0. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table"
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
})


# ==============================================================================
# 1. CAMINHOS
# ==============================================================================

raw_dir <- file.path(
  "data",
  "raw",
  "ROSMAP"
)

counts_input <- file.path(
  raw_dir,
  "ROSMAP_batch2_gene_all_counts_matrix_clean.txt"
)

metadata_input <- file.path(
  "results",
  "candidate_ROSMAP",
  "metadata",
  "metadata_ROSMAP.tsv"
)

candidate_dir <- file.path(
  "results",
  "candidate_ROSMAP"
)

counts_dir <- file.path(
  candidate_dir,
  "counts"
)

dir.create(
  counts_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Arquivos de saída
# ------------------------------------------------------------------------------

counts_output <- file.path(
  counts_dir,
  "counts_ROSMAP_raw_main.tsv.gz"
)

sample_qc_output <- file.path(
  counts_dir,
  "counts_ROSMAP_sample_qc.tsv"
)

feature_audit_output <- file.path(
  counts_dir,
  "counts_ROSMAP_feature_audit.tsv"
)

technical_features_output <- file.path(
  counts_dir,
  "counts_ROSMAP_technical_features_removed.tsv"
)

nonstandard_feature_output <- file.path(
  counts_dir,
  "counts_ROSMAP_nonstandard_feature_ids.tsv"
)

dimensions_output <- file.path(
  counts_dir,
  "counts_ROSMAP_dimensions.tsv"
)

report_output <- file.path(
  counts_dir,
  "02_prepare_counts_ROSMAP.txt"
)


# ==============================================================================
# 2. CONSTANTES E EXPECTATIVAS PRÉ-DEFINIDAS
# ==============================================================================

EXPECTED_MAIN_N <- 144L
EXPECTED_CONTROL_N <- 67L
EXPECTED_AD_N <- 77L

TECHNICAL_FEATURES_EXPECTED <- c(
  "N_unmapped",
  "N_multimapping",
  "N_noFeature",
  "N_ambiguous"
)

ENSEMBL_VERSIONED_PATTERN <- "^ENSG[0-9]+\\.[0-9]+$"
ENSEMBL_UNVERSIONED_PATTERN <- "^ENSG[0-9]+$"


# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================

assert_file <- function(path) {

  if (!file.exists(path)) {

    stop(
      "Arquivo não encontrado:\n",
      path
    )
  }

  invisible(TRUE)
}


assert_columns <- function(
  dt,
  required,
  object_name
) {

  missing <- setdiff(
    required,
    names(dt)
  )

  if (length(missing) > 0) {

    stop(
      object_name,
      ": colunas obrigatórias ausentes: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}


fmt_int <- function(x) {

  format(
    as.integer(
      round(
        as.numeric(x)
      )
    ),
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


fmt_num <- function(
  x,
  digits = 2
) {

  format(
    round(
      as.numeric(x),
      digits
    ),
    nsmall = digits,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


count_fractional_values <- function(x) {

  sum(
    abs(
      x -
        round(
          x
        )
    ) > 1e-8,
    na.rm = TRUE
  )
}


# ==============================================================================
# 4. VALIDAR ARQUIVOS DE ENTRADA
# ==============================================================================

assert_file(
  counts_input
)

assert_file(
  metadata_input
)


# ==============================================================================
# 5. LER E VALIDAR METADATA PRINCIPAL
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "02 — PREPARAÇÃO DOS COUNTS ROSMAP\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Lendo metadata principal...\n"
)

metadata <- fread(
  metadata_input,
  na.strings = c(
    "NA",
    ""
  )
)

assert_columns(
  metadata,
  c(
    "sample_id",
    "individual_id",
    "group",
    "cogdx"
  ),
  "metadata_ROSMAP.tsv"
)

metadata[
  ,
  sample_id := trimws(
    as.character(
      sample_id
    )
  )
]

metadata[
  ,
  individual_id := trimws(
    as.character(
      individual_id
    )
  )
]

metadata[
  ,
  group := trimws(
    as.character(
      group
    )
  )
]


# ------------------------------------------------------------------------------
# 5.1 Integridade do metadata
# ------------------------------------------------------------------------------

if (
  anyNA(
    metadata$sample_id
  ) ||
    any(
      metadata$sample_id == ""
    )
) {

  stop(
    "metadata principal contém sample_id ausente/vazio."
  )
}

if (
  anyDuplicated(
    metadata$sample_id
  )
) {

  stop(
    "metadata principal contém sample_id duplicado."
  )
}

if (
  anyNA(
    metadata$individual_id
  ) ||
    any(
      metadata$individual_id == ""
    )
) {

  stop(
    "metadata principal contém individual_id ausente/vazio."
  )
}

if (
  anyDuplicated(
    metadata$individual_id
  )
) {

  stop(
    "metadata principal contém individual_id duplicado."
  )
}

if (
  !all(
    metadata$group %in% c(
      "Control",
      "AD"
    )
  )
) {

  stop(
    "metadata principal contém valores inesperados em group."
  )
}


# ------------------------------------------------------------------------------
# 5.2 Tamanho esperado
# ------------------------------------------------------------------------------

if (
  nrow(
    metadata
  ) != EXPECTED_MAIN_N
) {

  stop(
    "Número inesperado de amostras no metadata principal. Esperado: ",
    EXPECTED_MAIN_N,
    "; observado: ",
    nrow(
      metadata
    ),
    "."
  )
}

group_counts <- metadata[
  ,
  .N,
  by = group
]

n_control <- group_counts[
  group == "Control",
  N
]

n_ad <- group_counts[
  group == "AD",
  N
]

if (
  length(
    n_control
  ) != 1 ||
    n_control != EXPECTED_CONTROL_N
) {

  stop(
    "Número inesperado de controles. Esperado: ",
    EXPECTED_CONTROL_N,
    "."
  )
}

if (
  length(
    n_ad
  ) != 1 ||
    n_ad != EXPECTED_AD_N
) {

  stop(
    "Número inesperado de casos AD. Esperado: ",
    EXPECTED_AD_N,
    "."
  )
}

cat(
  "Metadata principal validado: ",
  fmt_int(
    nrow(
      metadata
    )
  ),
  " amostras (",
  fmt_int(
    n_ad
  ),
  " AD; ",
  fmt_int(
    n_control
  ),
  " Control).\n\n",
  sep = ""
)


# ==============================================================================
# 6. LER CABEÇALHO DA MATRIZ ORIGINAL
# ==============================================================================

cat(
  "Lendo cabeçalho da matriz original...\n"
)

counts_header <- fread(
  counts_input,
  nrows = 0
)

header_names <- names(
  counts_header
)

if (
  length(
    header_names
  ) < 2
) {

  stop(
    "A matriz de counts possui menos de duas colunas."
  )
}

feature_col <- header_names[1]

counts_sample_ids_all <- header_names[
  -1
]


# ------------------------------------------------------------------------------
# 6.1 Validar que todas as 144 amostras existem na matriz
# ------------------------------------------------------------------------------

missing_in_counts <- setdiff(
  metadata$sample_id,
  counts_sample_ids_all
)

if (
  length(
    missing_in_counts
  ) > 0
) {

  stop(
    "Há sample_id do metadata ausentes na matriz de counts: ",
    paste(
      missing_in_counts,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# 6.2 Validar que as 144 amostras são DLPFC pela nomenclatura
# ------------------------------------------------------------------------------

non_dlpfc_main <- metadata$sample_id[
  !grepl(
    "-DLPFC$",
    metadata$sample_id
  )
]

if (
  length(
    non_dlpfc_main
  ) > 0
) {

  stop(
    "Existem amostras principais sem sufixo -DLPFC: ",
    paste(
      non_dlpfc_main,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# 6.3 Confirmar que o metadata já está na ordem relativa dos counts
# ------------------------------------------------------------------------------

expected_order_from_counts <- counts_sample_ids_all[
  counts_sample_ids_all %in% metadata$sample_id
]

if (
  !identical(
    as.character(
      metadata$sample_id
    ),
    as.character(
      expected_order_from_counts
    )
  )
) {

  stop(
    paste0(
      "A ordem das amostras no metadata principal não corresponde à ordem ",
      "relativa das mesmas amostras na matriz original."
    )
  )
}

cat(
  "As 144 amostras existem na matriz e a ordem counts ↔ metadata é compatível.\n\n"
)


# ==============================================================================
# 7. LER SOMENTE feature + 144 AMOSTRAS PRINCIPAIS
# ==============================================================================

cat(
  "Lendo apenas feature + 144 amostras principais...\n"
)

columns_to_read <- c(
  feature_col,
  metadata$sample_id
)

counts <- fread(
  counts_input,
  select = columns_to_read,
  na.strings = c(
    "NA",
    ""
  )
)

if (
  ncol(
    counts
  ) != (
    nrow(
      metadata
    ) + 1L
  )
) {

  stop(
    "Número inesperado de colunas após seleção da matriz principal."
  )
}

setnames(
  counts,
  1,
  "feature"
)

counts[
  ,
  feature := trimws(
    as.character(
      feature
    )
  )
]


# ------------------------------------------------------------------------------
# 7.1 Validar ordem das colunas
# ------------------------------------------------------------------------------

if (
  !identical(
    names(
      counts
    )[-1],
    metadata$sample_id
  )
) {

  stop(
    "A ordem das colunas da matriz selecionada não corresponde ao metadata."
  )
}


# ==============================================================================
# 8. AUDITAR FEATURES ANTES DE REMOVER AS LINHAS TÉCNICAS
# ==============================================================================

n_rows_input <- nrow(
  counts
)

n_missing_feature_ids <- sum(
  is.na(
    counts$feature
  ) |
    counts$feature == ""
)

if (
  n_missing_feature_ids > 0
) {

  stop(
    "Foram encontrados feature IDs ausentes/vazios."
  )
}

n_feature_duplicates_before <- sum(
  duplicated(
    counts$feature
  )
)

if (
  n_feature_duplicates_before > 0
) {

  duplicated_ids <- unique(
    counts$feature[
      duplicated(
        counts$feature
      ) |
        duplicated(
          counts$feature,
          fromLast = TRUE
        )
    ]
  )

  stop(
    "Foram encontrados feature IDs duplicados na matriz original. Exemplos: ",
    paste(
      head(
        duplicated_ids,
        20
      ),
      collapse = ", "
    )
  )
}


# ==============================================================================
# 9. IDENTIFICAR E REMOVER APENAS FEATURES TÉCNICAS N_*
# ==============================================================================

technical_mask <- grepl(
  "^N_",
  counts$feature
)

technical_features <- counts[
  technical_mask,
  .(
    feature
  )
]

technical_ids_observed <- technical_features$feature


# ------------------------------------------------------------------------------
# 9.1 Registrar qualquer diferença em relação às quatro linhas esperadas
# ------------------------------------------------------------------------------

technical_missing_expected <- setdiff(
  TECHNICAL_FEATURES_EXPECTED,
  technical_ids_observed
)

technical_unexpected <- setdiff(
  technical_ids_observed,
  TECHNICAL_FEATURES_EXPECTED
)

if (
  length(
    technical_missing_expected
  ) > 0
) {

  warning(
    "Features técnicas esperadas não encontradas: ",
    paste(
      technical_missing_expected,
      collapse = ", "
    )
  )
}

if (
  length(
    technical_unexpected
  ) > 0
) {

  warning(
    "Features N_* adicionais encontradas: ",
    paste(
      technical_unexpected,
      collapse = ", "
    )
  )
}

fwrite(
  technical_features,
  technical_features_output,
  sep = "\t"
)


# ------------------------------------------------------------------------------
# 9.2 Remover todas as linhas N_*
# ------------------------------------------------------------------------------

counts_gene <- counts[
  !technical_mask
]

n_technical_removed <- sum(
  technical_mask
)

n_gene_rows_after_technical <- nrow(
  counts_gene
)

cat(
  "Linhas totais antes da remoção técnica: ",
  fmt_int(
    n_rows_input
  ),
  "\n",
  sep = ""
)

cat(
  "Linhas técnicas N_* removidas: ",
  fmt_int(
    n_technical_removed
  ),
  "\n",
  sep = ""
)

cat(
  "Features gênicas restantes: ",
  fmt_int(
    n_gene_rows_after_technical
  ),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 10. AUDITAR IDENTIFICADORES GÊNICOS
# ==============================================================================

if (
  anyDuplicated(
    counts_gene$feature
  )
) {

  stop(
    "Existem feature IDs duplicados após remover as linhas técnicas."
  )
}

is_versioned_ensembl <- grepl(
  ENSEMBL_VERSIONED_PATTERN,
  counts_gene$feature
)

is_unversioned_ensembl <- grepl(
  ENSEMBL_UNVERSIONED_PATTERN,
  counts_gene$feature
)

is_nonstandard <- !(
  is_versioned_ensembl |
    is_unversioned_ensembl
)

n_versioned_ensembl <- sum(
  is_versioned_ensembl
)

n_unversioned_ensembl <- sum(
  is_unversioned_ensembl
)

n_nonstandard <- sum(
  is_nonstandard
)


# ------------------------------------------------------------------------------
# 10.1 NÃO excluir os IDs não padrão.
#
# Eles serão apenas registrados para investigação no 03_harmonize_genes.
# ------------------------------------------------------------------------------

nonstandard_features <- counts_gene[
  is_nonstandard,
  .(
    feature
  )
]

fwrite(
  nonstandard_features,
  nonstandard_feature_output,
  sep = "\t"
)


# ==============================================================================
# 11. CONVERTER COUNTS PARA MATRIZ NUMÉRICA PARA VALIDAÇÃO
# ==============================================================================

cat(
  "Validando valores da matriz...\n"
)

count_matrix <- as.matrix(
  counts_gene[
    ,
    -1,
    with = FALSE
  ]
)

storage.mode(
  count_matrix
) <- "numeric"


# ==============================================================================
# 12. VALIDAR INTEGRIDADE DOS COUNTS
# ==============================================================================

n_na <- sum(
  is.na(
    count_matrix
  )
)

n_nonfinite <- sum(
  !is.finite(
    count_matrix
  ),
  na.rm = TRUE
)

n_negative <- sum(
  count_matrix < 0,
  na.rm = TRUE
)

n_fractional <- count_fractional_values(
  count_matrix
)

if (
  n_na > 0
) {

  stop(
    "Foram encontrados valores NA na matriz principal."
  )
}

if (
  n_nonfinite > 0
) {

  stop(
    "Foram encontrados valores não finitos na matriz principal."
  )
}

if (
  n_negative > 0
) {

  stop(
    "Foram encontrados counts negativos na matriz principal."
  )
}

if (
  n_fractional > 0
) {

  stop(
    paste0(
      "Foram encontrados counts fracionários. ",
      "Este script NÃO arredonda valores automaticamente."
    )
  )
}


# ==============================================================================
# 13. AUDITAR GENES E AMOSTRAS ALL-ZERO
# ==============================================================================

gene_totals <- rowSums(
  count_matrix
)

sample_totals <- colSums(
  count_matrix
)

all_zero_gene_mask <- gene_totals == 0
all_zero_sample_mask <- sample_totals == 0

n_all_zero_genes <- sum(
  all_zero_gene_mask
)

n_all_zero_samples <- sum(
  all_zero_sample_mask
)

if (
  n_all_zero_samples > 0
) {

  zero_samples <- colnames(
    count_matrix
  )[
    all_zero_sample_mask
  ]

  stop(
    "Foram encontradas amostras totalmente zeradas: ",
    paste(
      zero_samples,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# Genes all-zero são AUDITADOS, mas NÃO removidos neste script.
# A filtragem será feita no 04_qc_counts.
# ------------------------------------------------------------------------------

cat(
  "Genes all-zero na análise principal: ",
  fmt_int(
    n_all_zero_genes
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras all-zero: ",
  fmt_int(
    n_all_zero_samples
  ),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 14. CALCULAR MÉTRICAS POR AMOSTRA
# ==============================================================================

library_size <- colSums(
  count_matrix
)

detected_genes_gt0 <- colSums(
  count_matrix > 0
)

detected_genes_ge10 <- colSums(
  count_matrix >= 10
)


sample_qc <- data.table(
  sample_id = colnames(
    count_matrix
  ),
  group = metadata$group,
  cogdx = metadata$cogdx,
  individual_id = metadata$individual_id,
  library_size = as.numeric(
    library_size
  ),
  detected_genes_gt0 = as.integer(
    detected_genes_gt0
  ),
  detected_genes_ge10 = as.integer(
    detected_genes_ge10
  )
)


# ------------------------------------------------------------------------------
# Acrescentar covariáveis úteis se presentes.
# ------------------------------------------------------------------------------

optional_metadata_cols <- c(
  "sequencingBatch",
  "msex",
  "RIN",
  "pmi",
  "age_death_raw",
  "braaksc",
  "ceradsc"
)

for (col in optional_metadata_cols) {

  if (col %in% names(metadata)) {

    sample_qc[
      ,
      (col) := metadata[[col]]
    ]
  }
}

fwrite(
  sample_qc,
  sample_qc_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 15. AUDITORIA DE FEATURES
# ==============================================================================

feature_audit <- data.table(
  metric = c(
    "rows_input_including_technical",
    "technical_features_removed",
    "gene_features_after_technical_removal",
    "missing_feature_ids",
    "duplicated_feature_ids_before",
    "duplicated_feature_ids_after",
    "versioned_Ensembl_ids",
    "unversioned_Ensembl_ids",
    "nonstandard_gene_feature_ids",
    "all_zero_genes_main_144",
    "NA_values",
    "nonfinite_values",
    "negative_values",
    "fractional_values",
    "all_zero_samples"
  ),
  value = c(
    n_rows_input,
    n_technical_removed,
    n_gene_rows_after_technical,
    n_missing_feature_ids,
    n_feature_duplicates_before,
    sum(
      duplicated(
        counts_gene$feature
      )
    ),
    n_versioned_ensembl,
    n_unversioned_ensembl,
    n_nonstandard,
    n_all_zero_genes,
    n_na,
    n_nonfinite,
    n_negative,
    n_fractional,
    n_all_zero_samples
  )
)

fwrite(
  feature_audit,
  feature_audit_output,
  sep = "\t"
)


# ==============================================================================
# 16. DIMENSÕES DA MATRIZ PREPARADA
# ==============================================================================

dimensions <- data.table(
  object = c(
    "original_header_samples",
    "metadata_main_samples",
    "prepared_matrix_samples",
    "prepared_gene_features"
  ),
  n = c(
    length(
      counts_sample_ids_all
    ),
    nrow(
      metadata
    ),
    ncol(
      count_matrix
    ),
    nrow(
      count_matrix
    )
  )
)

fwrite(
  dimensions,
  dimensions_output,
  sep = "\t"
)


# ==============================================================================
# 17. VALIDAÇÃO FINAL COUNTS ↔ METADATA
# ==============================================================================

if (
  !identical(
    colnames(
      count_matrix
    ),
    metadata$sample_id
  )
) {

  stop(
    "Validação final falhou: colnames(count_matrix) != metadata$sample_id."
  )
}

if (
  ncol(
    count_matrix
  ) != nrow(
    metadata
  )
) {

  stop(
    "Validação final falhou: número de amostras counts != metadata."
  )
}

if (
  ncol(
    count_matrix
  ) != EXPECTED_MAIN_N
) {

  stop(
    "Validação final falhou: a matriz preparada não possui 144 amostras."
  )
}


# ==============================================================================
# 18. SALVAR MATRIZ PREPARADA
# ==============================================================================
#
# A matriz é salva comprimida para reduzir espaço em disco.
#
# Formato:
#
#   feature    sample1    sample2    ...    sample144
#
# Os feature IDs permanecem EXATAMENTE como no ROSMAP original.
# ==============================================================================

cat(
  "Salvando matriz principal preparada...\n"
)

# ------------------------------------------------------------------------------
# data.table::fwrite(compress = "gzip") depende de o pacote data.table ter sido
# compilado com suporte à biblioteca zlib. No ambiente atual, esse suporte pode
# não estar disponível.
#
# Para tornar o pipeline independente dessa configuração, salvamos primeiro um
# TSV normal com fwrite() e depois comprimimos com o executável gzip do sistema.
#
# Isso NÃO altera os dados; muda apenas o método de compressão do arquivo.
# ------------------------------------------------------------------------------

counts_output_uncompressed <- sub(
  "\\.gz$",
  "",
  counts_output
)

fwrite(
  counts_gene,
  counts_output_uncompressed,
  sep = "\t",
  quote = FALSE
)

gzip_bin <- Sys.which(
  "gzip"
)

if (
  identical(
    gzip_bin,
    ""
  )
) {

  stop(
    paste0(
      "A matriz foi salva sem compressão em:\n  ",
      counts_output_uncompressed,
      "\n\n",
      "O executável 'gzip' não foi encontrado no sistema. ",
      "Comprima o arquivo manualmente antes de continuar."
    )
  )
}

gzip_status <- system2(
  command = gzip_bin,
  args = c(
    "-f",
    counts_output_uncompressed
  )
)

if (
  gzip_status != 0
) {

  stop(
    "Falha ao comprimir a matriz com gzip. Código de saída: ",
    gzip_status
  )
}

if (
  !file.exists(
    counts_output
  )
) {

  stop(
    "A compressão terminou sem erro, mas o arquivo esperado não foi encontrado: ",
    counts_output
  )
}

cat(
  "Matriz comprimida com gzip do sistema: ",
  counts_output,
  "\n",
  sep = ""
)


# ==============================================================================
# 19. RELATÓRIO TEXTUAL
# ==============================================================================

sink(
  report_output
)

cat(
  "============================================================\n"
)
cat(
  "FASE 01c — 02_prepare_counts_ROSMAP.R\n"
)
cat(
  "PREPARAÇÃO ESTRUTURAL DOS COUNTS\n"
)
cat(
  "============================================================\n\n"
)


# ------------------------------------------------------------------------------
# 1. Entradas
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "1. ENTRADAS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Counts original:\n  ",
  counts_input,
  "\n\n",
  sep = ""
)

cat(
  "Metadata principal:\n  ",
  metadata_input,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 2. Amostras
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "2. AMOSTRAS SELECIONADAS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Amostras totais no arquivo original: ",
  fmt_int(
    length(
      counts_sample_ids_all
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras da análise principal: ",
  fmt_int(
    nrow(
      metadata
    )
  ),
  "\n",
  sep = ""
)

cat(
  "AD: ",
  fmt_int(
    n_ad
  ),
  "\n",
  sep = ""
)

cat(
  "Control: ",
  fmt_int(
    n_control
  ),
  "\n",
  sep = ""
)

cat(
  "Todas possuem sufixo -DLPFC: TRUE\n"
)

cat(
  "Ordem counts ↔ metadata preservada: TRUE\n\n"
)


# ------------------------------------------------------------------------------
# 3. Features técnicas
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "3. FEATURES TÉCNICAS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Linhas totais antes da remoção técnica: ",
  fmt_int(
    n_rows_input
  ),
  "\n",
  sep = ""
)

cat(
  "Linhas N_* removidas: ",
  fmt_int(
    n_technical_removed
  ),
  "\n",
  sep = ""
)

if (
  length(
    technical_ids_observed
  ) > 0
) {

  cat(
    "IDs removidos:\n"
  )

  for (id in technical_ids_observed) {

    cat(
      "  ",
      id,
      "\n",
      sep = ""
    )
  }
}

cat(
  "\nFeatures gênicas restantes: ",
  fmt_int(
    n_gene_rows_after_technical
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 4. Identificadores
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "4. IDENTIFICADORES GÊNICOS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Ensembl versionados: ",
  fmt_int(
    n_versioned_ensembl
  ),
  "\n",
  sep = ""
)

cat(
  "Ensembl sem versão: ",
  fmt_int(
    n_unversioned_ensembl
  ),
  "\n",
  sep = ""
)

cat(
  "IDs não classificados pelo regex Ensembl: ",
  fmt_int(
    n_nonstandard
  ),
  "\n",
  sep = ""
)

cat(
  "Feature IDs duplicados: 0\n"
)

cat(
  "Feature IDs ausentes: 0\n\n"
)

cat(
  "Nenhum ID gênico foi convertido, truncado ou colapsado nesta etapa.\n"
)
cat(
  "A harmonização será realizada apenas no 03_harmonize_genes_ROSMAP.R.\n\n"
)


# ------------------------------------------------------------------------------
# 5. Integridade numérica
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "5. INTEGRIDADE NUMÉRICA\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "NA: ",
  fmt_int(
    n_na
  ),
  "\n",
  sep = ""
)

cat(
  "Não finitos: ",
  fmt_int(
    n_nonfinite
  ),
  "\n",
  sep = ""
)

cat(
  "Negativos: ",
  fmt_int(
    n_negative
  ),
  "\n",
  sep = ""
)

cat(
  "Fracionários: ",
  fmt_int(
    n_fractional
  ),
  "\n",
  sep = ""
)

cat(
  "Genes all-zero: ",
  fmt_int(
    n_all_zero_genes
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras all-zero: ",
  fmt_int(
    n_all_zero_samples
  ),
  "\n\n",
  sep = ""
)

cat(
  "OBSERVAÇÃO: genes all-zero NÃO foram removidos nesta etapa.\n"
)
cat(
  "A filtragem de baixa expressão será executada no 04_qc_counts.\n\n"
)


# ------------------------------------------------------------------------------
# 6. Library size
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "6. LIBRARY SIZE — 144 AMOSTRAS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Mínimo: ",
  fmt_int(
    min(
      library_size
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Mediana: ",
  fmt_int(
    median(
      library_size
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Máximo: ",
  fmt_int(
    max(
      library_size
    )
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 7. Genes detectados
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "7. GENES DETECTADOS POR AMOSTRA\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Critério > 0 count:\n"
)

cat(
  "  mínimo = ",
  fmt_int(
    min(
      detected_genes_gt0
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  mediana = ",
  fmt_int(
    median(
      detected_genes_gt0
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  máximo = ",
  fmt_int(
    max(
      detected_genes_gt0
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Critério >= 10 counts:\n"
)

cat(
  "  mínimo = ",
  fmt_int(
    min(
      detected_genes_ge10
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  mediana = ",
  fmt_int(
    median(
      detected_genes_ge10
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  máximo = ",
  fmt_int(
    max(
      detected_genes_ge10
    )
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 8. Critério futuro de baixa expressão
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "8. FILTRAGEM DE BAIXA EXPRESSÃO — AINDA NÃO APLICADA\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Menor grupo da análise principal: Control, n = ",
  fmt_int(
    n_control
  ),
  ".\n",
  sep = ""
)

cat(
  "O critério planejado para 04_qc_counts é:\n"
)
cat(
  "  count >= 10 em pelo menos ",
  fmt_int(
    n_control
  ),
  " amostras.\n\n",
  sep = ""
)

cat(
  "Esse filtro NÃO foi aplicado pelo 02_prepare_counts.\n\n"
)


# ------------------------------------------------------------------------------
# 9. Arquivos gerados
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "9. ARQUIVOS GERADOS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Matriz principal preparada:\n  ",
  counts_output,
  "\n\n",
  sep = ""
)

cat(
  "QC por amostra:\n  ",
  sample_qc_output,
  "\n\n",
  sep = ""
)

cat(
  "Auditoria de features:\n  ",
  feature_audit_output,
  "\n\n",
  sep = ""
)

cat(
  "Features técnicas removidas:\n  ",
  technical_features_output,
  "\n\n",
  sep = ""
)

cat(
  "IDs gênicos não classificados pelo regex Ensembl:\n  ",
  nonstandard_feature_output,
  "\n\n",
  sep = ""
)

cat(
  "Dimensões:\n  ",
  dimensions_output,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 10. Conclusão
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "10. CONCLUSÃO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "- Foram selecionadas somente as 144 amostras DLPFC da análise principal.\n"
)

cat(
  "- A ordem das amostras é idêntica à ordem definida no metadata.\n"
)

cat(
  "- Foram removidas apenas as linhas técnicas N_*.\n"
)

cat(
  "- Os IDs Ensembl versionados foram preservados sem alteração.\n"
)

cat(
  "- Não foram detectados NA, valores negativos, não finitos ou fracionários.\n"
)

cat(
  "- Genes all-zero foram apenas auditados e permanecem na matriz.\n"
)

cat(
  "- Nenhum filtro de baixa expressão foi aplicado nesta etapa.\n"
)

cat(
  "- A matriz está pronta para 03_harmonize_genes_ROSMAP.R.\n\n"
)

cat(
  "============================================================\n"
)
cat(
  "02_prepare_counts_ROSMAP.R CONCLUÍDO\n"
)
cat(
  "============================================================\n"
)

sink()


# ==============================================================================
# 20. RESUMO NO TERMINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "COUNTS ROSMAP PREPARADOS\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Amostras: ",
  fmt_int(
    ncol(
      count_matrix
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  AD: ",
  fmt_int(
    n_ad
  ),
  "\n",
  sep = ""
)

cat(
  "  Control: ",
  fmt_int(
    n_control
  ),
  "\n\n",
  sep = ""
)

cat(
  "Features gênicas: ",
  fmt_int(
    nrow(
      count_matrix
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes all-zero: ",
  fmt_int(
    n_all_zero_genes
  ),
  "\n",
  sep = ""
)

cat(
  "Ensembl versionados: ",
  fmt_int(
    n_versioned_ensembl
  ),
  "\n",
  sep = ""
)

cat(
  "IDs não padrão: ",
  fmt_int(
    n_nonstandard
  ),
  "\n\n",
  sep = ""
)

cat(
  "Validação numérica:\n"
)

cat(
  "  NA = ",
  fmt_int(
    n_na
  ),
  "\n",
  sep = ""
)

cat(
  "  negativos = ",
  fmt_int(
    n_negative
  ),
  "\n",
  sep = ""
)

cat(
  "  fracionários = ",
  fmt_int(
    n_fractional
  ),
  "\n\n",
  sep = ""
)

cat(
  "Matriz preparada:\n  ",
  counts_output,
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
  "Próxima etapa, após revisar este relatório:\n"
)
cat(
  "  03_harmonize_genes_ROSMAP.R\n"
)
