#!/usr/bin/env Rscript

# ==============================================================================
# 02_prepare_counts_GSE125583.R
#
# FASE 01b — PREPARAÇÃO DA MATRIZ DE COUNTS DO DATASET CANDIDATO DE AD
# GEO: GSE125583
#
# OBJETIVO GERAL
# ------------------------------------------------------------------------------
# Construir a matriz de counts canônica que será utilizada nas etapas seguintes
# do pipeline, preservando os raw counts padronizados pelo NCBI sem aplicar
# qualquer transformação biológica ou estatística nesta etapa.
#
# A entrada principal é:
#
#   GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz
#
# A etapa 00 demonstrou que essa matriz possui:
#   - 39.376 GeneIDs;
#   - 289 amostras;
#   - raw counts inteiros;
#   - ausência de valores negativos, fracionários ou ausentes;
#   - o mesmo universo bruto de GeneIDs de GSE64810 (HD) e GSE68719 (PD).
#
# A etapa 01 produziu um metadata canônico com:
#   - 219 amostras AD;
#   - 70 controles;
#   - correspondência completa entre sample_id e as colunas da matriz.
#
#
# O QUE ESTE SCRIPT FAZ
# ------------------------------------------------------------------------------
# 1. Lê a matriz de raw counts do NCBI.
# 2. Lê o metadata preparado na etapa 01.
# 3. Valida novamente:
#      - existência da coluna GeneID;
#      - GeneIDs ausentes ou duplicados;
#      - sample IDs duplicados;
#      - correspondência counts <-> metadata;
#      - counts numéricos;
#      - ausência de NA, Inf, valores negativos e valores fracionários.
# 4. Reordena as colunas das amostras para seguir exatamente a ordem do metadata.
# 5. Confirma que nenhuma amostra e nenhum gene foram removidos.
# 6. Confirma que nenhum count foi alterado.
# 7. Salva uma matriz canônica:
#
#      GeneID | GSM... | GSM... | ...
#
# 8. Gera um relatório textual com auditorias antes e depois da preparação.
#
#
# O QUE ESTE SCRIPT NÃO FAZ
# ------------------------------------------------------------------------------
# - NÃO remove genes totalmente zerados.
# - NÃO aplica o filtro count >= 10.
# - NÃO normaliza por library size.
# - NÃO calcula CPM, TPM ou RPKM.
# - NÃO aplica VST/rlog.
# - NÃO executa DESeq2.
# - NÃO remove amostras.
# - NÃO arredonda counts.
# - NÃO altera GeneIDs.
#
# Os genes de baixa expressão e genes totalmente zerados serão tratados apenas
# na etapa de QC/pré-filtro, para manter uma separação clara entre:
#
#   preparação da matriz  -> etapa 02
#   harmonização/auditoria -> etapa 03
#   filtro para DE/QC      -> etapa 04
#
#
# ENTRADAS
# ------------------------------------------------------------------------------
# data/raw/GSE125583/GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz
#
# results/candidate_GSE125583/metadata/metadata_GSE125583.tsv
#
#
# SAÍDAS
# ------------------------------------------------------------------------------
# results/candidate_GSE125583/counts/counts_GSE125583.tsv
#
# results/candidate_GSE125583/counts/02_prepare_counts_GSE125583.txt
# ==============================================================================


# ==============================================================================
# 1. PACOTES
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
# 2. CAMINHOS
# ==============================================================================

raw_counts_file <- file.path(
  "data",
  "raw",
  "GSE125583",
  "GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)

metadata_file <- file.path(
  "results",
  "candidate_GSE125583",
  "metadata",
  "metadata_GSE125583.tsv"
)

output_dir <- file.path(
  "results",
  "candidate_GSE125583",
  "counts"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

prepared_counts_file <- file.path(
  output_dir,
  "counts_GSE125583.tsv"
)

audit_file <- file.path(
  output_dir,
  "02_prepare_counts_GSE125583.txt"
)


# ==============================================================================
# 3. CONSTANTES ESPERADAS
# ==============================================================================
#
# Estes valores foram observados e validados nas etapas 00 e 01.
# Eles servem como "travas de segurança": se a entrada mudar silenciosamente,
# o script interrompe a execução.

EXPECTED_N_GENES <- 39376L
EXPECTED_N_SAMPLES <- 289L
EXPECTED_N_AD <- 219L
EXPECTED_N_CONTROL <- 70L


# ==============================================================================
# 4. FUNÇÕES AUXILIARES
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 Validar existência do arquivo
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
  dt,
  dataset_name
) {

  if (!"GeneID" %in% names(dt)) {
    stop(
      dataset_name,
      ": coluna GeneID não encontrada."
    )
  }

  if (names(dt)[1] != "GeneID") {
    stop(
      dataset_name,
      ": GeneID deve ser a primeira coluna."
    )
  }

  gene_id <- as.character(
    dt$GeneID
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
    "group"
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
      "Existem sample_ids ausentes no metadata."
    )
  }

  if (anyDuplicated(metadata$sample_id)) {
    stop(
      "Existem sample_ids duplicados no metadata."
    )
  }

  if (any(is.na(metadata$group))) {
    stop(
      "Existem grupos ausentes no metadata."
    )
  }

  unexpected_groups <- setdiff(
    unique(
      metadata$group
    ),
    c(
      "AD",
      "Control"
    )
  )

  if (length(unexpected_groups) > 0) {
    stop(
      "Grupos inesperados no metadata: ",
      paste(
        unexpected_groups,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 4.4 Validar correspondência counts <-> metadata
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

  if (anyDuplicated(count_samples)) {
    stop(
      dataset_name,
      ": existem sample IDs duplicados na matriz de counts."
    )
  }

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
      ": amostras presentes apenas nos counts: ",
      paste(
        only_counts,
        collapse = ", "
      )
    )
  }

  if (length(only_metadata) > 0) {
    stop(
      dataset_name,
      ": amostras presentes apenas no metadata: ",
      paste(
        only_metadata,
        collapse = ", "
      )
    )
  }

  if (length(count_samples) != length(metadata_samples)) {
    stop(
      dataset_name,
      ": número diferente de amostras em counts e metadata."
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 4.5 Validar integridade dos counts
# ------------------------------------------------------------------------------

audit_count_matrix <- function(
  counts,
  dataset_name
) {

  sample_columns <- setdiff(
    names(counts),
    "GeneID"
  )

  non_numeric_columns <- sample_columns[
    !vapply(
      counts[
        ,
        ..sample_columns
      ],
      is.numeric,
      logical(1)
    )
  ]

  if (length(non_numeric_columns) > 0) {
    stop(
      dataset_name,
      ": colunas de counts não numéricas: ",
      paste(
        non_numeric_columns,
        collapse = ", "
      )
    )
  }

  count_matrix <- as.matrix(
    counts[
      ,
      ..sample_columns
    ]
  )

  n_missing <- sum(
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

  fractional_mask <- (
    !is.na(count_matrix) &
      is.finite(count_matrix) &
      abs(
        count_matrix -
          round(
            count_matrix
          )
      ) > 1e-8
  )

  n_fractional <- sum(
    fractional_mask
  )

  if (n_missing > 0) {
    stop(
      dataset_name,
      ": existem ",
      n_missing,
      " counts ausentes."
    )
  }

  if (n_nonfinite > 0) {
    stop(
      dataset_name,
      ": existem ",
      n_nonfinite,
      " counts não finitos."
    )
  }

  if (n_negative > 0) {
    stop(
      dataset_name,
      ": existem ",
      n_negative,
      " counts negativos."
    )
  }

  if (n_fractional > 0) {
    stop(
      dataset_name,
      ": existem ",
      n_fractional,
      " counts fracionários."
    )
  }

  library_sizes <- colSums(
    count_matrix
  )

  genes_detected_per_sample <- colSums(
    count_matrix > 0
  )

  row_totals <- rowSums(
    count_matrix
  )

  list(
    n_missing = n_missing,
    n_nonfinite = n_nonfinite,
    n_negative = n_negative,
    n_fractional = n_fractional,

    n_zero_counts = sum(
      count_matrix == 0
    ),

    zero_fraction = mean(
      count_matrix == 0
    ),

    all_zero_genes = sum(
      row_totals == 0
    ),

    all_zero_samples = sum(
      library_sizes == 0
    ),

    total_counts = sum(
      count_matrix
    ),

    library_sizes = library_sizes,

    genes_detected_per_sample = genes_detected_per_sample
  )
}


# ------------------------------------------------------------------------------
# 4.6 Formatar inteiros
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


# ------------------------------------------------------------------------------
# 4.7 Formatar números decimais
# ------------------------------------------------------------------------------

format_num <- function(
  x,
  digits = 4
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


# ==============================================================================
# 5. VALIDAR ARQUIVOS DE ENTRADA
# ==============================================================================

validate_file(
  raw_counts_file
)

validate_file(
  metadata_file
)


# ==============================================================================
# 6. LER METADATA PREPARADO NA ETAPA 01
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
    "Número inesperado de amostras AD no metadata. Esperado: ",
    EXPECTED_N_AD,
    "; observado: ",
    n_ad,
    "."
  )
}

if (n_control != EXPECTED_N_CONTROL) {
  stop(
    "Número inesperado de controles no metadata. Esperado: ",
    EXPECTED_N_CONTROL,
    "; observado: ",
    n_control,
    "."
  )
}


# ==============================================================================
# 7. LER MATRIZ BRUTA DO NCBI
# ==============================================================================

counts_raw <- fread(
  raw_counts_file,
  showProgress = TRUE
)

# GeneID será tratado explicitamente como identificador textual.
counts_raw[
  ,
  GeneID := as.character(
    GeneID
  )
]

validate_gene_id(
  counts_raw,
  "GSE125583"
)

if (nrow(counts_raw) != EXPECTED_N_GENES) {
  stop(
    "Número inesperado de genes. Esperado: ",
    EXPECTED_N_GENES,
    "; observado: ",
    nrow(counts_raw),
    "."
  )
}

raw_sample_ids <- setdiff(
  names(counts_raw),
  "GeneID"
)

if (length(raw_sample_ids) != EXPECTED_N_SAMPLES) {
  stop(
    "Número inesperado de amostras na matriz. Esperado: ",
    EXPECTED_N_SAMPLES,
    "; observado: ",
    length(raw_sample_ids),
    "."
  )
}


# ==============================================================================
# 8. VALIDAR CORRESPONDÊNCIA COUNTS <-> METADATA
# ==============================================================================

validate_sample_mapping(
  counts = counts_raw,
  metadata = metadata,
  dataset_name = "GSE125583"
)

same_order_before <- identical(
  raw_sample_ids,
  metadata$sample_id
)


# ==============================================================================
# 9. AUDITAR COUNTS ANTES DA PREPARAÇÃO
# ==============================================================================

audit_before <- audit_count_matrix(
  counts_raw,
  "GSE125583 — entrada"
)


# ==============================================================================
# 10. CONSTRUIR MATRIZ CANÔNICA
# ==============================================================================
#
# A única operação permitida nesta etapa é organizar as colunas:
#
#   GeneID + amostras na ordem exata do metadata.
#
# Nenhum count é transformado.

canonical_columns <- c(
  "GeneID",
  metadata$sample_id
)

counts_prepared <- counts_raw[
  ,
  ..canonical_columns
]


# ==============================================================================
# 11. VALIDAR A MATRIZ PREPARADA
# ==============================================================================

validate_gene_id(
  counts_prepared,
  "GSE125583 preparado"
)

validate_sample_mapping(
  counts = counts_prepared,
  metadata = metadata,
  dataset_name = "GSE125583 preparado"
)

prepared_sample_ids <- setdiff(
  names(counts_prepared),
  "GeneID"
)

if (!identical(
  prepared_sample_ids,
  metadata$sample_id
)) {
  stop(
    "A ordem das colunas da matriz preparada não coincide com o metadata."
  )
}

if (nrow(counts_prepared) != nrow(counts_raw)) {
  stop(
    "O número de genes mudou durante a preparação."
  )
}

if (ncol(counts_prepared) != ncol(counts_raw)) {
  stop(
    "O número de colunas mudou durante a preparação."
  )
}

if (!identical(
  counts_prepared$GeneID,
  counts_raw$GeneID
)) {
  stop(
    "A ordem ou o conteúdo de GeneID mudou durante a preparação."
  )
}


# ==============================================================================
# 12. CONFIRMAR QUE NENHUM COUNT FOI ALTERADO
# ==============================================================================
#
# Como a preparação pode reordenar colunas, comparamos cada amostra pelo nome.
# O conteúdo de cada coluna deve permanecer exatamente igual ao arquivo bruto.

count_values_preserved <- all(
  vapply(
    metadata$sample_id,
    function(sample_id) {

      identical(
        counts_prepared[[sample_id]],
        counts_raw[[sample_id]]
      )
    },
    logical(1)
  )
)

if (!count_values_preserved) {
  stop(
    "Pelo menos uma coluna de counts foi alterada durante a preparação."
  )
}


# ==============================================================================
# 13. AUDITAR COUNTS APÓS A PREPARAÇÃO
# ==============================================================================

audit_after <- audit_count_matrix(
  counts_prepared,
  "GSE125583 — preparado"
)


# ==============================================================================
# 14. VALIDAÇÕES DE CONSERVAÇÃO
# ==============================================================================
#
# Todas estas métricas devem ser idênticas antes e depois da preparação.

if (
  audit_before$total_counts !=
    audit_after$total_counts
) {
  stop(
    "A soma total dos counts mudou durante a preparação."
  )
}

if (
  audit_before$n_zero_counts !=
    audit_after$n_zero_counts
) {
  stop(
    "O número de zeros mudou durante a preparação."
  )
}

if (
  audit_before$all_zero_genes !=
    audit_after$all_zero_genes
) {
  stop(
    "O número de genes totalmente zerados mudou durante a preparação."
  )
}

if (
  !identical(
    sort(
      audit_before$library_sizes
    ),
    sort(
      audit_after$library_sizes
    )
  )
) {
  stop(
    "Library sizes mudaram durante a preparação."
  )
}


# ==============================================================================
# 15. SALVAR MATRIZ CANÔNICA
# ==============================================================================

fwrite(
  counts_prepared,
  prepared_counts_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 16. REABRIR CABEÇALHO DO ARQUIVO SALVO
# ==============================================================================
#
# Esta verificação final confirma que o arquivo escrito em disco possui a
# estrutura esperada.

saved_header <- fread(
  prepared_counts_file,
  nrows = 0,
  showProgress = FALSE
)

expected_saved_columns <- c(
  "GeneID",
  metadata$sample_id
)

if (!identical(
  names(saved_header),
  expected_saved_columns
)) {
  stop(
    "O arquivo salvo não possui as colunas esperadas."
  )
}


# ==============================================================================
# 17. GERAR RELATÓRIO TEXTUAL
# ==============================================================================

sink(
  audit_file
)

cat("============================================================\n")
cat("02 — PREPARAÇÃO DOS COUNTS GSE125583\n")
cat("============================================================\n\n")


cat("------------------------------------------------------------\n")
cat("1. ENTRADA\n")
cat("------------------------------------------------------------\n")

cat(
  "Arquivo bruto: ",
  raw_counts_file,
  "\n",
  sep = ""
)

cat(
  "Genes: ",
  format_int(
    nrow(
      counts_raw
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras: ",
  format_int(
    length(
      raw_sample_ids
    )
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
cat("2. VALIDAÇÃO COUNTS <-> METADATA\n")
cat("------------------------------------------------------------\n")

cat(
  "Mesmas amostras em counts e metadata: TRUE\n"
)

cat(
  "Mesma ordem antes da preparação: ",
  same_order_before,
  "\n",
  sep = ""
)

cat(
  "Mesma ordem após a preparação: ",
  identical(
    prepared_sample_ids,
    metadata$sample_id
  ),
  "\n",
  sep = ""
)

cat(
  "Counts preservados exatamente: ",
  count_values_preserved,
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("3. INTEGRIDADE DOS RAW COUNTS\n")
cat("------------------------------------------------------------\n")

cat(
  "Counts ausentes: ",
  format_int(
    audit_before$n_missing
  ),
  "\n",
  sep = ""
)

cat(
  "Counts não finitos: ",
  format_int(
    audit_before$n_nonfinite
  ),
  "\n",
  sep = ""
)

cat(
  "Counts negativos: ",
  format_int(
    audit_before$n_negative
  ),
  "\n",
  sep = ""
)

cat(
  "Counts fracionários: ",
  format_int(
    audit_before$n_fractional
  ),
  "\n",
  sep = ""
)

cat(
  "Genes totalmente zerados: ",
  format_int(
    audit_before$all_zero_genes
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras totalmente zeradas: ",
  format_int(
    audit_before$all_zero_samples
  ),
  "\n",
  sep = ""
)

cat(
  "Proporção global de zeros: ",
  format_num(
    100 *
      audit_before$zero_fraction,
    digits = 2
  ),
  "%\n",
  sep = ""
)

cat(
  "Soma total dos counts: ",
  format_int(
    audit_before$total_counts
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("4. LIBRARY SIZE\n")
cat("------------------------------------------------------------\n")

cat(
  "Mínimo: ",
  format_int(
    min(
      audit_before$library_sizes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Mediana: ",
  format_int(
    median(
      audit_before$library_sizes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Máximo: ",
  format_int(
    max(
      audit_before$library_sizes
    )
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("5. GENES DETECTADOS POR AMOSTRA\n")
cat("------------------------------------------------------------\n")

cat(
  "Mínimo: ",
  format_int(
    min(
      audit_before$genes_detected_per_sample
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Mediana: ",
  format_int(
    median(
      audit_before$genes_detected_per_sample
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Máximo: ",
  format_int(
    max(
      audit_before$genes_detected_per_sample
    )
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("6. CONSERVAÇÃO DA MATRIZ\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes antes: ",
  format_int(
    nrow(
      counts_raw
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes depois: ",
  format_int(
    nrow(
      counts_prepared
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras antes: ",
  format_int(
    length(
      raw_sample_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras depois: ",
  format_int(
    length(
      prepared_sample_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes removidos: ",
  format_int(
    nrow(
      counts_raw
    ) -
      nrow(
        counts_prepared
      )
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras removidas: ",
  format_int(
    length(
      raw_sample_ids
    ) -
      length(
        prepared_sample_ids
      )
  ),
  "\n",
  sep = ""
)

cat(
  "Soma total antes: ",
  format_int(
    audit_before$total_counts
  ),
  "\n",
  sep = ""
)

cat(
  "Soma total depois: ",
  format_int(
    audit_after$total_counts
  ),
  "\n",
  sep = ""
)

cat(
  "Zeros antes: ",
  format_int(
    audit_before$n_zero_counts
  ),
  "\n",
  sep = ""
)

cat(
  "Zeros depois: ",
  format_int(
    audit_after$n_zero_counts
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("7. SAÍDA\n")
cat("------------------------------------------------------------\n")

cat(
  "Matriz preparada: ",
  prepared_counts_file,
  "\n",
  sep = ""
)

cat(
  "Estrutura: GeneID + ",
  format_int(
    length(
      prepared_sample_ids
    )
  ),
  " amostras\n",
  sep = ""
)

cat(
  "Ordem das amostras = metadata: TRUE\n\n"
)


cat("------------------------------------------------------------\n")
cat("8. INTERPRETAÇÃO DA ETAPA\n")
cat("------------------------------------------------------------\n")

cat(
  "- A matriz NCBI foi preservada sem alteração dos valores de counts.\n"
)

cat(
  "- Nenhum gene foi removido nesta etapa.\n"
)

cat(
  "- Nenhuma amostra foi removida nesta etapa.\n"
)

cat(
  "- Nenhuma normalização foi aplicada.\n"
)

cat(
  "- Nenhum filtro de baixa expressão foi aplicado.\n"
)

cat(
  "- Genes totalmente zerados foram deliberadamente preservados para que a remoção ocorra formalmente na etapa de QC/pré-filtro.\n"
)

cat(
  "- A única operação estrutural permitida foi alinhar a ordem das amostras ao metadata.\n"
)

cat("\n")
cat("============================================================\n")
cat("02_prepare_counts_GSE125583.R CONCLUÍDO\n")
cat("============================================================\n")

sink()


# ==============================================================================
# 18. RESUMO NO TERMINAL
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("02 — COUNTS GSE125583 CONCLUÍDO\n")
cat("============================================================\n")

cat(
  "Genes: ",
  format_int(
    nrow(
      counts_prepared
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras: ",
  format_int(
    length(
      prepared_sample_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes removidos: 0\n"
)

cat(
  "Amostras removidas: 0\n"
)

cat(
  "Genes totalmente zerados preservados: ",
  format_int(
    audit_after$all_zero_genes
  ),
  "\n",
  sep = ""
)

cat(
  "Counts preservados exatamente: ",
  count_values_preserved,
  "\n",
  sep = ""
)

cat(
  "Counts <-> metadata alinhados: ",
  identical(
    prepared_sample_ids,
    metadata$sample_id
  ),
  "\n",
  sep = ""
)

cat(
  "Matriz preparada: ",
  prepared_counts_file,
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
