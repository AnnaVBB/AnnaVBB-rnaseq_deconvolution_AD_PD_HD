# ==============================================================================
# 04_qc_counts_ROSMAP.R
# ==============================================================================
#
# Fase 01c — ROSMAP
# Controle de qualidade, filtragem de baixa expressão e PCA exploratória
#
# Este script deve ser executado APÓS:
#
#   00_audit_ROSMAP.R
#   01_prepare_metadata_ROSMAP.R
#   02_prepare_counts_ROSMAP.R
#   03_harmonize_genes_ROSMAP.R
#
# Entradas principais
# -------------------
#
# Matriz harmonizada por NCBI GeneID:
#
#   results/candidate_ROSMAP/harmonization/
#   counts_ROSMAP_GeneID_main.tsv.gz
#
# Metadata principal:
#
#   results/candidate_ROSMAP/metadata/
#   metadata_ROSMAP.tsv
#
# Anotação gênica:
#
#   results/candidate_ROSMAP/harmonization/
#   gene_annotation_ROSMAP_GeneID.tsv
#
# Objetivos
# ---------
# 1. Validar novamente a correspondência 1:1 entre counts e metadata.
#
# 2. Auditar a matriz harmonizada antes da filtragem:
#
#       - GeneID duplicados/ausentes;
#       - NA;
#       - valores não finitos;
#       - valores negativos;
#       - valores fracionários;
#       - genes all-zero;
#       - amostras all-zero;
#       - library size;
#       - genes detectados por amostra.
#
# 3. Aplicar o filtro de baixa expressão PRÉ-DEFINIDO:
#
#       count >= 10 em pelo menos 67 amostras
#
#    O número 67 corresponde ao tamanho do menor grupo experimental:
#
#       Control = 67
#       AD      = 77
#
# 4. Salvar a matriz filtrada de RAW COUNTS que será usada no DESeq2.
#
# 5. Quantificar o impacto do filtro:
#
#       - genes mantidos/removidos;
#       - counts mantidos/removidos;
#       - métricas de detecção por amostra antes/depois.
#
# 6. Realizar transformação VST SOMENTE para controle de qualidade:
#
#       varianceStabilizingTransformation / vst
#
#    IMPORTANTE:
#    A VST NÃO será usada como entrada da análise diferencial.
#    O DESeq2 será executado posteriormente sobre raw counts filtrados.
#
# 7. Produzir PCA exploratória não supervisionada (blind = TRUE), usando os
#    500 genes de maior variância após VST.
#
# 8. Avaliar visualmente:
#
#       - group (AD vs Control);
#       - sequencingBatch.
#
# 9. NÃO excluir nenhuma amostra automaticamente com base em PCA, library size
#    ou genes detectados. Qualquer exclusão posterior exigiria evidência técnica
#    independente e documentação explícita.
#
# Saídas principais
# -----------------
#
#   results/candidate_ROSMAP/qc/
#
#     counts_ROSMAP_GeneID_filtered.tsv.gz
#     gene_filter_audit_ROSMAP.tsv
#     gene_filter_summary_ROSMAP.tsv
#     sample_qc_ROSMAP.tsv
#     qc_group_tests_ROSMAP.tsv
#     size_factors_ROSMAP.tsv
#     pca_scores_ROSMAP.tsv
#     pca_variance_ROSMAP.tsv
#     vst_ROSMAP_filtered.rds
#     gene_annotation_ROSMAP_filtered.tsv
#     04_qc_counts_ROSMAP.txt
#
#     plots/
#       library_size_by_group_ROSMAP.png
#       detected_genes_by_group_ROSMAP.png
#       PCA_ROSMAP_group.png
#       PCA_ROSMAP_sequencingBatch.png
#
# ==============================================================================


# ==============================================================================
# 0. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table",
  "DESeq2",
  "SummarizedExperiment",
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

if (
  length(
    missing_packages
  ) > 0
) {

  stop(
    paste0(
      "Pacotes necessários não instalados: ",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n\n",
      "Este script não instala pacotes automaticamente. ",
      "Use o ambiente do projeto antes de continuar."
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
})


# ==============================================================================
# 1. PARÂMETROS PRÉ-DEFINIDOS
# ==============================================================================

# ------------------------------------------------------------------------------
# Tamanho dos grupos
# ------------------------------------------------------------------------------

EXPECTED_MAIN_N <- 144L
EXPECTED_AD_N <- 77L
EXPECTED_CONTROL_N <- 67L


# ------------------------------------------------------------------------------
# Filtro de baixa expressão
#
# Manter GeneID se:
#
#   count >= MIN_COUNT
#
# em pelo menos:
#
#   MIN_SAMPLES
#
# amostras.
#
# MIN_SAMPLES = 67 corresponde ao menor grupo experimental.
# ------------------------------------------------------------------------------

MIN_COUNT <- 10L
MIN_SAMPLES <- EXPECTED_CONTROL_N


# ------------------------------------------------------------------------------
# PCA
#
# Segue a lógica usual de PCA exploratória com genes de maior variância.
# ------------------------------------------------------------------------------

PCA_NTOP <- 500L


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

candidate_dir <- file.path(
  "results",
  "candidate_ROSMAP"
)

harmonization_dir <- file.path(
  candidate_dir,
  "harmonization"
)

metadata_dir <- file.path(
  candidate_dir,
  "metadata"
)

qc_dir <- file.path(
  candidate_dir,
  "qc"
)

plot_dir <- file.path(
  qc_dir,
  "plots"
)

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Entradas
# ------------------------------------------------------------------------------

counts_input <- file.path(
  harmonization_dir,
  "counts_ROSMAP_GeneID_main.tsv.gz"
)

metadata_input <- file.path(
  metadata_dir,
  "metadata_ROSMAP.tsv"
)

annotation_input <- file.path(
  harmonization_dir,
  "gene_annotation_ROSMAP_GeneID.tsv"
)


# ------------------------------------------------------------------------------
# Saídas
# ------------------------------------------------------------------------------

counts_filtered_output <- file.path(
  qc_dir,
  "counts_ROSMAP_GeneID_filtered.tsv.gz"
)

gene_filter_audit_output <- file.path(
  qc_dir,
  "gene_filter_audit_ROSMAP.tsv"
)

gene_filter_summary_output <- file.path(
  qc_dir,
  "gene_filter_summary_ROSMAP.tsv"
)

sample_qc_output <- file.path(
  qc_dir,
  "sample_qc_ROSMAP.tsv"
)

group_tests_output <- file.path(
  qc_dir,
  "qc_group_tests_ROSMAP.tsv"
)

size_factors_output <- file.path(
  qc_dir,
  "size_factors_ROSMAP.tsv"
)

pca_scores_output <- file.path(
  qc_dir,
  "pca_scores_ROSMAP.tsv"
)

pca_variance_output <- file.path(
  qc_dir,
  "pca_variance_ROSMAP.tsv"
)

vst_output <- file.path(
  qc_dir,
  "vst_ROSMAP_filtered.rds"
)

annotation_filtered_output <- file.path(
  qc_dir,
  "gene_annotation_ROSMAP_filtered.tsv"
)

report_output <- file.path(
  qc_dir,
  "04_qc_counts_ROSMAP.txt"
)


# ------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------

library_size_plot_output <- file.path(
  plot_dir,
  "library_size_by_group_ROSMAP.png"
)

detected_genes_plot_output <- file.path(
  plot_dir,
  "detected_genes_by_group_ROSMAP.png"
)

pca_group_plot_output <- file.path(
  plot_dir,
  "PCA_ROSMAP_group.png"
)

pca_batch_plot_output <- file.path(
  plot_dir,
  "PCA_ROSMAP_sequencingBatch.png"
)


# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================

assert_file <- function(path) {

  if (
    !file.exists(
      path
    )
  ) {

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
    names(
      dt
    )
  )

  if (
    length(
      missing
    ) > 0
  ) {

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
    round(
      as.numeric(
        x
      ),
      digits = 0
    ),
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


fmt_num <- function(
  x,
  digits = 4
) {

  format(
    round(
      as.numeric(
        x
      ),
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


safe_wilcox <- function(
  x,
  group
) {

  ok <- !is.na(
    x
  ) &
    !is.na(
      group
    )

  x <- x[
    ok
  ]

  group <- group[
    ok
  ]

  if (
    length(
      unique(
        group
      )
    ) != 2
  ) {

    return(
      NA_real_
    )
  }

  suppressWarnings(
    wilcox.test(
      x ~ group,
      exact = FALSE
    )$p.value
  )
}


# ------------------------------------------------------------------------------
# Salvar TSV comprimido usando gzip do sistema.
#
# Isso evita dependência do suporte zlib na compilação local de data.table.
# ------------------------------------------------------------------------------

write_tsv_gzip <- function(
  dt,
  path_gz
) {

  if (
    !grepl(
      "\\.gz$",
      path_gz
    )
  ) {

    stop(
      "write_tsv_gzip exige um caminho terminado em .gz."
    )
  }

  path_plain <- sub(
    "\\.gz$",
    "",
    path_gz
  )

  fwrite(
    dt,
    path_plain,
    sep = "\t",
    quote = FALSE,
    na = "NA"
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
        "gzip não encontrado no sistema. ",
        "Arquivo sem compressão salvo em:\n",
        path_plain
      )
    )
  }

  status <- system2(
    command = gzip_bin,
    args = c(
      "-f",
      path_plain
    )
  )

  if (
    status != 0
  ) {

    stop(
      "Falha ao comprimir com gzip. Código de saída: ",
      status
    )
  }

  if (
    !file.exists(
      path_gz
    )
  ) {

    stop(
      "Arquivo comprimido esperado não foi encontrado: ",
      path_gz
    )
  }

  invisible(
    path_gz
  )
}


# ==============================================================================
# 4. VALIDAR ARQUIVOS
# ==============================================================================

assert_file(
  counts_input
)

assert_file(
  metadata_input
)

assert_file(
  annotation_input
)


# ==============================================================================
# 5. LER METADATA
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "04 — QC E FILTRAGEM DOS COUNTS ROSMAP\n"
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
  sample_id := trimws(
    as.character(
      sample_id
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

metadata[
  ,
  sequencingBatch := trimws(
    as.character(
      sequencingBatch
    )
  )
]


# ------------------------------------------------------------------------------
# Validar metadata
# ------------------------------------------------------------------------------

if (
  nrow(
    metadata
  ) != EXPECTED_MAIN_N
) {

  stop(
    "Metadata principal deveria conter ",
    EXPECTED_MAIN_N,
    " amostras; observado: ",
    nrow(
      metadata
    ),
    "."
  )
}

if (
  anyDuplicated(
    metadata$sample_id
  )
) {

  stop(
    "metadata contém sample_id duplicado."
  )
}

if (
  anyDuplicated(
    metadata$individual_id
  )
) {

  stop(
    "metadata contém individual_id duplicado."
  )
}

group_counts <- metadata[
  ,
  .N,
  by = group
]

n_ad <- group_counts[
  group == "AD",
  N
]

n_control <- group_counts[
  group == "Control",
  N
]

if (
  length(
    n_ad
  ) != 1 ||
    n_ad != EXPECTED_AD_N
) {

  stop(
    "Número inesperado de amostras AD."
  )
}

if (
  length(
    n_control
  ) != 1 ||
    n_control != EXPECTED_CONTROL_N
) {

  stop(
    "Número inesperado de controles."
  )
}

cat(
  "Metadata validado: ",
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
# 6. LER MATRIZ HARMONIZADA
# ==============================================================================

cat(
  "Lendo matriz harmonizada por GeneID...\n"
)

counts <- fread(
  counts_input,
  na.strings = c(
    "NA",
    ""
  )
)

assert_columns(
  counts,
  "GeneID",
  "counts_ROSMAP_GeneID_main.tsv.gz"
)

counts[
  ,
  GeneID := trimws(
    as.character(
      GeneID
    )
  )
]

count_sample_ids <- names(
  counts
)[
  -1
]


# ------------------------------------------------------------------------------
# Correspondência counts ↔ metadata
# ------------------------------------------------------------------------------

if (
  length(
    count_sample_ids
  ) != EXPECTED_MAIN_N
) {

  stop(
    "A matriz harmonizada não possui 144 amostras."
  )
}

if (
  !identical(
    count_sample_ids,
    metadata$sample_id
  )
) {

  stop(
    paste0(
      "A ordem das colunas da matriz harmonizada não corresponde ",
      "exatamente à ordem do metadata."
    )
  )
}

if (
  anyNA(
    counts$GeneID
  ) ||
    any(
      counts$GeneID == ""
    )
) {

  stop(
    "Existem GeneID ausentes/vazios."
  )
}

if (
  anyDuplicated(
    counts$GeneID
  )
) {

  stop(
    "Existem GeneID duplicados na matriz harmonizada."
  )
}

n_genes_before <- nrow(
  counts
)

cat(
  "Matriz harmonizada: ",
  fmt_int(
    n_genes_before
  ),
  " GeneIDs × ",
  fmt_int(
    length(
      count_sample_ids
    )
  ),
  " amostras.\n\n",
  sep = ""
)


# ==============================================================================
# 7. CONVERTER PARA MATRIZ NUMÉRICA E VALIDAR
# ==============================================================================

count_matrix <- as.matrix(
  counts[
    ,
    -1,
    with = FALSE
  ]
)

storage.mode(
  count_matrix
) <- "numeric"

rownames(
  count_matrix
) <- counts$GeneID


# ------------------------------------------------------------------------------
# Integridade antes do filtro
# ------------------------------------------------------------------------------

n_na_before <- sum(
  is.na(
    count_matrix
  )
)

n_nonfinite_before <- sum(
  !is.finite(
    count_matrix
  ),
  na.rm = TRUE
)

n_negative_before <- sum(
  count_matrix < 0,
  na.rm = TRUE
)

n_fractional_before <- count_fractional_values(
  count_matrix
)

n_all_zero_genes_before <- sum(
  rowSums(
    count_matrix
  ) == 0
)

n_all_zero_samples_before <- sum(
  colSums(
    count_matrix
  ) == 0
)

if (
  n_na_before > 0
) {

  stop(
    "Foram encontrados NA antes da filtragem."
  )
}

if (
  n_nonfinite_before > 0
) {

  stop(
    "Foram encontrados valores não finitos antes da filtragem."
  )
}

if (
  n_negative_before > 0
) {

  stop(
    "Foram encontrados counts negativos antes da filtragem."
  )
}

if (
  n_fractional_before > 0
) {

  stop(
    "Foram encontrados counts fracionários antes da filtragem."
  )
}

if (
  n_all_zero_samples_before > 0
) {

  stop(
    "Existe amostra all-zero antes da filtragem."
  )
}


# ==============================================================================
# 8. MÉTRICAS PRÉ-FILTRO
# ==============================================================================

library_size_before <- colSums(
  count_matrix
)

detected_gt0_before <- colSums(
  count_matrix > 0
)

detected_ge10_before <- colSums(
  count_matrix >= MIN_COUNT
)

total_counts_before <- sum(
  count_matrix
)


# ==============================================================================
# 9. FILTRO DE BAIXA EXPRESSÃO
# ==============================================================================
#
# Critério congelado:
#
#   count >= 10 em pelo menos 67 amostras
#
# Não depende do diagnóstico de cada gene.
# O filtro é aplicado de forma idêntica a toda a matriz.
# ==============================================================================

cat(
  "Aplicando filtro de baixa expressão:\n"
)

cat(
  "  count >= ",
  MIN_COUNT,
  " em pelo menos ",
  MIN_SAMPLES,
  " amostras.\n\n",
  sep = ""
)

n_samples_ge_min_count <- rowSums(
  count_matrix >= MIN_COUNT
)

keep_gene <- n_samples_ge_min_count >= MIN_SAMPLES

n_genes_after <- sum(
  keep_gene
)

n_genes_removed <- sum(
  !keep_gene
)

if (
  n_genes_after == 0
) {

  stop(
    "O filtro removeu todos os genes. Verifique os parâmetros."
  )
}


# ------------------------------------------------------------------------------
# Auditoria gene a gene
# ------------------------------------------------------------------------------

gene_filter_audit <- data.table(
  GeneID = rownames(
    count_matrix
  ),
  n_samples_count_ge_10 = as.integer(
    n_samples_ge_min_count
  ),
  total_count = as.numeric(
    rowSums(
      count_matrix
    )
  ),
  mean_count = as.numeric(
    rowMeans(
      count_matrix
    )
  ),
  median_count = as.numeric(
    apply(
      count_matrix,
      1,
      median
    )
  ),
  all_zero_before_filter = rowSums(
    count_matrix
  ) == 0,
  keep_low_expression_filter = keep_gene
)

fwrite(
  gene_filter_audit,
  gene_filter_audit_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 10. MATRIZ FILTRADA
# ==============================================================================

filtered_matrix <- count_matrix[
  keep_gene,
  ,
  drop = FALSE
]

n_all_zero_genes_after <- sum(
  rowSums(
    filtered_matrix
  ) == 0
)

n_all_zero_samples_after <- sum(
  colSums(
    filtered_matrix
  ) == 0
)

if (
  n_all_zero_genes_after > 0
) {

  stop(
    paste0(
      "Genes all-zero permaneceram após um filtro que exige count >= 10. ",
      "Isso indica inconsistência lógica."
    )
  )
}

if (
  n_all_zero_samples_after > 0
) {

  stop(
    "Existe amostra all-zero após a filtragem."
  )
}


# ------------------------------------------------------------------------------
# Integridade após filtro
# ------------------------------------------------------------------------------

n_na_after <- sum(
  is.na(
    filtered_matrix
  )
)

n_nonfinite_after <- sum(
  !is.finite(
    filtered_matrix
  ),
  na.rm = TRUE
)

n_negative_after <- sum(
  filtered_matrix < 0,
  na.rm = TRUE
)

n_fractional_after <- count_fractional_values(
  filtered_matrix
)

if (
  n_na_after > 0 ||
    n_nonfinite_after > 0 ||
    n_negative_after > 0 ||
    n_fractional_after > 0
) {

  stop(
    "Problema de integridade numérica após filtragem."
  )
}


# ==============================================================================
# 11. IMPACTO DO FILTRO SOBRE COUNTS
# ==============================================================================

total_counts_after <- sum(
  filtered_matrix
)

total_counts_removed <- (
  total_counts_before -
    total_counts_after
)

percent_counts_retained <- (
  100 *
    total_counts_after /
    total_counts_before
)

percent_counts_removed <- (
  100 -
    percent_counts_retained
)


# ==============================================================================
# 12. MÉTRICAS PÓS-FILTRO POR AMOSTRA
# ==============================================================================

library_size_after <- colSums(
  filtered_matrix
)

detected_gt0_after <- colSums(
  filtered_matrix > 0
)

detected_ge10_after <- colSums(
  filtered_matrix >= MIN_COUNT
)


# ==============================================================================
# 13. CONSTRUIR QC POR AMOSTRA
# ==============================================================================

sample_qc <- data.table(
  sample_id = metadata$sample_id,
  individual_id = metadata$individual_id,
  group = metadata$group,
  cogdx = metadata$cogdx,
  sequencingBatch = metadata$sequencingBatch,
  msex = metadata$msex,
  RIN = metadata$RIN,
  pmi = metadata$pmi,

  library_size_before_filter = as.numeric(
    library_size_before
  ),

  library_size_after_filter = as.numeric(
    library_size_after
  ),

  library_size_percent_retained = (
    100 *
      library_size_after /
      library_size_before
  ),

  detected_genes_gt0_before = as.integer(
    detected_gt0_before
  ),

  detected_genes_gt0_after = as.integer(
    detected_gt0_after
  ),

  detected_genes_ge10_before = as.integer(
    detected_ge10_before
  ),

  detected_genes_ge10_after = as.integer(
    detected_ge10_after
  )
)

fwrite(
  sample_qc,
  sample_qc_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 14. TESTES DESCRITIVOS ENTRE GRUPOS
# ==============================================================================
#
# Esses testes NÃO são usados para decidir exclusão de amostras.
# Servem apenas como caracterização de possíveis diferenças globais de
# profundidade/detecção entre AD e Control.
# ==============================================================================

library_wilcox_p <- safe_wilcox(
  sample_qc$library_size_before_filter,
  sample_qc$group
)

detected_wilcox_p <- safe_wilcox(
  sample_qc$detected_genes_gt0_before,
  sample_qc$group
)

library_after_wilcox_p <- safe_wilcox(
  sample_qc$library_size_after_filter,
  sample_qc$group
)

detected_after_wilcox_p <- safe_wilcox(
  sample_qc$detected_genes_gt0_after,
  sample_qc$group
)

group_tests <- data.table(
  metric = c(
    "library_size_before_filter",
    "detected_genes_gt0_before",
    "library_size_after_filter",
    "detected_genes_gt0_after"
  ),
  test = "Wilcoxon rank-sum",
  p_value = c(
    library_wilcox_p,
    detected_wilcox_p,
    library_after_wilcox_p,
    detected_after_wilcox_p
  )
)

fwrite(
  group_tests,
  group_tests_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 15. SALVAR MATRIZ FILTRADA
# ==============================================================================

filtered_counts_dt <- as.data.table(
  filtered_matrix,
  keep.rownames = "GeneID"
)

filtered_counts_dt[
  ,
  GeneID := as.character(
    GeneID
  )
]

if (
  anyDuplicated(
    filtered_counts_dt$GeneID
  )
) {

  stop(
    "GeneID duplicado na matriz filtrada."
  )
}

if (
  !identical(
    names(
      filtered_counts_dt
    )[
      -1
    ],
    metadata$sample_id
  )
) {

  stop(
    "Ordem counts ↔ metadata mudou após filtragem."
  )
}

cat(
  "Salvando matriz filtrada...\n"
)

write_tsv_gzip(
  filtered_counts_dt,
  counts_filtered_output
)


# ==============================================================================
# 16. ANOTAÇÃO FILTRADA
# ==============================================================================

annotation <- fread(
  annotation_input,
  na.strings = c(
    "NA",
    ""
  )
)

assert_columns(
  annotation,
  "GeneID",
  "gene_annotation_ROSMAP_GeneID.tsv"
)

annotation[
  ,
  GeneID := as.character(
    GeneID
  )
]

if (
  anyDuplicated(
    annotation$GeneID
  )
) {

  stop(
    "A tabela de anotação contém GeneID duplicado."
  )
}

annotation_filtered <- annotation[
  match(
    rownames(
      filtered_matrix
    ),
    GeneID
  )
]

if (
  anyNA(
    annotation_filtered$GeneID
  )
) {

  stop(
    "Há GeneID filtrado sem linha correspondente na tabela de anotação."
  )
}

if (
  !identical(
    annotation_filtered$GeneID,
    rownames(
      filtered_matrix
    )
  )
) {

  stop(
    "Ordem da anotação filtrada não corresponde à matriz filtrada."
  )
}

fwrite(
  annotation_filtered,
  annotation_filtered_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 17. PREPARAR DESeqDataSet APENAS PARA VST/QC
# ==============================================================================
#
# IMPORTANTE:
# Este objeto NÃO é o modelo de expressão diferencial.
#
# Usamos design = ~ 1 porque:
#
#   - a PCA é exploratória;
#   - blind = TRUE;
#   - queremos evitar usar diagnóstico para definir a transformação.
#
# O modelo diferencial real será construído no script 05.
# ==============================================================================

cat(
  "Preparando VST para PCA exploratória...\n"
)

# DESeq2 exige counts inteiros.
# Já comprovamos que não há valores fracionários, portanto a coerção para
# integer não modifica os valores.
filtered_matrix_integer <- filtered_matrix

storage.mode(
  filtered_matrix_integer
) <- "integer"

coldata_qc <- as.data.frame(
  metadata
)

rownames(
  coldata_qc
) <- metadata$sample_id

dds_qc <- DESeqDataSetFromMatrix(
  countData = filtered_matrix_integer,
  colData = coldata_qc,
  design = ~ 1
)

dds_qc <- estimateSizeFactors(
  dds_qc
)

size_factors <- sizeFactors(
  dds_qc
)

size_factors_dt <- data.table(
  sample_id = names(
    size_factors
  ),
  size_factor = as.numeric(
    size_factors
  ),
  group = metadata$group,
  sequencingBatch = metadata$sequencingBatch
)

if (
  !identical(
    size_factors_dt$sample_id,
    metadata$sample_id
  )
) {

  stop(
    "A ordem dos size factors não corresponde ao metadata."
  )
}

fwrite(
  size_factors_dt,
  size_factors_output,
  sep = "\t"
)


# ------------------------------------------------------------------------------
# VST blind = TRUE
# ------------------------------------------------------------------------------

vsd <- vst(
  dds_qc,
  blind = TRUE
)

vst_matrix <- SummarizedExperiment::assay(
  vsd
)

if (
  !identical(
    colnames(
      vst_matrix
    ),
    metadata$sample_id
  )
) {

  stop(
    "A ordem das amostras mudou durante a VST."
  )
}

saveRDS(
  vst_matrix,
  vst_output
)


# ==============================================================================
# 18. PCA — TOP GENES DE MAIOR VARIÂNCIA
# ==============================================================================

row_variance <- apply(
  vst_matrix,
  1,
  var
)

if (
  anyNA(
    row_variance
  )
) {

  stop(
    "Foram produzidas variâncias NA na matriz VST."
  )
}

n_top_actual <- min(
  PCA_NTOP,
  length(
    row_variance
  )
)

top_variable_indices <- order(
  row_variance,
  decreasing = TRUE
)[
  seq_len(
    n_top_actual
  )
]

pca_input <- t(
  vst_matrix[
    top_variable_indices,
    ,
    drop = FALSE
  ]
)

pca <- prcomp(
  pca_input,
  center = TRUE,
  scale. = FALSE
)

variance_explained <- (
  100 *
    pca$sdev^2 /
    sum(
      pca$sdev^2
    )
)

pca_scores <- as.data.table(
  pca$x,
  keep.rownames = "sample_id"
)

pca_scores <- merge(
  pca_scores,
  metadata[
    ,
    .(
      sample_id,
      individual_id,
      group,
      cogdx,
      sequencingBatch,
      msex,
      RIN,
      pmi
    )
  ],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------------------------
# Restaurar ordem do metadata após merge.
# ------------------------------------------------------------------------------

pca_scores <- pca_scores[
  match(
    metadata$sample_id,
    sample_id
  )
]

if (
  anyNA(
    pca_scores$sample_id
  )
) {

  stop(
    "Falha ao associar PCA ao metadata."
  )
}

fwrite(
  pca_scores,
  pca_scores_output,
  sep = "\t",
  na = "NA"
)


pca_variance <- data.table(
  component = paste0(
    "PC",
    seq_along(
      variance_explained
    )
  ),
  percent_variance = as.numeric(
    variance_explained
  )
)

fwrite(
  pca_variance,
  pca_variance_output,
  sep = "\t"
)


# ==============================================================================
# 19. PLOTS — LIBRARY SIZE
# ==============================================================================

library_size_plot <- ggplot(
  sample_qc,
  aes(
    x = group,
    y = library_size_before_filter,
    fill = group
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_jitter(
    width = 0.15,
    height = 0,
    alpha = 0.6,
    size = 1.8
  ) +
  labs(
    title = "ROSMAP — Library size por grupo",
    subtitle = "Raw counts harmonizados antes do filtro de baixa expressão",
    x = NULL,
    y = "Library size"
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = library_size_plot_output,
  plot = library_size_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ==============================================================================
# 20. PLOTS — GENES DETECTADOS
# ==============================================================================

detected_genes_plot <- ggplot(
  sample_qc,
  aes(
    x = group,
    y = detected_genes_gt0_before,
    fill = group
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_jitter(
    width = 0.15,
    height = 0,
    alpha = 0.6,
    size = 1.8
  ) +
  labs(
    title = "ROSMAP — Genes detectados por grupo",
    subtitle = "Número de GeneIDs com count > 0 antes da filtragem",
    x = NULL,
    y = "GeneIDs detectados"
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = detected_genes_plot_output,
  plot = detected_genes_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ==============================================================================
# 21. PLOT PCA — GROUP
# ==============================================================================

pc1_label <- paste0(
  "PC1 (",
  fmt_num(
    variance_explained[1],
    2
  ),
  "%)"
)

pc2_label <- paste0(
  "PC2 (",
  fmt_num(
    variance_explained[2],
    2
  ),
  "%)"
)

pca_group_plot <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = group
  )
) +
  geom_point(
    size = 2.7,
    alpha = 0.8
  ) +
  labs(
    title = "ROSMAP — PCA por diagnóstico",
    subtitle = paste0(
      "VST blind=TRUE; top ",
      n_top_actual,
      " genes de maior variância"
    ),
    x = pc1_label,
    y = pc2_label,
    color = "Group"
  ) +
  theme_bw(
    base_size = 13
  )

ggsave(
  filename = pca_group_plot_output,
  plot = pca_group_plot,
  width = 7.5,
  height = 6,
  dpi = 300
)


# ==============================================================================
# 22. PLOT PCA — sequencingBatch
# ==============================================================================

pca_batch_plot <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = sequencingBatch
  )
) +
  geom_point(
    size = 2.7,
    alpha = 0.8
  ) +
  labs(
    title = "ROSMAP — PCA por sequencing batch",
    subtitle = paste0(
      "VST blind=TRUE; top ",
      n_top_actual,
      " genes de maior variância"
    ),
    x = pc1_label,
    y = pc2_label,
    color = "sequencingBatch"
  ) +
  theme_bw(
    base_size = 13
  )

ggsave(
  filename = pca_batch_plot_output,
  plot = pca_batch_plot,
  width = 7.5,
  height = 6,
  dpi = 300
)


# ==============================================================================
# 23. RESUMO DO FILTRO
# ==============================================================================

gene_filter_summary <- data.table(
  metric = c(
    "genes_before_filter",
    "all_zero_genes_before_filter",
    "genes_kept_count_ge_10_in_ge_67_samples",
    "genes_removed_by_filter",
    "genes_after_filter",
    "all_zero_genes_after_filter",
    "samples_before_filter",
    "samples_after_filter",
    "all_zero_samples_before_filter",
    "all_zero_samples_after_filter",
    "total_counts_before_filter",
    "total_counts_after_filter",
    "total_counts_removed_by_filter",
    "percent_counts_retained",
    "percent_counts_removed",
    "filter_min_count",
    "filter_min_samples"
  ),
  value = c(
    n_genes_before,
    n_all_zero_genes_before,
    n_genes_after,
    n_genes_removed,
    n_genes_after,
    n_all_zero_genes_after,
    ncol(
      count_matrix
    ),
    ncol(
      filtered_matrix
    ),
    n_all_zero_samples_before,
    n_all_zero_samples_after,
    total_counts_before,
    total_counts_after,
    total_counts_removed,
    percent_counts_retained,
    percent_counts_removed,
    MIN_COUNT,
    MIN_SAMPLES
  )
)

fwrite(
  gene_filter_summary,
  gene_filter_summary_output,
  sep = "\t"
)


# ==============================================================================
# 24. RELATÓRIO TEXTUAL
# ==============================================================================

sink(
  report_output
)

cat(
  "============================================================\n"
)
cat(
  "FASE 01c — 04_qc_counts_ROSMAP.R\n"
)
cat(
  "QC, FILTRAGEM DE BAIXA EXPRESSÃO E PCA\n"
)
cat(
  "============================================================\n\n"
)


# ------------------------------------------------------------------------------
# 1. Entrada
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "1. ENTRADA\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Matriz harmonizada:\n  ",
  counts_input,
  "\n\n",
  sep = ""
)

cat(
  "Metadata:\n  ",
  metadata_input,
  "\n\n",
  sep = ""
)

cat(
  "GeneIDs antes do filtro: ",
  fmt_int(
    n_genes_before
  ),
  "\n",
  sep = ""
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
  "Ordem counts ↔ metadata: TRUE\n\n"
)


# ------------------------------------------------------------------------------
# 2. Integridade antes do filtro
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "2. INTEGRIDADE PRÉ-FILTRO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "GeneID duplicados: 0\n"
)

cat(
  "NA: ",
  fmt_int(
    n_na_before
  ),
  "\n",
  sep = ""
)

cat(
  "Não finitos: ",
  fmt_int(
    n_nonfinite_before
  ),
  "\n",
  sep = ""
)

cat(
  "Negativos: ",
  fmt_int(
    n_negative_before
  ),
  "\n",
  sep = ""
)

cat(
  "Fracionários: ",
  fmt_int(
    n_fractional_before
  ),
  "\n",
  sep = ""
)

cat(
  "Genes all-zero: ",
  fmt_int(
    n_all_zero_genes_before
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras all-zero: ",
  fmt_int(
    n_all_zero_samples_before
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 3. Filtro
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "3. FILTRAGEM DE BAIXA EXPRESSÃO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Critério pré-definido:\n"
)

cat(
  "  count >= ",
  MIN_COUNT,
  " em pelo menos ",
  MIN_SAMPLES,
  " amostras.\n\n",
  sep = ""
)

cat(
  "Justificativa:\n"
)
cat(
  "  67 corresponde ao tamanho do menor grupo experimental (Control).\n\n"
)

cat(
  "GeneIDs antes do filtro: ",
  fmt_int(
    n_genes_before
  ),
  "\n",
  sep = ""
)

cat(
  "GeneIDs mantidos: ",
  fmt_int(
    n_genes_after
  ),
  "\n",
  sep = ""
)

cat(
  "GeneIDs removidos: ",
  fmt_int(
    n_genes_removed
  ),
  "\n",
  sep = ""
)

cat(
  "Percentual de genes mantidos: ",
  fmt_num(
    100 *
      n_genes_after /
      n_genes_before,
    2
  ),
  "%\n",
  sep = ""
)

cat(
  "Percentual de genes removidos: ",
  fmt_num(
    100 *
      n_genes_removed /
      n_genes_before,
    2
  ),
  "%\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 4. Counts preservados
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "4. IMPACTO DO FILTRO SOBRE COUNTS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Counts totais antes do filtro: ",
  fmt_int(
    total_counts_before
  ),
  "\n",
  sep = ""
)

cat(
  "Counts totais após o filtro: ",
  fmt_int(
    total_counts_after
  ),
  "\n",
  sep = ""
)

cat(
  "Counts removidos: ",
  fmt_int(
    total_counts_removed
  ),
  "\n",
  sep = ""
)

cat(
  "Percentual de counts mantidos: ",
  fmt_num(
    percent_counts_retained,
    4
  ),
  "%\n",
  sep = ""
)

cat(
  "Percentual de counts removidos: ",
  fmt_num(
    percent_counts_removed,
    4
  ),
  "%\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 5. Library size
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "5. LIBRARY SIZE\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Antes do filtro:\n"
)

cat(
  "  mínimo = ",
  fmt_int(
    min(
      library_size_before
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  mediana = ",
  fmt_int(
    median(
      library_size_before
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  máximo = ",
  fmt_int(
    max(
      library_size_before
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  Wilcoxon AD vs Control p = ",
  fmt_num(
    library_wilcox_p,
    6
  ),
  "\n\n",
  sep = ""
)

cat(
  "Após o filtro:\n"
)

cat(
  "  mínimo = ",
  fmt_int(
    min(
      library_size_after
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  mediana = ",
  fmt_int(
    median(
      library_size_after
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  máximo = ",
  fmt_int(
    max(
      library_size_after
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  Wilcoxon AD vs Control p = ",
  fmt_num(
    library_after_wilcox_p,
    6
  ),
  "\n\n",
  sep = ""
)

cat(
  "Os testes são descritivos e NÃO foram utilizados para excluir amostras.\n\n"
)


# ------------------------------------------------------------------------------
# 6. Genes detectados
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "6. GENES DETECTADOS POR AMOSTRA\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Antes do filtro, count > 0:\n"
)

cat(
  "  mínimo = ",
  fmt_int(
    min(
      detected_gt0_before
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  mediana = ",
  fmt_int(
    median(
      detected_gt0_before
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  máximo = ",
  fmt_int(
    max(
      detected_gt0_before
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  Wilcoxon AD vs Control p = ",
  fmt_num(
    detected_wilcox_p,
    6
  ),
  "\n\n",
  sep = ""
)

cat(
  "Após o filtro, count > 0:\n"
)

cat(
  "  mínimo = ",
  fmt_int(
    min(
      detected_gt0_after
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  mediana = ",
  fmt_int(
    median(
      detected_gt0_after
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  máximo = ",
  fmt_int(
    max(
      detected_gt0_after
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  Wilcoxon AD vs Control p = ",
  fmt_num(
    detected_after_wilcox_p,
    6
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 7. Integridade pós-filtro
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "7. INTEGRIDADE PÓS-FILTRO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "GeneIDs finais: ",
  fmt_int(
    nrow(
      filtered_matrix
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras finais: ",
  fmt_int(
    ncol(
      filtered_matrix
    )
  ),
  "\n",
  sep = ""
)

cat(
  "GeneID duplicados: 0\n"
)

cat(
  "NA: ",
  fmt_int(
    n_na_after
  ),
  "\n",
  sep = ""
)

cat(
  "Não finitos: ",
  fmt_int(
    n_nonfinite_after
  ),
  "\n",
  sep = ""
)

cat(
  "Negativos: ",
  fmt_int(
    n_negative_after
  ),
  "\n",
  sep = ""
)

cat(
  "Fracionários: ",
  fmt_int(
    n_fractional_after
  ),
  "\n",
  sep = ""
)

cat(
  "Genes all-zero: ",
  fmt_int(
    n_all_zero_genes_after
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras all-zero: ",
  fmt_int(
    n_all_zero_samples_after
  ),
  "\n",
  sep = ""
)

cat(
  "Ordem counts ↔ metadata preservada: TRUE\n\n"
)


# ------------------------------------------------------------------------------
# 8. VST
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "8. VST PARA CONTROLE DE QUALIDADE\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "A transformação VST foi aplicada APENAS para PCA/QC.\n"
)

cat(
  "blind = TRUE\n"
)

cat(
  "design do objeto de QC = ~ 1\n"
)

cat(
  "O DESeq2 diferencial posterior utilizará raw counts filtrados, não VST.\n\n"
)

cat(
  "Size factor:\n"
)

cat(
  "  mínimo = ",
  fmt_num(
    min(
      size_factors
    ),
    4
  ),
  "\n",
  sep = ""
)

cat(
  "  mediana = ",
  fmt_num(
    median(
      size_factors
    ),
    4
  ),
  "\n",
  sep = ""
)

cat(
  "  máximo = ",
  fmt_num(
    max(
      size_factors
    ),
    4
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 9. PCA
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "9. PCA EXPLORATÓRIA\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Genes utilizados na PCA: top ",
  fmt_int(
    n_top_actual
  ),
  " por variância após VST.\n",
  sep = ""
)

cat(
  "PC1: ",
  fmt_num(
    variance_explained[1],
    2
  ),
  "%\n",
  sep = ""
)

cat(
  "PC2: ",
  fmt_num(
    variance_explained[2],
    2
  ),
  "%\n\n",
  sep = ""
)

cat(
  "Foram geradas PCA coloridas por:\n"
)
cat(
  "  - group\n"
)
cat(
  "  - sequencingBatch\n\n"
)

cat(
  "Nenhuma amostra é excluída automaticamente com base na PCA.\n"
)
cat(
  "Qualquer exclusão futura exigiria evidência técnica independente.\n\n"
)


# ------------------------------------------------------------------------------
# 10. Arquivos
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "10. ARQUIVOS GERADOS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Counts filtrados:\n  ",
  counts_filtered_output,
  "\n\n",
  sep = ""
)

cat(
  "Auditoria gene a gene do filtro:\n  ",
  gene_filter_audit_output,
  "\n\n",
  sep = ""
)

cat(
  "Resumo do filtro:\n  ",
  gene_filter_summary_output,
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
  "Testes descritivos AD vs Control:\n  ",
  group_tests_output,
  "\n\n",
  sep = ""
)

cat(
  "Size factors:\n  ",
  size_factors_output,
  "\n\n",
  sep = ""
)

cat(
  "Scores da PCA:\n  ",
  pca_scores_output,
  "\n\n",
  sep = ""
)

cat(
  "Variância da PCA:\n  ",
  pca_variance_output,
  "\n\n",
  sep = ""
)

cat(
  "Matriz VST para QC:\n  ",
  vst_output,
  "\n\n",
  sep = ""
)

cat(
  "Anotação filtrada:\n  ",
  annotation_filtered_output,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 11. Conclusão
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "11. CONCLUSÃO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "- O filtro pré-definido foi aplicado sem usar informação de significância.\n"
)

cat(
  "- O critério foi count >= 10 em pelo menos 67 amostras.\n"
)

cat(
  "- Nenhuma amostra foi removida pelo filtro de genes.\n"
)

cat(
  "- A matriz filtrada preserva raw counts inteiros e a ordem do metadata.\n"
)

cat(
  "- VST foi utilizada exclusivamente para QC/PCA.\n"
)

cat(
  "- Nenhuma exclusão de amostra foi realizada automaticamente com base em PCA.\n"
)

cat(
  "- A matriz filtrada está pronta para 05_deseq2_ROSMAP.R.\n\n"
)

cat(
  "============================================================\n"
)
cat(
  "04_qc_counts_ROSMAP.R CONCLUÍDO\n"
)
cat(
  "============================================================\n"
)

sink()


# ==============================================================================
# 25. RESUMO NO TERMINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "QC ROSMAP CONCLUÍDO\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Filtro:\n"
)

cat(
  "  count >= ",
  MIN_COUNT,
  " em >= ",
  MIN_SAMPLES,
  " amostras\n\n",
  sep = ""
)

cat(
  "GeneIDs antes: ",
  fmt_int(
    n_genes_before
  ),
  "\n",
  sep = ""
)

cat(
  "GeneIDs depois: ",
  fmt_int(
    n_genes_after
  ),
  "\n",
  sep = ""
)

cat(
  "GeneIDs removidos: ",
  fmt_int(
    n_genes_removed
  ),
  "\n\n",
  sep = ""
)

cat(
  "Counts mantidos: ",
  fmt_num(
    percent_counts_retained,
    4
  ),
  "%\n\n",
  sep = ""
)

cat(
  "PC1: ",
  fmt_num(
    variance_explained[1],
    2
  ),
  "%\n",
  sep = ""
)

cat(
  "PC2: ",
  fmt_num(
    variance_explained[2],
    2
  ),
  "%\n\n",
  sep = ""
)

cat(
  "Matriz filtrada:\n  ",
  counts_filtered_output,
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
  "Próxima etapa, após revisar QC e PCA:\n"
)
cat(
  "  05_deseq2_ROSMAP.R\n"
)
