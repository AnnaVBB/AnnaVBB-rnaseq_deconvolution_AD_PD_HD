#!/usr/bin/env Rscript

# ==============================================================================
# 03_build_reference_Mathys2019.R
# ==============================================================================
#
# FASE 02 — DECONVOLUÇÃO CELULAR
# ETAPA 03 — CONSTRUÇÃO DA REFERÊNCIA MATHYS 2019
#
# Objetivo
# --------
# Construir a referência de expressão a partir dos 24 doadores "no-pathology"
# selecionados no script 01 e dos genes harmonizados para NCBI GeneID no
# script 02.
#
# A etapa produz:
#
#   1. matriz nucleus-level harmonizada:
#        GeneID x núcleos no-pathology
#
#   2. pseudobulk raw counts:
#        GeneID x (doador x tipo celular)
#
#   3. pseudobulk CPM:
#        GeneID x (doador x tipo celular)
#
#   4. metadados alinhados de genes, núcleos e pseudobulks;
#
#   5. auditorias de conservação de counts, genes removidos, cobertura
#      donor x cell type e overlap com o bulk ROSMAP.
#
#
# DEFINIÇÃO DA REFERÊNCIA
# -----------------------
# A referência principal contém:
#
#   - 24 doadores no-pathology;
#   - 35.140 núcleos, conforme etapa 01;
#   - 8 tipos celulares:
#       Ex, In, Ast, Oli, Opc, Mic, End, Per
#
# End e Per não estão presentes em todos os doadores. Combinações ausentes
# NÃO são transformadas em pseudobulks artificiais de zeros.
#
#
# HARMONIZAÇÃO GÊNICA
# -------------------
# O script usa SOMENTE linhas classificadas como resolvidas pelo script 02:
#
#   direct_symbol
#   rescued_unique_alias
#
# Mappings ambíguos e genes não resolvidos permanecem excluídos.
#
# Se múltiplas linhas originais convergirem para o mesmo NCBI GeneID, os counts
# dessas linhas são SOMADOS. Essa operação é auditada e a conservação das
# contagens é verificada.
#
#
# PSEUDOBULK
# ----------
# Para cada GeneID g, doador d e tipo celular t:
#
#     count(g,d,t) = soma dos counts de todos os núcleos do grupo (d,t)
#
# A unidade de replicação biológica é o DOADOR, não cada núcleo.
#
# Somente combinações donor x cell type realmente observadas são criadas.
#
#
# CPM
# ---
# Para cada perfil pseudobulk:
#
#     CPM = counts / library_size * 1e6
#
# Raw counts e CPM são salvos separadamente.
#
#
# IMPORTANTE
# ----------
# Este script NÃO:
#
#   - deriva marcadores;
#   - roda deconvolução;
#   - força um universo comum entre AD, HD e PD;
#   - preenche genes ausentes com zero em datasets bulk;
#   - cria pseudobulks zero para donor x cell type ausentes;
#   - remove os 3 doadores que se sobrepõem ao bulk ROSMAP da referência
#     principal.
#
# A análise de sensibilidade sem os doadores sobrepostos será construída
# posteriormente a partir da flag já preservada no metadata.
#
#
# INPUTS PRINCIPAIS
# -----------------
#
# data/reference/mathys2019/raw/
#   filtered_count_matrix.mtx
#
# data/reference/mathys2019/phase02_prepared/
#   mathys2019_pathology_free_cell_metadata.tsv
#   mathys2019_pathology_free_donor_metadata.tsv
#   mathys2019_pathology_free_matrix_columns.rds
#
# data/reference/mathys2019/phase02_prepared/gene_harmonization/
#   mathys2019_gene_mapping.tsv
#
#
# OUTPUTS PRINCIPAIS
# ------------------
#
# data/reference/mathys2019/phase02_prepared/reference/
#   mathys2019_pathology_free_geneid_counts_sparse.rds
#   mathys2019_pathology_free_gene_metadata.tsv
#   mathys2019_pathology_free_cell_metadata_aligned.tsv
#   mathys2019_pathology_free_pseudobulk_counts_sparse.rds
#   mathys2019_pathology_free_pseudobulk_cpm_sparse.rds
#   mathys2019_pathology_free_pseudobulk_metadata.tsv
#
# ==============================================================================


# ==============================================================================
# 0. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table",
  "Matrix"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {

  stop(
    "Pacotes necessários não instalados: ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\nInstale somente os pacotes ausentes no ambiente R atual e execute ",
    "novamente. Não altere o renv automaticamente."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})


# ==============================================================================
# 1. PARÂMETROS PREDEFINIDOS
# ==============================================================================

EXPECTED_RAW_GENE_ROWS <- 17926L
EXPECTED_RAW_NUCLEI <- 70634L
EXPECTED_PATHOLOGY_FREE_DONORS <- 24L
EXPECTED_PATHOLOGY_FREE_NUCLEI <- 35140L
EXPECTED_RESOLVED_GENE_ROWS <- 17457L
EXPECTED_UNIQUE_RESOLVED_GENEIDS <- 17450L

# Derivado do QC da etapa 01:
# 24 x 8 = 192 combinações possíveis
# - 6 Per ausentes
# - 5 End ausentes
# = 181 perfis observados.
EXPECTED_PSEUDOBULK_PROFILES <- 181L

EXPECTED_CELL_TYPES <- c(
  "Ex",
  "In",
  "Ast",
  "Oli",
  "Opc",
  "Mic",
  "End",
  "Per"
)

# Ordem apenas para organização de outputs.
CELL_TYPE_ORDER <- c(
  "Ex",
  "In",
  "Ast",
  "Oli",
  "Opc",
  "Mic",
  "End",
  "Per"
)

COUNT_TOLERANCE <- 1e-8
CPM_SUM_TOLERANCE <- 1e-4


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

raw_dir <- file.path(
  "data",
  "reference",
  "mathys2019",
  "raw"
)

prepared_dir <- file.path(
  "data",
  "reference",
  "mathys2019",
  "phase02_prepared"
)

gene_harmonization_dir <- file.path(
  prepared_dir,
  "gene_harmonization"
)

reference_dir <- file.path(
  prepared_dir,
  "reference"
)

result_dir <- file.path(
  "results",
  "deconvolution",
  "mathys2019",
  "reference_build"
)

dir.create(
  reference_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------

matrix_file <- file.path(
  raw_dir,
  "filtered_count_matrix.mtx"
)

cell_metadata_file <- file.path(
  prepared_dir,
  "mathys2019_pathology_free_cell_metadata.tsv"
)

donor_metadata_file <- file.path(
  prepared_dir,
  "mathys2019_pathology_free_donor_metadata.tsv"
)

selected_columns_file <- file.path(
  prepared_dir,
  "mathys2019_pathology_free_matrix_columns.rds"
)

gene_mapping_file <- file.path(
  gene_harmonization_dir,
  "mathys2019_gene_mapping.tsv"
)

step01_status_file <- file.path(
  "results",
  "deconvolution",
  "mathys2019",
  "metadata",
  "01_prepare_metadata_status.tsv"
)

step02_status_file <- file.path(
  "results",
  "deconvolution",
  "mathys2019",
  "gene_harmonization",
  "02_harmonization_status.tsv"
)


# ------------------------------------------------------------------------------
# Outputs nucleus-level
# ------------------------------------------------------------------------------

nucleus_counts_file <- file.path(
  reference_dir,
  "mathys2019_pathology_free_geneid_counts_sparse.rds"
)

gene_metadata_output_file <- file.path(
  reference_dir,
  "mathys2019_pathology_free_gene_metadata.tsv"
)

cell_metadata_output_file <- file.path(
  reference_dir,
  "mathys2019_pathology_free_cell_metadata_aligned.tsv"
)


# ------------------------------------------------------------------------------
# Outputs pseudobulk
# ------------------------------------------------------------------------------

pseudobulk_counts_file <- file.path(
  reference_dir,
  "mathys2019_pathology_free_pseudobulk_counts_sparse.rds"
)

pseudobulk_cpm_file <- file.path(
  reference_dir,
  "mathys2019_pathology_free_pseudobulk_cpm_sparse.rds"
)

pseudobulk_metadata_file <- file.path(
  reference_dir,
  "mathys2019_pathology_free_pseudobulk_metadata.tsv"
)


# ------------------------------------------------------------------------------
# Logs / QC
# ------------------------------------------------------------------------------

log_file <- file.path(
  result_dir,
  "03_build_reference_Mathys2019.txt"
)

session_file <- file.path(
  result_dir,
  "03_build_reference_Mathys2019_sessionInfo.txt"
)


# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================

clean_chr <- function(x) {

  x <- trimws(
    as.character(x)
  )

  x[
    is.na(x) |
      x == "" |
      toupper(x) == "NA"
  ] <- NA_character_

  x
}


parse_logical <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  y <- toupper(
    clean_chr(x)
  )

  out <- rep(
    NA,
    length(y)
  )

  out[
    y %in%
      c(
        "TRUE",
        "T",
        "1",
        "YES",
        "Y"
      )
  ] <- TRUE

  out[
    y %in%
      c(
        "FALSE",
        "F",
        "0",
        "NO",
        "N"
      )
  ] <- FALSE

  out
}


collapse_values <- function(x) {

  x <- clean_chr(x)

  x <- sort(
    unique(
      x[
        !is.na(x)
      ]
    )
  )

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste(
    x,
    collapse = ";"
  )
}


make_check <- function(
  check,
  observed,
  expected,
  status,
  note = NA_character_
) {

  data.table(
    check = check,
    observed = as.character(observed),
    expected = as.character(expected),
    status = status,
    note = note
  )
}


validate_status_file <- function(
  file,
  label
) {

  if (!file.exists(file)) {

    stop(
      label,
      ": arquivo de status não encontrado:\n",
      file
    )
  }

  dt <- fread(
    file,
    na.strings = c(
      "NA",
      ""
    )
  )

  if (!"status" %in% names(dt)) {

    stop(
      label,
      ": coluna 'status' não encontrada."
    )
  }

  if (
    any(
      dt$status == "FAIL",
      na.rm = TRUE
    )
  ) {

    stop(
      label,
      ": existe pelo menos um check FAIL."
    )
  }

  invisible(TRUE)
}


safe_sparse_sum <- function(x) {

  if (!inherits(x, "sparseMatrix")) {

    stop(
      "safe_sparse_sum recebeu objeto que não é sparseMatrix."
    )
  }

  sum(
    x@x
  )
}


check_sparse_counts <- function(
  x,
  label
) {

  if (!inherits(x, "sparseMatrix")) {

    stop(
      label,
      ": objeto não é sparseMatrix."
    )
  }

  if (
    any(
      !is.finite(
        x@x
      )
    )
  ) {

    stop(
      label,
      ": existem valores não finitos."
    )
  }

  if (
    any(
      x@x < 0
    )
  ) {

    stop(
      label,
      ": existem counts negativos."
    )
  }

  if (
    any(
      abs(
        x@x -
          round(
            x@x
          )
      ) >
        COUNT_TOLERANCE
    )
  ) {

    stop(
      label,
      ": existem valores fracionários inesperados."
    )
  }

  invisible(TRUE)
}


# ==============================================================================
# 4. FUNÇÃO PRINCIPAL
# ==============================================================================

run_build_reference <- function() {


  # ============================================================================
  # 4.1 CABEÇALHO E PRÉ-REQUISITOS
  # ============================================================================

  cat(
    "============================================================\n"
  )

  cat(
    "FASE 02 — 03_build_reference_Mathys2019.R\n"
  )

  cat(
    "CONSTRUÇÃO DA REFERÊNCIA PATHOLOGY-FREE\n"
  )

  cat(
    "============================================================\n\n"
  )

  validate_status_file(
    step01_status_file,
    "Etapa 01"
  )

  validate_status_file(
    step02_status_file,
    "Etapa 02"
  )

  cat(
    "Etapas 01 e 02 validadas: nenhum check FAIL.\n"
  )

  required_files <- c(
    raw_matrix = matrix_file,
    cell_metadata = cell_metadata_file,
    donor_metadata = donor_metadata_file,
    selected_columns = selected_columns_file,
    gene_mapping = gene_mapping_file
  )

  missing_files <- required_files[
    !file.exists(
      required_files
    )
  ]

  if (length(missing_files) > 0L) {

    cat(
      "\nArquivos obrigatórios ausentes:\n"
    )

    print(
      missing_files
    )

    stop(
      "Construção da referência interrompida."
    )
  }

  input_audit <- data.table(
    role = names(
      required_files
    ),
    file_name = basename(
      required_files
    ),
    size_MB = round(
      file.info(
        required_files
      )$size /
        1024^2,
      3
    )
  )

  fwrite(
    input_audit,
    file.path(
      result_dir,
      "03_input_file_audit.tsv"
    ),
    sep = "\t"
  )

  cat(
    "\nInputs encontrados:\n"
  )

  print(
    input_audit
  )


  # ============================================================================
  # 4.2 CARREGAR METADATA E ÍNDICES
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "1. METADATA E ÍNDICES PATHOLOGY-FREE\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  cell_meta <- fread(
    cell_metadata_file,
    na.strings = c(
      "NA",
      ""
    )
  )

  donor_meta <- fread(
    donor_metadata_file,
    na.strings = c(
      "NA",
      ""
    )
  )

  selected_columns <- readRDS(
    selected_columns_file
  )

  required_cell_columns <- c(
    "matrix_col_index",
    "TAG",
    "donor_id",
    "cell_type",
    "overlap_with_bulk_ROSMAP"
  )

  missing_cell_columns <- setdiff(
    required_cell_columns,
    names(
      cell_meta
    )
  )

  if (length(missing_cell_columns) > 0L) {

    stop(
      "Metadata de núcleos sem colunas obrigatórias: ",
      paste(
        missing_cell_columns,
        collapse = ", "
      )
    )
  }

  cell_meta[
    ,
    matrix_col_index :=
      as.integer(
        matrix_col_index
      )
  ]

  cell_meta[
    ,
    TAG :=
      clean_chr(
        TAG
      )
  ]

  cell_meta[
    ,
    donor_id :=
      clean_chr(
        donor_id
      )
  ]

  cell_meta[
    ,
    cell_type :=
      clean_chr(
        cell_type
      )
  ]

  cell_meta[
    ,
    overlap_with_bulk_ROSMAP :=
      parse_logical(
        overlap_with_bulk_ROSMAP
      )
  ]

  if (
    anyNA(
      cell_meta$matrix_col_index
    )
  ) {

    stop(
      "Existem matrix_col_index ausentes."
    )
  }

  if (
    anyDuplicated(
      cell_meta$matrix_col_index
    )
  ) {

    stop(
      "Existem matrix_col_index duplicados."
    )
  }

  if (
    anyNA(
      cell_meta$TAG
    ) |
      anyDuplicated(
        cell_meta$TAG
      )
  ) {

    stop(
      "TAG deve ser completo e único no metadata pathology-free."
    )
  }

  selected_columns <- as.integer(
    selected_columns
  )

  if (
    length(
      selected_columns
    ) !=
      EXPECTED_PATHOLOGY_FREE_NUCLEI
  ) {

    stop(
      "Esperávamos ",
      EXPECTED_PATHOLOGY_FREE_NUCLEI,
      " índices de núcleos; encontrados: ",
      length(
        selected_columns
      ),
      "."
    )
  }

  if (
    anyDuplicated(
      selected_columns
    )
  ) {

    stop(
      "O vetor de colunas selecionadas contém duplicatas."
    )
  }

  # Alinhar metadata EXATAMENTE à ordem das colunas que será usada na matriz.
  idx_meta <- match(
    selected_columns,
    cell_meta$matrix_col_index
  )

  if (
    anyNA(
      idx_meta
    )
  ) {

    stop(
      "Há índices da matriz sem correspondência no metadata pathology-free."
    )
  }

  cell_meta <- cell_meta[
    idx_meta
  ]

  if (
    !identical(
      cell_meta$matrix_col_index,
      selected_columns
    )
  ) {

    stop(
      "Falha ao alinhar metadata às colunas selecionadas."
    )
  }

  n_selected_donors <- uniqueN(
    cell_meta$donor_id
  )

  n_selected_nuclei <- nrow(
    cell_meta
  )

  observed_cell_types <- sort(
    unique(
      cell_meta$cell_type
    )
  )

  cat(
    "Núcleos selecionados: ",
    n_selected_nuclei,
    "\n",
    sep = ""
  )

  cat(
    "Doadores selecionados: ",
    n_selected_donors,
    "\n",
    sep = ""
  )

  cat(
    "Tipos celulares: ",
    paste(
      observed_cell_types,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )

  if (
    n_selected_nuclei !=
      EXPECTED_PATHOLOGY_FREE_NUCLEI
  ) {

    stop(
      "Número inesperado de núcleos pathology-free."
    )
  }

  if (
    n_selected_donors !=
      EXPECTED_PATHOLOGY_FREE_DONORS
  ) {

    stop(
      "Número inesperado de doadores pathology-free."
    )
  }

  if (
    !setequal(
      observed_cell_types,
      EXPECTED_CELL_TYPES
    )
  ) {

    stop(
      "Os 8 tipos celulares esperados não estão todos representados."
    )
  }


  # ============================================================================
  # 4.3 CARREGAR MAPA GÊNICO
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "2. MAPA GÊNICO HARMONIZADO\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  gene_map <- fread(
    gene_mapping_file,
    na.strings = c(
      "NA",
      ""
    )
  )

  required_gene_columns <- c(
    "matrix_row_index",
    "original_symbol",
    "gene_id",
    "mapping_method",
    "resolved",
    "geneid_n_original_rows",
    "geneid_collision"
  )

  missing_gene_columns <- setdiff(
    required_gene_columns,
    names(
      gene_map
    )
  )

  if (
    length(
      missing_gene_columns
    ) > 0L
  ) {

    stop(
      "Mapa gênico sem colunas obrigatórias: ",
      paste(
        missing_gene_columns,
        collapse = ", "
      )
    )
  }

  gene_map[
    ,
    matrix_row_index :=
      as.integer(
        matrix_row_index
      )
  ]

  gene_map[
    ,
    gene_id :=
      clean_chr(
        gene_id
      )
  ]

  gene_map[
    ,
    original_symbol :=
      clean_chr(
        original_symbol
      )
  ]

  gene_map[
    ,
    resolved :=
      parse_logical(
        resolved
      )
  ]

  gene_map[
    ,
    geneid_collision :=
      parse_logical(
        geneid_collision
      )
  ]

  if (
    nrow(
      gene_map
    ) !=
      EXPECTED_RAW_GENE_ROWS
  ) {

    stop(
      "Esperávamos ",
      EXPECTED_RAW_GENE_ROWS,
      " linhas no mapa gênico; encontradas: ",
      nrow(
        gene_map
      ),
      "."
    )
  }

  if (
    !identical(
      gene_map$matrix_row_index,
      seq_len(
        nrow(
          gene_map
        )
      )
    )
  ) {

    stop(
      "matrix_row_index do mapa gênico não corresponde às linhas originais."
    )
  }

  if (
    anyNA(
      gene_map$resolved
    )
  ) {

    stop(
      "A coluna resolved contém NA."
    )
  }

  resolved_map <- gene_map[
    resolved == TRUE
  ]

  unresolved_map <- gene_map[
    resolved == FALSE
  ]

  if (
    nrow(
      resolved_map
    ) !=
      EXPECTED_RESOLVED_GENE_ROWS
  ) {

    stop(
      "Esperávamos ",
      EXPECTED_RESOLVED_GENE_ROWS,
      " linhas gênicas resolvidas; encontradas: ",
      nrow(
        resolved_map
      ),
      "."
    )
  }

  if (
    uniqueN(
      resolved_map$gene_id
    ) !=
      EXPECTED_UNIQUE_RESOLVED_GENEIDS
  ) {

    stop(
      "Esperávamos ",
      EXPECTED_UNIQUE_RESOLVED_GENEIDS,
      " GeneIDs únicos resolvidos; encontrados: ",
      uniqueN(
        resolved_map$gene_id
      ),
      "."
    )
  }

  if (
    any(
      is.na(
        resolved_map$gene_id
      )
    )
  ) {

    stop(
      "Há linhas resolvidas sem GeneID."
    )
  }

  if (
    any(
      !grepl(
        "^[0-9]+$",
        resolved_map$gene_id
      )
    )
  ) {

    stop(
      "Há GeneID não numérico entre as linhas resolvidas."
    )
  }

  cat(
    "Linhas gênicas resolvidas: ",
    nrow(
      resolved_map
    ),
    "\n",
    sep = ""
  )

  cat(
    "GeneIDs únicos antes do filtro de expressão: ",
    uniqueN(
      resolved_map$gene_id
    ),
    "\n",
    sep = ""
  )

  cat(
    "Linhas não resolvidas/ambíguas: ",
    nrow(
      unresolved_map
    ),
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.4 CARREGAR MATRIZ RAW
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "3. CARREGAMENTO DA MATRIZ RAW\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  cat(
    "Lendo Matrix Market. Esta é a etapa mais pesada do script.\n"
  )

  raw_counts <- Matrix::readMM(
    matrix_file
  )

  if (
    nrow(
      raw_counts
    ) !=
      EXPECTED_RAW_GENE_ROWS |
      ncol(
        raw_counts
      ) !=
      EXPECTED_RAW_NUCLEI
  ) {

    stop(
      "Dimensão inesperada da matriz raw: ",
      nrow(
        raw_counts
      ),
      " x ",
      ncol(
        raw_counts
      ),
      ". Esperado: ",
      EXPECTED_RAW_GENE_ROWS,
      " x ",
      EXPECTED_RAW_NUCLEI,
      "."
    )
  }

  # Converter uma única vez para formato eficiente para subset por coluna.
  raw_counts <- as(
    raw_counts,
    "dgCMatrix"
  )

  check_sparse_counts(
    raw_counts,
    "Matriz raw"
  )

  cat(
    "Matriz raw: ",
    nrow(
      raw_counts
    ),
    " genes x ",
    ncol(
      raw_counts
    ),
    " núcleos\n",
    sep = ""
  )

  cat(
    "Non-zero entries: ",
    length(
      raw_counts@x
    ),
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.5 SUBSET DOS 35.140 NÚCLEOS PATHOLOGY-FREE
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "4. SUBSET DOS NÚCLEOS PATHOLOGY-FREE\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  if (
    any(
      selected_columns < 1L |
        selected_columns >
          ncol(
            raw_counts
          )
    )
  ) {

    stop(
      "Existem índices de colunas fora dos limites da matriz."
    )
  }

  selected_all_genes <- raw_counts[
    ,
    selected_columns,
    drop = FALSE
  ]

  # Liberar a matriz completa tão cedo quanto possível.
  rm(
    raw_counts
  )

  invisible(
    gc()
  )

  if (
    ncol(
      selected_all_genes
    ) !=
      EXPECTED_PATHOLOGY_FREE_NUCLEI
  ) {

    stop(
      "Falha no subset dos núcleos pathology-free."
    )
  }

  check_sparse_counts(
    selected_all_genes,
    "Matriz pathology-free antes da harmonização gênica"
  )

  total_counts_selected_all_genes <- safe_sparse_sum(
    selected_all_genes
  )

  row_totals_original <- Matrix::rowSums(
    selected_all_genes
  )

  cat(
    "Dimensão após selecionar núcleos: ",
    nrow(
      selected_all_genes
    ),
    " x ",
    ncol(
      selected_all_genes
    ),
    "\n",
    sep = ""
  )

  cat(
    "Counts totais nos 35.140 núcleos: ",
    format(
      total_counts_selected_all_genes,
      scientific = FALSE,
      trim = TRUE
    ),
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.6 QUANTIFICAR COUNTS PERDIDOS PELA HARMONIZAÇÃO
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "5. IMPACTO DA HARMONIZAÇÃO SOBRE COUNTS\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  resolved_rows <- resolved_map$matrix_row_index

  unresolved_rows <- unresolved_map$matrix_row_index

  resolved_counts_before_collapse <- sum(
    row_totals_original[
      resolved_rows
    ]
  )

  unresolved_counts_total <- sum(
    row_totals_original[
      unresolved_rows
    ]
  )

  unresolved_counts_percent <- if (
    total_counts_selected_all_genes > 0
  ) {

    100 *
      unresolved_counts_total /
      total_counts_selected_all_genes

  } else {

    NA_real_
  }

  cat(
    "Counts em genes resolvidos: ",
    format(
      resolved_counts_before_collapse,
      scientific = FALSE,
      trim = TRUE
    ),
    "\n",
    sep = ""
  )

  cat(
    "Counts em genes não resolvidos/ambíguos: ",
    format(
      unresolved_counts_total,
      scientific = FALSE,
      trim = TRUE
    ),
    "\n",
    sep = ""
  )

  cat(
    "Percentual de counts removidos pela harmonização: ",
    round(
      unresolved_counts_percent,
      6
    ),
    "%\n",
    sep = ""
  )

  # Auditoria por categoria de mapping.
  mapping_count_loss <- copy(
    gene_map[
      ,
      .(
        matrix_row_index,
        mapping_method,
        resolved
      )
    ]
  )

  mapping_count_loss[
    ,
    total_counts_selected_nuclei :=
      as.numeric(
        row_totals_original[
          matrix_row_index
        ]
      )
  ]

  mapping_count_loss_summary <- mapping_count_loss[
    ,
    .(
      n_gene_rows = .N,
      total_counts = sum(
        total_counts_selected_nuclei
      )
    ),
    by = .(
      mapping_method,
      resolved
    )
  ][
    order(
      -total_counts
    )
  ]

  mapping_count_loss_summary[
    ,
    percent_of_all_counts :=
      100 *
      total_counts /
      total_counts_selected_all_genes
  ]

  fwrite(
    mapping_count_loss_summary,
    file.path(
      result_dir,
      "03_counts_by_mapping_method.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(
    mapping_count_loss_summary
  )


  # ============================================================================
  # 4.7 SUBSET DAS LINHAS RESOLVIDAS
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "6. SUBSET DOS GENES RESOLVIDOS\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  resolved_counts <- selected_all_genes[
    resolved_rows,
    ,
    drop = FALSE
  ]

  # selected_all_genes não é mais necessário depois deste ponto.
  rm(
    selected_all_genes
  )

  invisible(
    gc()
  )

  check_sparse_counts(
    resolved_counts,
    "Matriz de genes resolvidos"
  )

  if (
    nrow(
      resolved_counts
    ) !=
      nrow(
        resolved_map
      )
  ) {

    stop(
      "Número de linhas da matriz resolvida não corresponde ao mapa resolvido."
    )
  }

  total_resolved_matrix <- safe_sparse_sum(
    resolved_counts
  )

  if (
    abs(
      total_resolved_matrix -
        resolved_counts_before_collapse
    ) >
      COUNT_TOLERANCE
  ) {

    stop(
      "Counts dos genes resolvidos mudaram durante o subset."
    )
  }

  cat(
    "Matriz resolvida antes do colapso: ",
    nrow(
      resolved_counts
    ),
    " x ",
    ncol(
      resolved_counts
    ),
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.8 COLAPSAR LINHAS QUE CONVERGEM PARA O MESMO GeneID
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "7. COLAPSO PARA GeneID ÚNICO\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  # Ordenação numérica do GeneID para que o output tenha ordem determinística.
  unique_gene_ids <- sort(
    unique(
      as.numeric(
        resolved_map$gene_id
      )
    )
  )

  unique_gene_ids <- as.character(
    unique_gene_ids
  )

  if (
    length(
      unique_gene_ids
    ) !=
      EXPECTED_UNIQUE_RESOLVED_GENEIDS
  ) {

    stop(
      "Número inesperado de GeneIDs antes do colapso."
    )
  }

  gene_group_index <- match(
    resolved_map$gene_id,
    unique_gene_ids
  )

  if (
    anyNA(
      gene_group_index
    )
  ) {

    stop(
      "Falha ao construir índice de colapso por GeneID."
    )
  }

  # Matriz de agregação extremamente esparsa:
  #
  #   rows = GeneIDs únicos
  #   cols = linhas gênicas resolvidas
  #
  # Cada coluna possui exatamente um valor 1.
  gene_aggregation_matrix <- Matrix::sparseMatrix(
    i = gene_group_index,
    j = seq_len(
      nrow(
        resolved_map
      )
    ),
    x = rep(
      1,
      nrow(
        resolved_map
      )
    ),
    dims = c(
      length(
        unique_gene_ids
      ),
      nrow(
        resolved_map
      )
    )
  )

  geneid_counts <- gene_aggregation_matrix %*%
    resolved_counts

  geneid_counts <- as(
    geneid_counts,
    "dgCMatrix"
  )

  rm(
    gene_aggregation_matrix,
    resolved_counts
  )

  invisible(
    gc()
  )

  check_sparse_counts(
    geneid_counts,
    "Matriz colapsada por GeneID"
  )

  total_counts_after_geneid_collapse <- safe_sparse_sum(
    geneid_counts
  )

  if (
    abs(
      total_counts_after_geneid_collapse -
        total_resolved_matrix
    ) >
      COUNT_TOLERANCE
  ) {

    stop(
      "Counts não foram conservados durante o colapso para GeneID."
    )
  }

  rownames(
    geneid_counts
  ) <- unique_gene_ids

  colnames(
    geneid_counts
  ) <- cell_meta$TAG

  cat(
    "Matriz após colapso: ",
    nrow(
      geneid_counts
    ),
    " GeneIDs x ",
    ncol(
      geneid_counts
    ),
    " núcleos\n",
    sep = ""
  )

  cat(
    "Counts conservados no colapso: SIM\n"
  )


  # ============================================================================
  # 4.9 REMOVER GeneIDs SEM EXPRESSÃO NOS 24 NO-PATHOLOGY
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "8. GeneIDs ALL-ZERO NOS 24 NO-PATHOLOGY\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  geneid_totals <- Matrix::rowSums(
    geneid_counts
  )

  zero_gene_idx <- which(
    geneid_totals == 0
  )

  nonzero_gene_idx <- which(
    geneid_totals > 0
  )

  n_zero_geneids <- length(
    zero_gene_idx
  )

  cat(
    "GeneIDs com count total zero: ",
    n_zero_geneids,
    "\n",
    sep = ""
  )

  if (
    n_zero_geneids > 0L
  ) {

    zero_gene_table <- data.table(
      gene_id = rownames(
        geneid_counts
      )[
        zero_gene_idx
      ]
    )

    fwrite(
      zero_gene_table,
      file.path(
        result_dir,
        "03_all_zero_geneids_removed.tsv"
      ),
      sep = "\t"
    )

    geneid_counts <- geneid_counts[
      nonzero_gene_idx,
      ,
      drop = FALSE
    ]

    geneid_totals <- geneid_totals[
      nonzero_gene_idx
    ]

  } else {

    fwrite(
      data.table(
        note = "Nenhum GeneID all-zero após seleção dos 24 no-pathology."
      ),
      file.path(
        result_dir,
        "03_all_zero_geneids_removed.tsv"
      ),
      sep = "\t"
    )
  }

  final_gene_ids <- rownames(
    geneid_counts
  )

  total_counts_after_zero_filter <- safe_sparse_sum(
    geneid_counts
  )

  if (
    abs(
      total_counts_after_zero_filter -
        total_counts_after_geneid_collapse
    ) >
      COUNT_TOLERANCE
  ) {

    stop(
      "A remoção de genes all-zero alterou os counts totais."
    )
  }

  cat(
    "GeneIDs finais com expressão > 0: ",
    nrow(
      geneid_counts
    ),
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.10 CONSTRUIR METADATA DE GENES
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "9. METADATA DOS GeneIDs FINAIS\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  gene_metadata <- resolved_map[
    gene_id %in%
      final_gene_ids,
    .(
      n_original_rows = .N,
      n_original_symbols = uniqueN(
        original_symbol
      ),
      original_symbols = collapse_values(
        original_symbol
      ),
      mapping_methods = collapse_values(
        mapping_method
      ),
      annotation_symbols = if (
        "annotation_symbol" %in%
          names(
            resolved_map
          )
      ) {
        collapse_values(
          annotation_symbol
        )
      } else {
        NA_character_
      },
      annotation_ensembl = if (
        "annotation_ensembl" %in%
          names(
            resolved_map
          )
      ) {
        collapse_values(
          annotation_ensembl
        )
      } else {
        NA_character_
      },
      annotation_gene_names = if (
        "annotation_gene_name" %in%
          names(
            resolved_map
          )
      ) {
        collapse_values(
          annotation_gene_name
        )
      } else {
        NA_character_
      }
    ),
    by = gene_id
  ]

  gene_metadata[
    ,
    total_counts_pathology_free_nuclei :=
      as.numeric(
        geneid_totals[
          match(
            gene_id,
            final_gene_ids
          )
        ]
      )
  ]

  gene_metadata[
    ,
    collapsed_multiple_input_rows :=
      n_original_rows > 1L
  ]

  # Mesma ordem da matriz.
  gene_metadata <- gene_metadata[
    match(
      final_gene_ids,
      gene_id
    )
  ]

  if (
    anyNA(
      gene_metadata$gene_id
    )
  ) {

    stop(
      "Falha ao alinhar gene metadata à matriz."
    )
  }

  if (
    !identical(
      gene_metadata$gene_id,
      final_gene_ids
    )
  ) {

    stop(
      "gene_metadata não está na mesma ordem da matriz."
    )
  }

  fwrite(
    gene_metadata,
    gene_metadata_output_file,
    sep = "\t",
    na = "NA"
  )

  cat(
    "Gene metadata: ",
    nrow(
      gene_metadata
    ),
    " GeneIDs\n",
    sep = ""
  )


  # ============================================================================
  # 4.11 SALVAR MATRIZ NUCLEUS-LEVEL HARMONIZADA
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "10. SALVAR MATRIZ NUCLEUS-LEVEL\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  # Salvar metadata alinhado exatamente à ordem das colunas.
  fwrite(
    cell_meta,
    cell_metadata_output_file,
    sep = "\t",
    na = "NA"
  )

  # compress = FALSE:
  # evita custo de compressão e deixa o carregamento futuro mais rápido.
  saveRDS(
    geneid_counts,
    nucleus_counts_file,
    compress = FALSE
  )

  cat(
    "Matriz nucleus-level salva:\n  ",
    nucleus_counts_file,
    "\n",
    sep = ""
  )

  cat(
    "Dimensão: ",
    nrow(
      geneid_counts
    ),
    " x ",
    ncol(
      geneid_counts
    ),
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.12 CONSTRUIR GRUPOS PSEUDOBULK
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "11. DEFINIÇÃO DOS GRUPOS PSEUDOBULK\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  # Cada núcleo pertence a exatamente um grupo donor x cell type.
  cell_meta[
    ,
    pseudobulk_id :=
      paste(
        donor_id,
        cell_type,
        sep = "__"
      )
  ]

  if (
    anyNA(
      cell_meta$pseudobulk_id
    )
  ) {

    stop(
      "Existem núcleos sem pseudobulk_id."
    )
  }

  # Metadata observado: NÃO criar combinações ausentes.
  pseudobulk_meta <- cell_meta[
    ,
    .(
      n_nuclei = .N,
      overlap_with_bulk_ROSMAP = {
        x <- unique(
          overlap_with_bulk_ROSMAP
        )

        x <- x[
          !is.na(
            x
          )
        ]

        if (length(x) == 1L) {
          x
        } else if (length(x) == 0L) {
          NA
        } else {
          stop(
            "Flag de overlap inconsistente dentro de um donor x cell type."
          )
        }
      }
    ),
    by = .(
      pseudobulk_id,
      donor_id,
      cell_type
    )
  ]

  pseudobulk_meta[
    ,
    cell_type_order :=
      match(
        cell_type,
        CELL_TYPE_ORDER
      )
  ]

  if (
    anyNA(
      pseudobulk_meta$cell_type_order
    )
  ) {

    stop(
      "Tipo celular inesperado no pseudobulk metadata."
    )
  }

  setorder(
    pseudobulk_meta,
    donor_id,
    cell_type_order
  )

  pseudobulk_meta[
    ,
    cell_type_order := NULL
  ]

  n_pseudobulks <- nrow(
    pseudobulk_meta
  )

  cat(
    "Perfis pseudobulk observados: ",
    n_pseudobulks,
    "\n",
    sep = ""
  )

  if (
    n_pseudobulks !=
      EXPECTED_PSEUDOBULK_PROFILES
  ) {

    stop(
      "Esperávamos ",
      EXPECTED_PSEUDOBULK_PROFILES,
      " perfis pseudobulk com base no QC da etapa 01; encontrados: ",
      n_pseudobulks,
      "."
    )
  }

  group_index <- match(
    cell_meta$pseudobulk_id,
    pseudobulk_meta$pseudobulk_id
  )

  if (
    anyNA(
      group_index
    )
  ) {

    stop(
      "Há núcleos sem grupo pseudobulk."
    )
  }

  if (
    length(
      group_index
    ) !=
      ncol(
        geneid_counts
      )
  ) {

    stop(
      "Número de grupos por núcleo não corresponde às colunas da matriz."
    )
  }

  # Auditoria explícita das combinações ausentes.
  donor_cell_complete <- merge(
    CJ(
      donor_id = sort(
        unique(
          cell_meta$donor_id
        )
      ),
      cell_type = EXPECTED_CELL_TYPES,
      unique = TRUE
    ),
    pseudobulk_meta[
      ,
      .(
        donor_id,
        cell_type,
        pseudobulk_id,
        n_nuclei
      )
    ],
    by = c(
      "donor_id",
      "cell_type"
    ),
    all.x = TRUE,
    sort = TRUE
  )

  donor_cell_complete[
    ,
    observed :=
      !is.na(
        pseudobulk_id
      )
  ]

  fwrite(
    donor_cell_complete,
    file.path(
      result_dir,
      "03_donor_celltype_pseudobulk_coverage.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  absent_profiles <- donor_cell_complete[
    observed == FALSE
  ]

  cat(
    "Combinações donor x cell type ausentes: ",
    nrow(
      absent_profiles
    ),
    "\n",
    sep = ""
  )

  if (
    nrow(
      absent_profiles
    ) > 0L
  ) {

    cat(
      "As combinações ausentes NÃO serão representadas por vetores de zero.\n"
    )
  }


  # ============================================================================
  # 4.13 AGREGAR NÚCLEOS EM PSEUDOBULKS
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "12. AGREGAÇÃO PSEUDOBULK POR SOMA\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  # Matriz núcleos x pseudobulks.
  # Cada linha possui exatamente um 1.
  cell_to_pseudobulk <- Matrix::sparseMatrix(
    i = seq_len(
      ncol(
        geneid_counts
      )
    ),
    j = group_index,
    x = rep(
      1,
      ncol(
        geneid_counts
      )
    ),
    dims = c(
      ncol(
        geneid_counts
      ),
      n_pseudobulks
    )
  )

  pseudobulk_counts <- geneid_counts %*%
    cell_to_pseudobulk

  pseudobulk_counts <- as(
    pseudobulk_counts,
    "dgCMatrix"
  )

  rm(
    cell_to_pseudobulk
  )

  invisible(
    gc()
  )

  rownames(
    pseudobulk_counts
  ) <- rownames(
    geneid_counts
  )

  colnames(
    pseudobulk_counts
  ) <- pseudobulk_meta$pseudobulk_id

  check_sparse_counts(
    pseudobulk_counts,
    "Matriz pseudobulk raw counts"
  )

  total_nucleus_counts_final <- safe_sparse_sum(
    geneid_counts
  )

  total_pseudobulk_counts <- safe_sparse_sum(
    pseudobulk_counts
  )

  if (
    abs(
      total_nucleus_counts_final -
        total_pseudobulk_counts
    ) >
      COUNT_TOLERANCE
  ) {

    stop(
      "Counts não foram conservados na agregação nucleus -> pseudobulk."
    )
  }

  cat(
    "Dimensão pseudobulk raw counts: ",
    nrow(
      pseudobulk_counts
    ),
    " genes x ",
    ncol(
      pseudobulk_counts
    ),
    " perfis\n",
    sep = ""
  )

  cat(
    "Counts conservados nucleus -> pseudobulk: SIM\n"
  )


  # ============================================================================
  # 4.14 LIBRARY SIZE E CPM
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "13. LIBRARY SIZE E CPM\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  library_sizes <- Matrix::colSums(
    pseudobulk_counts
  )

  if (
    any(
      !is.finite(
        library_sizes
      )
    )
  ) {

    stop(
      "Library size não finito."
    )
  }

  if (
    any(
      library_sizes <= 0
    )
  ) {

    stop(
      "Existe pseudobulk com library size <= 0."
    )
  }

  pseudobulk_meta[
    ,
    library_size :=
      as.numeric(
        library_sizes[
          match(
            pseudobulk_id,
            names(
              library_sizes
            )
          )
        ]
      )
  ]

  if (
    anyNA(
      pseudobulk_meta$library_size
    )
  ) {

    stop(
      "Falha ao associar library size ao metadata pseudobulk."
    )
  }

  # Escalonamento de cada coluna por 1e6 / library_size.
  cpm_scaling <- Matrix::Diagonal(
    x = 1e6 /
      library_sizes
  )

  pseudobulk_cpm <- pseudobulk_counts %*%
    cpm_scaling

  pseudobulk_cpm <- as(
    pseudobulk_cpm,
    "dgCMatrix"
  )

  rm(
    cpm_scaling
  )

  invisible(
    gc()
  )

  rownames(
    pseudobulk_cpm
  ) <- rownames(
    pseudobulk_counts
  )

  colnames(
    pseudobulk_cpm
  ) <- colnames(
    pseudobulk_counts
  )

  cpm_column_sums <- Matrix::colSums(
    pseudobulk_cpm
  )

  max_cpm_sum_deviation <- max(
    abs(
      cpm_column_sums -
        1e6
    )
  )

  if (
    max_cpm_sum_deviation >
      CPM_SUM_TOLERANCE
  ) {

    stop(
      "As colunas CPM não somam aproximadamente 1e6. ",
      "Desvio máximo: ",
      max_cpm_sum_deviation
    )
  }

  cat(
    "Library size mínima: ",
    format(
      min(
        library_sizes
      ),
      scientific = FALSE,
      trim = TRUE
    ),
    "\n",
    sep = ""
  )

  cat(
    "Library size mediana: ",
    format(
      median(
        library_sizes
      ),
      scientific = FALSE,
      trim = TRUE
    ),
    "\n",
    sep = ""
  )

  cat(
    "Library size máxima: ",
    format(
      max(
        library_sizes
      ),
      scientific = FALSE,
      trim = TRUE
    ),
    "\n",
    sep = ""
  )

  cat(
    "Desvio máximo de soma CPM em relação a 1e6: ",
    format(
      max_cpm_sum_deviation,
      scientific = TRUE
    ),
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.15 AUDITAR PERFIS POR TIPO CELULAR
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "14. QC PSEUDOBULK POR TIPO CELULAR\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  pseudobulk_type_summary <- pseudobulk_meta[
    ,
    .(
      # Casts explícitos evitam que data.table encontre tipos diferentes
      # entre grupos (por exemplo, median() pode retornar integer em um
      # grupo e double em outro).
      n_pseudobulks = as.integer(.N),
      n_donors = as.integer(
        uniqueN(
          donor_id
        )
      ),
      total_nuclei = as.numeric(
        sum(
          n_nuclei
        )
      ),
      min_nuclei = as.numeric(
        min(
          n_nuclei
        )
      ),
      median_nuclei = as.numeric(
        median(
          n_nuclei
        )
      ),
      max_nuclei = as.numeric(
        max(
          n_nuclei
        )
      ),
      min_library_size = as.numeric(
        min(
          library_size
        )
      ),
      median_library_size = as.numeric(
        median(
          library_size
        )
      ),
      max_library_size = as.numeric(
        max(
          library_size
        )
      )
    ),
    by = cell_type
  ]

  pseudobulk_type_summary[
    ,
    cell_type_order :=
      match(
        cell_type,
        CELL_TYPE_ORDER
      )
  ]

  setorder(
    pseudobulk_type_summary,
    cell_type_order
  )

  pseudobulk_type_summary[
    ,
    cell_type_order := NULL
  ]

  fwrite(
    pseudobulk_type_summary,
    file.path(
      result_dir,
      "03_pseudobulk_celltype_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(
    pseudobulk_type_summary
  )


  # ============================================================================
  # 4.16 AUDITAR OVERLAP ROSMAP
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "15. OVERLAP COM BULK ROSMAP\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  overlap_donors <- unique(
    cell_meta[
      overlap_with_bulk_ROSMAP == TRUE,
      donor_id
    ]
  )

  nonoverlap_donors <- unique(
    cell_meta[
      overlap_with_bulk_ROSMAP == FALSE,
      donor_id
    ]
  )

  n_overlap_donors <- length(
    overlap_donors
  )

  n_nonoverlap_donors <- length(
    nonoverlap_donors
  )

  sensitivity_summary <- data.table(
    reference_version = c(
      "main_pathology_free",
      "sensitivity_no_ROSMAP_overlap"
    ),
    n_donors = c(
      uniqueN(
        cell_meta$donor_id
      ),
      n_nonoverlap_donors
    ),
    n_nuclei = c(
      nrow(
        cell_meta
      ),
      cell_meta[
        overlap_with_bulk_ROSMAP == FALSE,
        .N
      ]
    ),
    note = c(
      "Referência principal; inclui todos os 24 no-pathology.",
      "Definição preparada para análise de sensibilidade; não foi usada para substituir a referência principal nesta etapa."
    )
  )

  fwrite(
    sensitivity_summary,
    file.path(
      result_dir,
      "03_reference_sensitivity_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  cat(
    "Doadores overlap com bulk ROSMAP: ",
    n_overlap_donors,
    "\n",
    sep = ""
  )

  cat(
    "Doadores disponíveis para sensibilidade sem overlap: ",
    n_nonoverlap_donors,
    "\n",
    sep = ""
  )

  if (
    n_overlap_donors != 3L
  ) {

    warning(
      "Esperávamos 3 doadores com overlap conforme etapa 01; encontrados: ",
      n_overlap_donors,
      "."
    )
  }


  # ============================================================================
  # 4.17 SALVAR MATRIZES PSEUDOBULK E METADATA
  # ==============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "16. SALVAR REFERÊNCIA PSEUDOBULK\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  fwrite(
    pseudobulk_meta,
    pseudobulk_metadata_file,
    sep = "\t",
    na = "NA"
  )

  saveRDS(
    pseudobulk_counts,
    pseudobulk_counts_file,
    compress = FALSE
  )

  saveRDS(
    pseudobulk_cpm,
    pseudobulk_cpm_file,
    compress = FALSE
  )

  cat(
    "Raw counts:\n  ",
    pseudobulk_counts_file,
    "\n",
    sep = ""
  )

  cat(
    "CPM:\n  ",
    pseudobulk_cpm_file,
    "\n",
    sep = ""
  )

  cat(
    "Metadata:\n  ",
    pseudobulk_metadata_file,
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.18 RESUMO GLOBAL
  # ==============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "17. RESUMO GLOBAL\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  build_summary <- data.table(
    metric = c(
      "raw_gene_rows",
      "resolved_gene_rows",
      "unique_geneids_before_zero_filter",
      "all_zero_geneids_removed",
      "final_geneids",
      "raw_nuclei",
      "pathology_free_nuclei",
      "pathology_free_donors",
      "cell_types",
      "observed_pseudobulk_profiles",
      "missing_donor_celltype_combinations",
      "counts_total_before_gene_harmonization",
      "counts_removed_unresolved_or_ambiguous",
      "counts_removed_percent",
      "counts_final_harmonized",
      "ROSMAP_overlap_donors",
      "nonoverlap_donors_for_sensitivity"
    ),
    value = c(
      EXPECTED_RAW_GENE_ROWS,
      nrow(
        resolved_map
      ),
      EXPECTED_UNIQUE_RESOLVED_GENEIDS,
      n_zero_geneids,
      nrow(
        geneid_counts
      ),
      EXPECTED_RAW_NUCLEI,
      n_selected_nuclei,
      n_selected_donors,
      length(
        observed_cell_types
      ),
      n_pseudobulks,
      nrow(
        absent_profiles
      ),
      total_counts_selected_all_genes,
      unresolved_counts_total,
      round(
        unresolved_counts_percent,
        8
      ),
      total_nucleus_counts_final,
      n_overlap_donors,
      n_nonoverlap_donors
    )
  )

  fwrite(
    build_summary,
    file.path(
      result_dir,
      "03_reference_build_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(
    build_summary
  )


  # ============================================================================
  # 4.19 CHECKS FINAIS
  # ==============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "18. CHECKS FINAIS\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  checks <- rbindlist(
    list(

      make_check(
        "pathology_free_nuclei",
        ncol(
          geneid_counts
        ),
        EXPECTED_PATHOLOGY_FREE_NUCLEI,
        if (
          ncol(
            geneid_counts
          ) ==
            EXPECTED_PATHOLOGY_FREE_NUCLEI
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "pathology_free_donors",
        uniqueN(
          cell_meta$donor_id
        ),
        EXPECTED_PATHOLOGY_FREE_DONORS,
        if (
          uniqueN(
            cell_meta$donor_id
          ) ==
            EXPECTED_PATHOLOGY_FREE_DONORS
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "eight_cell_types_present",
        setequal(
          unique(
            cell_meta$cell_type
          ),
          EXPECTED_CELL_TYPES
        ),
        TRUE,
        if (
          setequal(
            unique(
              cell_meta$cell_type
            ),
            EXPECTED_CELL_TYPES
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "gene_collapse_count_conservation",
        abs(
          total_counts_after_geneid_collapse -
            total_resolved_matrix
        ) <=
          COUNT_TOLERANCE,
        TRUE,
        if (
          abs(
            total_counts_after_geneid_collapse -
              total_resolved_matrix
          ) <=
            COUNT_TOLERANCE
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "zero_gene_filter_count_conservation",
        abs(
          total_counts_after_zero_filter -
            total_counts_after_geneid_collapse
        ) <=
          COUNT_TOLERANCE,
        TRUE,
        if (
          abs(
            total_counts_after_zero_filter -
              total_counts_after_geneid_collapse
          ) <=
            COUNT_TOLERANCE
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "pseudobulk_count_conservation",
        abs(
          total_pseudobulk_counts -
            total_nucleus_counts_final
        ) <=
          COUNT_TOLERANCE,
        TRUE,
        if (
          abs(
            total_pseudobulk_counts -
              total_nucleus_counts_final
          ) <=
            COUNT_TOLERANCE
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "observed_pseudobulk_profiles",
        n_pseudobulks,
        EXPECTED_PSEUDOBULK_PROFILES,
        if (
          n_pseudobulks ==
            EXPECTED_PSEUDOBULK_PROFILES
        ) {
          "PASS"
        } else {
          "FAIL"
        },
        "Combinações donor x cell type ausentes não são criadas artificialmente."
      ),

      make_check(
        "pseudobulk_library_sizes_positive",
        all(
          library_sizes > 0
        ),
        TRUE,
        if (
          all(
            library_sizes > 0
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "CPM_columns_sum_to_1e6",
        max_cpm_sum_deviation <=
          CPM_SUM_TOLERANCE,
        TRUE,
        if (
          max_cpm_sum_deviation <=
            CPM_SUM_TOLERANCE
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "ROSMAP_overlap_donors",
        n_overlap_donors,
        3L,
        if (
          n_overlap_donors ==
            3L
        ) {
          "PASS"
        } else {
          "WARN"
        },
        "Overlap não invalida a referência principal; será avaliado por análise de sensibilidade."
      )
    ),
    fill = TRUE
  )

  fwrite(
    checks,
    file.path(
      result_dir,
      "03_reference_build_status.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(
    checks
  )

  if (
    any(
      checks$status == "FAIL",
      na.rm = TRUE
    )
  ) {

    stop(
      "Existe pelo menos um check FAIL após a construção da referência."
    )
  }


  # ============================================================================
  # 4.20 MANIFESTO DOS OUTPUTS
  # ==============================================================================

  output_manifest <- data.table(
    object = c(
      "nucleus_geneid_counts",
      "gene_metadata",
      "aligned_cell_metadata",
      "pseudobulk_raw_counts",
      "pseudobulk_cpm",
      "pseudobulk_metadata"
    ),
    path = c(
      nucleus_counts_file,
      gene_metadata_output_file,
      cell_metadata_output_file,
      pseudobulk_counts_file,
      pseudobulk_cpm_file,
      pseudobulk_metadata_file
    ),
    rows = c(
      nrow(
        geneid_counts
      ),
      nrow(
        gene_metadata
      ),
      nrow(
        cell_meta
      ),
      nrow(
        pseudobulk_counts
      ),
      nrow(
        pseudobulk_cpm
      ),
      nrow(
        pseudobulk_meta
      )
    ),
    columns = c(
      ncol(
        geneid_counts
      ),
      ncol(
        gene_metadata
      ),
      ncol(
        cell_meta
      ),
      ncol(
        pseudobulk_counts
      ),
      ncol(
        pseudobulk_cpm
      ),
      ncol(
        pseudobulk_meta
      )
    )
  )

  fwrite(
    output_manifest,
    file.path(
      result_dir,
      "03_output_manifest.tsv"
    ),
    sep = "\t",
    na = "NA"
  )


  # ============================================================================
  # 4.21 SESSION INFO
  # ==============================================================================

  sink(
    session_file
  )

  print(
    sessionInfo()
  )

  sink()


  # ============================================================================
  # 4.22 CONCLUSÃO
  # ==============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "19. CONCLUSÃO\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  cat(
    "Referência nucleus-level:\n  ",
    nrow(
      geneid_counts
    ),
    " GeneIDs x ",
    ncol(
      geneid_counts
    ),
    " núcleos\n",
    sep = ""
  )

  cat(
    "Referência pseudobulk:\n  ",
    nrow(
      pseudobulk_counts
    ),
    " GeneIDs x ",
    ncol(
      pseudobulk_counts
    ),
    " perfis donor x cell type\n",
    sep = ""
  )

  cat(
    "Doadores: ",
    n_selected_donors,
    "\n",
    sep = ""
  )

  cat(
    "Tipos celulares: ",
    length(
      observed_cell_types
    ),
    "\n",
    sep = ""
  )

  cat(
    "Counts removidos por genes não resolvidos/ambíguos: ",
    round(
      unresolved_counts_percent,
      6
    ),
    "%\n",
    sep = ""
  )

  cat(
    "GeneIDs all-zero removidos: ",
    n_zero_geneids,
    "\n",
    sep = ""
  )

  cat(
    "Overlap com bulk ROSMAP: ",
    n_overlap_donors,
    " doadores\n",
    sep = ""
  )

  cat(
    "Sensibilidade futura sem overlap: ",
    n_nonoverlap_donors,
    " doadores\n",
    sep = ""
  )

  cat(
    "\nNenhum pseudobulk artificial de zeros foi criado.\n"
  )

  cat(
    "Raw counts nucleus-level e pseudobulk foram preservados.\n"
  )

  cat(
    "CPM foi calculado somente para os perfis pseudobulk.\n"
  )

  cat(
    "\nPróxima etapa:\n"
  )

  cat(
    "  04_derive_markers_Mathys2019.R\n"
  )

  invisible(TRUE)
}


# ==============================================================================
# 5. EXECUÇÃO COM LOG
# ==============================================================================

log_connection <- file(
  log_file,
  open = "wt"
)

sink(
  log_connection,
  split = TRUE
)

execution_error <- NULL

tryCatch(
  {
    run_build_reference()
  },
  error = function(e) {

    execution_error <<- e

    cat(
      "\n============================================================\n"
    )

    cat(
      "ERRO NA CONSTRUÇÃO DA REFERÊNCIA\n"
    )

    cat(
      "============================================================\n"
    )

    cat(
      conditionMessage(
        e
      ),
      "\n"
    )
  }
)

sink()

close(
  log_connection
)

if (!is.null(execution_error)) {

  stop(
    conditionMessage(
      execution_error
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 6. MENSAGEM FINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "03_build_reference_Mathys2019.R CONCLUÍDO\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Log:\n  ",
  log_file,
  "\n\n",
  sep = ""
)

cat(
  "Resumo:\n  ",
  file.path(
    result_dir,
    "03_reference_build_summary.tsv"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Status:\n  ",
  file.path(
    result_dir,
    "03_reference_build_status.tsv"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Referência nucleus-level:\n  ",
  nucleus_counts_file,
  "\n\n",
  sep = ""
)

cat(
  "Pseudobulk raw counts:\n  ",
  pseudobulk_counts_file,
  "\n\n",
  sep = ""
)

cat(
  "Pseudobulk CPM:\n  ",
  pseudobulk_cpm_file,
  "\n\n",
  sep = ""
)

cat(
  "Pseudobulk metadata:\n  ",
  pseudobulk_metadata_file,
  "\n",
  sep = ""
)
