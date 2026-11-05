#!/usr/bin/env Rscript

# ==============================================================================
# 04_derive_markers_Mathys2019.R
# ==============================================================================
#
# FASE 02 — DECONVOLUÇÃO CELULAR
# ETAPA 04 — DERIVAÇÃO DE MARCADORES DA REFERÊNCIA MATHYS 2019
#
# Objetivo
# --------
# Derivar marcadores celulares diretamente da referência pseudobulk construída
# no script 03, usando os 24 doadores no-pathology como replicatas biológicas.
#
# O desenho definido para esta etapa é:
#
#   - expressão = CPM pseudobulk;
#   - unidade de replicação = doador;
#   - 8 tipos celulares:
#       Ex, In, Ast, Oli, Opc, Mic, End, Per;
#   - expressão média calculada separadamente por tipo celular;
#   - especificidade definida como a fração do sinal médio do gene atribuída
#     a cada tipo celular;
#   - ranking determinístico por:
#       1) specificity decrescente;
#       2) mean_CPM decrescente;
#       3) GeneID crescente;
#   - conjunto principal:
#       top 300 marcadores por tipo;
#   - análises de sensibilidade:
#       top 100 e top 500 por tipo.
#
#
# DEFINIÇÃO DE ESPECIFICIDADE
# ---------------------------
# Para um gene g e um tipo celular t:
#
#                  mean_CPM(g,t)
# specificity = -----------------------
#               sum_k mean_CPM(g,k)
#
# em que k percorre os 8 tipos celulares.
#
# Para cada gene com expressão positiva em pelo menos um tipo:
#
#     sum_t specificity(g,t) = 1
#
#
# IMPORTANTE SOBRE DOADORES AUSENTES
# ----------------------------------
# End e Per não estão presentes em todos os 24 doadores:
#
#   End: 19 doadores
#   Per: 18 doadores
#
# Um donor x cell type ausente NÃO é convertido em expressão zero.
#
# Portanto:
#
#   mean_CPM(g, End)
#
# é calculado usando apenas os pseudobulks End realmente observados.
#
# Isso evita interpretar ausência de amostragem daquele tipo celular como
# ausência biológica de expressão.
#
#
# ROBUSTEZ DOS MARCADORES
# -----------------------
# As execuções exploratórias mostraram que specificity isolada favorecia genes
# de expressão muito baixa, enquanto um cutoff rígido de detection_fraction
# eliminava marcadores plausíveis de tipos raros.
#
# A seleção final combina três dimensões de forma contínua:
#
#   marker_score =
#       specificity *
#       log1p(mean_CPM) *
#       detection_fraction
#
# O denominador de detection_fraction é o número de doadores em que aquele tipo
# celular foi realmente observado; donor x cell type ausente NÃO entra como zero.
#
# Ranking final:
#   1) marker_score decrescente;
#   2) specificity decrescente;
#   3) mean_CPM decrescente;
#   4) GeneID crescente.
#
# Nenhum cutoff rígido adicional de mean_CPM ou detection_fraction é imposto.
# O marker_score é uma regra operacional deste projeto, não uma estatística
# inferencial nem um método universal de seleção de marcadores.
#
#
# O SCRIPT NÃO
# ------------
#   - usa informação de AD/HD/PD para escolher os marcadores;
#   - usa expressão bulk na derivação dos marcadores;
#   - remove genes por nome;
#   - exclui genes mitocondriais/ribossomais arbitrariamente;
#   - impõe cutoff rígido adicional de prevalência ou mean_CPM;
#   - trata donor x cell type ausente como zero;
#   - roda deconvolução.
#
#
# INPUTS
# ------
# data/reference/mathys2019/phase02_prepared/reference/
#   mathys2019_pathology_free_pseudobulk_cpm_sparse.rds
#   mathys2019_pathology_free_pseudobulk_metadata.tsv
#   mathys2019_pathology_free_gene_metadata.tsv
#
#
# OUTPUTS PRINCIPAIS
# ------------------
# data/reference/mathys2019/phase02_prepared/markers/
#
#   mathys2019_celltype_mean_cpm.tsv
#   mathys2019_celltype_specificity.tsv
#   mathys2019_marker_metrics_long.tsv
#
#   mathys2019_top100_markers.tsv
#   mathys2019_top300_markers.tsv
#   mathys2019_top500_markers.tsv
#
#   mathys2019_signature_top100_mean_cpm.tsv
#   mathys2019_signature_top300_mean_cpm.tsv
#   mathys2019_signature_top500_mean_cpm.tsv
#
# O conjunto top300 é o principal. Top100 e top500 são sensibilidades.
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
    "\nInstale somente os pacotes ausentes no ambiente R atual. ",
    "Não altere o renv automaticamente."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})


# ==============================================================================
# 1. PARÂMETROS PREDEFINIDOS
# ==============================================================================

EXPECTED_GENEIDS <- 17394L
EXPECTED_PSEUDOBULKS <- 181L
EXPECTED_DONORS <- 24L

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

MARKER_SET_SIZES <- c(
  100L,
  300L,
  500L
)

PRIMARY_MARKER_SET_SIZE <- 300L

# CPM > 0 equivale a count > 0 para esses pseudobulks, pois o CPM apenas
# reescala cada coluna por um fator positivo.
DETECTION_THRESHOLD_CPM <- 0

# Não há cutoff rígido de detection_fraction na seleção final.
MIN_POSITIVE_DETECTION_FRACTION <- 0

SPECIFICITY_SUM_TOLERANCE <- 1e-8
CPM_COLUMN_SUM_TOLERANCE <- 1e-4


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

prepared_dir <- file.path(
  "data",
  "reference",
  "mathys2019",
  "phase02_prepared"
)

reference_dir <- file.path(
  prepared_dir,
  "reference"
)

markers_dir <- file.path(
  prepared_dir,
  "markers"
)

result_dir <- file.path(
  "results",
  "deconvolution",
  "mathys2019",
  "marker_derivation"
)

dir.create(
  markers_dir,
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

pseudobulk_cpm_file <- file.path(
  reference_dir,
  "mathys2019_pathology_free_pseudobulk_cpm_sparse.rds"
)

pseudobulk_metadata_file <- file.path(
  reference_dir,
  "mathys2019_pathology_free_pseudobulk_metadata.tsv"
)

gene_metadata_file <- file.path(
  reference_dir,
  "mathys2019_pathology_free_gene_metadata.tsv"
)

step03_status_file <- file.path(
  "results",
  "deconvolution",
  "mathys2019",
  "reference_build",
  "03_reference_build_status.tsv"
)


# ------------------------------------------------------------------------------
# Outputs principais
# ------------------------------------------------------------------------------

mean_cpm_file <- file.path(
  markers_dir,
  "mathys2019_celltype_mean_cpm.tsv"
)

specificity_file <- file.path(
  markers_dir,
  "mathys2019_celltype_specificity.tsv"
)

marker_metrics_file <- file.path(
  markers_dir,
  "mathys2019_marker_metrics_long.tsv"
)

top100_file <- file.path(
  markers_dir,
  "mathys2019_top100_markers.tsv"
)

top300_file <- file.path(
  markers_dir,
  "mathys2019_top300_markers.tsv"
)

top500_file <- file.path(
  markers_dir,
  "mathys2019_top500_markers.tsv"
)

signature100_file <- file.path(
  markers_dir,
  "mathys2019_signature_top100_mean_cpm.tsv"
)

signature300_file <- file.path(
  markers_dir,
  "mathys2019_signature_top300_mean_cpm.tsv"
)

signature500_file <- file.path(
  markers_dir,
  "mathys2019_signature_top500_mean_cpm.tsv"
)

log_file <- file.path(
  result_dir,
  "04_derive_markers_Mathys2019.txt"
)

session_file <- file.path(
  result_dir,
  "04_derive_markers_Mathys2019_sessionInfo.txt"
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


# ------------------------------------------------------------------------------
# Construir um data.table de uma matriz densa com GeneID como primeira coluna.
# ------------------------------------------------------------------------------

matrix_to_gene_table <- function(
  x,
  gene_ids
) {

  if (
    nrow(x) !=
      length(gene_ids)
  ) {

    stop(
      "matrix_to_gene_table: número de linhas diferente do número de GeneIDs."
    )
  }

  out <- as.data.table(
    x
  )

  out[
    ,
    gene_id := gene_ids
  ]

  setcolorder(
    out,
    c(
      "gene_id",
      setdiff(
        names(out),
        "gene_id"
      )
    )
  )

  out
}


# ------------------------------------------------------------------------------
# Criar conjunto top-K de forma determinística.
#
# Ranking:
#   specificity DESC
#   mean_CPM DESC
#   GeneID numérico ASC
#
# O gene precisa ter expressão, especificidade, detecção e marker_score
# positivos. O ranking usa marker_score, specificity, mean_CPM e GeneID.
# ------------------------------------------------------------------------------

derive_top_k <- function(
  metrics,
  k
) {

  out_list <- vector(
    "list",
    length(
      CELL_TYPE_ORDER
    )
  )

  names(
    out_list
  ) <- CELL_TYPE_ORDER

  for (
    cell_type_i in
    CELL_TYPE_ORDER
  ) {

    x <- copy(
      metrics[
        cell_type ==
          cell_type_i &
          mean_CPM > 0 &
          specificity > 0 &
          detection_fraction >
            MIN_POSITIVE_DETECTION_FRACTION &
          marker_score > 0
      ]
    )

    if (
      nrow(x) <
        k
    ) {

      stop(
        "Tipo celular ",
        cell_type_i,
        " possui apenas ",
        nrow(x),
        " genes candidatos com marker_score > 0; são necessários ",
        k,
        " para o conjunto top",
        k,
        "."
      )
    }

    x[
      ,
      gene_id_numeric :=
        as.numeric(
          gene_id
        )
    ]

    if (
      anyNA(
        x$gene_id_numeric
      )
    ) {

      stop(
        "GeneID não numérico durante o ranking de ",
        cell_type_i,
        "."
      )
    }

    setorder(
      x,
      -marker_score,
      -specificity,
      -mean_CPM,
      gene_id_numeric
    )

    x[
      ,
      marker_rank := seq_len(.N)
    ]

    x <- x[
      marker_rank <=
        k
    ]

    x[
      ,
      marker_set_size :=
        as.integer(
          k
        )
    ]

    x[
      ,
      gene_id_numeric := NULL
    ]

    out_list[[cell_type_i]] <- x
  }

  out <- rbindlist(
    out_list,
    use.names = TRUE,
    fill = TRUE
  )

  out[
    ,
    cell_type_order :=
      match(
        cell_type,
        CELL_TYPE_ORDER
      )
  ]

  setorder(
    out,
    cell_type_order,
    marker_rank
  )

  out[
    ,
    cell_type_order := NULL
  ]

  out
}


# ------------------------------------------------------------------------------
# Criar signature matrix restrita à união dos marcadores top-K.
# ------------------------------------------------------------------------------

make_signature_table <- function(
  top_markers,
  mean_matrix,
  gene_ids
) {

  marker_gene_ids <- unique(
    top_markers$gene_id
  )

  idx <- match(
    marker_gene_ids,
    gene_ids
  )

  if (
    anyNA(idx)
  ) {

    stop(
      "Há marcador sem correspondência na matriz média."
    )
  }

  signature_matrix <- mean_matrix[
    idx,
    ,
    drop = FALSE
  ]

  out <- matrix_to_gene_table(
    signature_matrix,
    marker_gene_ids
  )

  out
}


# ------------------------------------------------------------------------------
# Auditoria de overlap entre listas de marcadores de tipos celulares.
# ------------------------------------------------------------------------------

pairwise_marker_overlap <- function(
  top_markers,
  k
) {

  rows <- list()
  counter <- 1L

  for (
    i in seq_along(
      CELL_TYPE_ORDER
    )
  ) {

    type_i <- CELL_TYPE_ORDER[i]

    genes_i <- unique(
      top_markers[
        cell_type ==
          type_i,
        gene_id
      ]
    )

    for (
      j in i:length(
        CELL_TYPE_ORDER
      )
    ) {

      type_j <- CELL_TYPE_ORDER[j]

      genes_j <- unique(
        top_markers[
          cell_type ==
            type_j,
        gene_id
      ]
    )

      n_intersection <- length(
        intersect(
          genes_i,
          genes_j
        )
      )

      n_union <- length(
        union(
          genes_i,
          genes_j
        )
      )

      jaccard <- if (
        n_union > 0L
      ) {
        n_intersection /
          n_union
      } else {
        NA_real_
      }

      rows[[counter]] <- data.table(
        marker_set_size = as.integer(
          k
        ),
        cell_type_1 = type_i,
        cell_type_2 = type_j,
        n_markers_1 = as.integer(
          length(
            genes_i
          )
        ),
        n_markers_2 = as.integer(
          length(
            genes_j
          )
        ),
        overlap_n = as.integer(
          n_intersection
        ),
        jaccard = as.numeric(
          jaccard
        )
      )

      counter <- counter + 1L
    }
  }

  rbindlist(
    rows,
    use.names = TRUE,
    fill = TRUE
  )
}


# ==============================================================================
# 4. FUNÇÃO PRINCIPAL
# ==============================================================================

run_marker_derivation <- function() {


  # ============================================================================
  # 4.1 CABEÇALHO E PRÉ-REQUISITOS
  # ============================================================================

  cat(
    "============================================================\n"
  )

  cat(
    "FASE 02 — 04_derive_markers_Mathys2019.R\n"
  )

  cat(
    "DERIVAÇÃO DE MARCADORES CELULARES\n"
  )

  cat(
    "============================================================\n\n"
  )

  validate_status_file(
    step03_status_file,
    "Etapa 03"
  )

  cat(
    "Etapa 03 validada: nenhum check FAIL.\n"
  )

  required_files <- c(
    pseudobulk_cpm = pseudobulk_cpm_file,
    pseudobulk_metadata = pseudobulk_metadata_file,
    gene_metadata = gene_metadata_file
  )

  missing_files <- required_files[
    !file.exists(
      required_files
    )
  ]

  if (
    length(
      missing_files
    ) > 0L
  ) {

    cat(
      "\nArquivos obrigatórios ausentes:\n"
    )

    print(
      missing_files
    )

    stop(
      "Derivação de marcadores interrompida."
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
      "04_input_file_audit.tsv"
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
  # 4.2 CARREGAR PSEUDOBULK CPM E METADATA
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "1. CARREGAMENTO DA REFERÊNCIA PSEUDOBULK\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  cpm <- readRDS(
    pseudobulk_cpm_file
  )

  if (
    !inherits(
      cpm,
      "sparseMatrix"
    )
  ) {

    stop(
      "Objeto CPM não é sparseMatrix."
    )
  }

  cpm <- as(
    cpm,
    "dgCMatrix"
  )

  pseudobulk_meta <- fread(
    pseudobulk_metadata_file,
    na.strings = c(
      "NA",
      ""
    )
  )

  gene_meta <- fread(
    gene_metadata_file,
    na.strings = c(
      "NA",
      ""
    )
  )

  required_pb_columns <- c(
    "pseudobulk_id",
    "donor_id",
    "cell_type",
    "n_nuclei",
    "library_size"
  )

  missing_pb_columns <- setdiff(
    required_pb_columns,
    names(
      pseudobulk_meta
    )
  )

  if (
    length(
      missing_pb_columns
    ) > 0L
  ) {

    stop(
      "Pseudobulk metadata sem colunas obrigatórias: ",
      paste(
        missing_pb_columns,
        collapse = ", "
      )
    )
  }

  if (
    !"gene_id" %in%
      names(
        gene_meta
      )
  ) {

    stop(
      "Gene metadata sem coluna gene_id."
    )
  }

  pseudobulk_meta[
    ,
    pseudobulk_id :=
      clean_chr(
        pseudobulk_id
      )
  ]

  pseudobulk_meta[
    ,
    donor_id :=
      clean_chr(
        donor_id
      )
  ]

  pseudobulk_meta[
    ,
    cell_type :=
      clean_chr(
        cell_type
      )
  ]

  gene_meta[
    ,
    gene_id :=
      clean_chr(
        gene_id
      )
  ]

  if (
    nrow(
      cpm
    ) !=
      EXPECTED_GENEIDS
  ) {

    stop(
      "Esperávamos ",
      EXPECTED_GENEIDS,
      " GeneIDs; encontrados: ",
      nrow(
        cpm
      ),
      "."
    )
  }

  if (
    ncol(
      cpm
    ) !=
      EXPECTED_PSEUDOBULKS
  ) {

    stop(
      "Esperávamos ",
      EXPECTED_PSEUDOBULKS,
      " pseudobulks; encontrados: ",
      ncol(
        cpm
      ),
      "."
    )
  }

  if (
    nrow(
      pseudobulk_meta
    ) !=
      EXPECTED_PSEUDOBULKS
  ) {

    stop(
      "Pseudobulk metadata não possui 181 linhas."
    )
  }

  if (
    nrow(
      gene_meta
    ) !=
      EXPECTED_GENEIDS
  ) {

    stop(
      "Gene metadata não possui 17.394 linhas."
    )
  }

  if (
    is.null(
      rownames(
        cpm
      )
    ) |
      is.null(
        colnames(
          cpm
        )
      )
  ) {

    stop(
      "Matriz CPM precisa possuir rownames GeneID e colnames pseudobulk_id."
    )
  }

  # Alinhar metadata à ordem da matriz.
  pb_idx <- match(
    colnames(
      cpm
    ),
    pseudobulk_meta$pseudobulk_id
  )

  if (
    anyNA(
      pb_idx
    )
  ) {

    stop(
      "Existem colunas da matriz CPM sem correspondência no pseudobulk metadata."
    )
  }

  pseudobulk_meta <- pseudobulk_meta[
    pb_idx
  ]

  if (
    !identical(
      pseudobulk_meta$pseudobulk_id,
      colnames(
        cpm
      )
    )
  ) {

    stop(
      "Falha ao alinhar pseudobulk metadata às colunas da matriz CPM."
    )
  }

  gene_idx <- match(
    rownames(
      cpm
    ),
    gene_meta$gene_id
  )

  if (
    anyNA(
      gene_idx
    )
  ) {

    stop(
      "Existem GeneIDs da matriz CPM sem correspondência no gene metadata."
    )
  }

  gene_meta <- gene_meta[
    gene_idx
  ]

  if (
    !identical(
      gene_meta$gene_id,
      rownames(
        cpm
      )
    )
  ) {

    stop(
      "Falha ao alinhar gene metadata às linhas da matriz CPM."
    )
  }

  if (
    any(
      !is.finite(
        cpm@x
      )
    )
  ) {

    stop(
      "Matriz CPM contém valores não finitos."
    )
  }

  if (
    any(
      cpm@x < 0
    )
  ) {

    stop(
      "Matriz CPM contém valores negativos."
    )
  }

  cpm_sums <- Matrix::colSums(
    cpm
  )

  max_cpm_deviation <- max(
    abs(
      cpm_sums -
        1e6
    )
  )

  if (
    max_cpm_deviation >
      CPM_COLUMN_SUM_TOLERANCE
  ) {

    stop(
      "As colunas CPM não somam aproximadamente 1e6. ",
      "Desvio máximo: ",
      max_cpm_deviation
    )
  }

  observed_types <- sort(
    unique(
      pseudobulk_meta$cell_type
    )
  )

  observed_donors <- uniqueN(
    pseudobulk_meta$donor_id
  )

  if (
    !setequal(
      observed_types,
      EXPECTED_CELL_TYPES
    )
  ) {

    stop(
      "Os 8 tipos celulares esperados não estão presentes."
    )
  }

  if (
    observed_donors !=
      EXPECTED_DONORS
  ) {

    stop(
      "Esperávamos 24 doadores na referência; encontrados: ",
      observed_donors,
      "."
    )
  }

  cat(
    "Matriz CPM: ",
    nrow(
      cpm
    ),
    " GeneIDs x ",
    ncol(
      cpm
    ),
    " pseudobulks\n",
    sep = ""
  )

  cat(
    "Doadores: ",
    observed_donors,
    "\n",
    sep = ""
  )

  cat(
    "Tipos celulares: ",
    paste(
      observed_types,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.3 SUPORTE DA REFERÊNCIA POR TIPO CELULAR
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "2. SUPORTE POR TIPO CELULAR\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  type_support_list <- vector(
    "list",
    length(
      CELL_TYPE_ORDER
    )
  )

  names(
    type_support_list
  ) <- CELL_TYPE_ORDER

  for (
    cell_type_i in
    CELL_TYPE_ORDER
  ) {

    dt_i <- pseudobulk_meta[
      cell_type ==
        cell_type_i
    ]

    type_support_list[[cell_type_i]] <- data.table(
      cell_type = cell_type_i,
      n_pseudobulks = as.integer(
        nrow(
          dt_i
        )
      ),
      n_donors = as.integer(
        uniqueN(
          dt_i$donor_id
        )
      ),
      total_nuclei = as.numeric(
        sum(
          dt_i$n_nuclei
        )
      ),
      min_nuclei = as.numeric(
        min(
          dt_i$n_nuclei
        )
      ),
      median_nuclei = as.numeric(
        median(
          dt_i$n_nuclei
        )
      ),
      max_nuclei = as.numeric(
        max(
          dt_i$n_nuclei
        )
      ),
      min_library_size = as.numeric(
        min(
          dt_i$library_size
        )
      ),
      median_library_size = as.numeric(
        median(
          dt_i$library_size
        )
      ),
      max_library_size = as.numeric(
        max(
          dt_i$library_size
        )
      ),
      rare_reference_type = cell_type_i %in%
        c(
          "End",
          "Per"
        )
    )
  }

  type_support <- rbindlist(
    type_support_list,
    use.names = TRUE,
    fill = TRUE
  )

  fwrite(
    type_support,
    file.path(
      result_dir,
      "04_celltype_reference_support.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(
    type_support
  )


  # ============================================================================
  # 4.4 CALCULAR MÉDIA, MEDIANA E DETECÇÃO POR TIPO
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "3. ESTATÍSTICAS DE EXPRESSÃO POR TIPO CELULAR\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  gene_ids <- rownames(
    cpm
  )

  n_genes <- length(
    gene_ids
  )

  n_types <- length(
    CELL_TYPE_ORDER
  )

  mean_cpm_matrix <- matrix(
    0,
    nrow = n_genes,
    ncol = n_types,
    dimnames = list(
      gene_ids,
      CELL_TYPE_ORDER
    )
  )

  median_cpm_matrix <- matrix(
    0,
    nrow = n_genes,
    ncol = n_types,
    dimnames = list(
      gene_ids,
      CELL_TYPE_ORDER
    )
  )

  detection_n_matrix <- matrix(
    0L,
    nrow = n_genes,
    ncol = n_types,
    dimnames = list(
      gene_ids,
      CELL_TYPE_ORDER
    )
  )

  detection_fraction_matrix <- matrix(
    0,
    nrow = n_genes,
    ncol = n_types,
    dimnames = list(
      gene_ids,
      CELL_TYPE_ORDER
    )
  )

  n_donors_by_type <- setNames(
    integer(
      n_types
    ),
    CELL_TYPE_ORDER
  )

  for (
    cell_type_i in
    CELL_TYPE_ORDER
  ) {

    idx <- which(
      pseudobulk_meta$cell_type ==
        cell_type_i
    )

    if (
      length(
        idx
      ) == 0L
    ) {

      stop(
        "Nenhum pseudobulk encontrado para ",
        cell_type_i,
        "."
      )
    }

    donor_ids_i <- pseudobulk_meta$donor_id[
      idx
    ]

    if (
      anyDuplicated(
        donor_ids_i
      )
    ) {

      stop(
        "Há mais de um pseudobulk do mesmo doador para ",
        cell_type_i,
        "."
      )
    }

    n_donors_i <- length(
      idx
    )

    n_donors_by_type[
      cell_type_i
    ] <- as.integer(
      n_donors_i
    )

    type_cpm <- cpm[
      ,
      idx,
      drop = FALSE
    ]

    # Cada donor x cell type observado recebe o mesmo peso.
    mean_i <- as.numeric(
      Matrix::rowMeans(
        type_cpm
      )
    )

    # Apenas para mediana convertemos o pequeno subconjunto por tipo em matriz
    # densa. O maior bloco possui 17.394 x 24, tamanho pequeno o suficiente.
    type_dense <- as.matrix(
      type_cpm
    )

    median_i <- as.numeric(
      apply(
        type_dense,
        1,
        stats::median
      )
    )

    detection_n_i <- as.integer(
      rowSums(
        type_dense >
          DETECTION_THRESHOLD_CPM
      )
    )

    detection_fraction_i <- as.numeric(
      detection_n_i /
        n_donors_i
    )

    mean_cpm_matrix[
      ,
      cell_type_i
    ] <- mean_i

    median_cpm_matrix[
      ,
      cell_type_i
    ] <- median_i

    detection_n_matrix[
      ,
      cell_type_i
    ] <- detection_n_i

    detection_fraction_matrix[
      ,
      cell_type_i
    ] <- detection_fraction_i

    rm(
      type_cpm,
      type_dense
    )

    invisible(
      gc()
    )

    cat(
      cell_type_i,
      ": ",
      n_donors_i,
      " doadores processados.\n",
      sep = ""
    )
  }


  # ============================================================================
  # 4.5 CALCULAR ESPECIFICIDADE
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "4. ESPECIFICIDADE POR GENE × TIPO CELULAR\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  total_mean_signal <- rowSums(
    mean_cpm_matrix
  )

  zero_total_signal <- which(
    total_mean_signal <= 0
  )

  if (
    length(
      zero_total_signal
    ) > 0L
  ) {

    stop(
      "Foram encontrados ",
      length(
        zero_total_signal
      ),
      " GeneIDs com soma de mean_CPM igual a zero. ",
      "Isso não era esperado após o filtro all-zero do script 03."
    )
  }

  specificity_matrix <- mean_cpm_matrix /
    total_mean_signal

  specificity_sum <- rowSums(
    specificity_matrix
  )

  max_specificity_sum_deviation <- max(
    abs(
      specificity_sum -
        1
    )
  )

  if (
    max_specificity_sum_deviation >
      SPECIFICITY_SUM_TOLERANCE
  ) {

    stop(
      "A soma das especificidades por gene difere de 1 acima da tolerância. ",
      "Desvio máximo: ",
      max_specificity_sum_deviation
    )
  }

  cat(
    "Desvio máximo de sum(specificity) em relação a 1: ",
    format(
      max_specificity_sum_deviation,
      scientific = TRUE
    ),
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.6 SALVAR MATRIZES MÉDIAS E DE ESPECIFICIDADE
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "5. MATRIZES DE REFERÊNCIA POR TIPO CELULAR\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  mean_cpm_dt <- matrix_to_gene_table(
    mean_cpm_matrix,
    gene_ids
  )

  specificity_dt <- matrix_to_gene_table(
    specificity_matrix,
    gene_ids
  )

  fwrite(
    mean_cpm_dt,
    mean_cpm_file,
    sep = "\t",
    na = "NA"
  )

  fwrite(
    specificity_dt,
    specificity_file,
    sep = "\t",
    na = "NA"
  )

  cat(
    "Matriz mean CPM salva:\n  ",
    mean_cpm_file,
    "\n",
    sep = ""
  )

  cat(
    "Matriz specificity salva:\n  ",
    specificity_file,
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.7 CONSTRUIR TABELA LONGA DE MÉTRICAS
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "6. TABELA LONGA DE MÉTRICAS DE MARCADOR\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  metric_list <- vector(
    "list",
    length(
      CELL_TYPE_ORDER
    )
  )

  names(
    metric_list
  ) <- CELL_TYPE_ORDER

  # Colunas úteis do gene metadata para anotação.
  gene_annotation_cols <- intersect(
    c(
      "gene_id",
      "original_symbols",
      "annotation_symbols",
      "annotation_gene_names",
      "mapping_methods",
      "collapsed_multiple_input_rows"
    ),
    names(
      gene_meta
    )
  )

  gene_annotation <- gene_meta[
    ,
    ..gene_annotation_cols
  ]

  for (
    cell_type_i in
    CELL_TYPE_ORDER
  ) {

    metrics_i <- data.table(
      gene_id = gene_ids,
      cell_type = cell_type_i,
      n_donors = as.integer(
        n_donors_by_type[
          cell_type_i
        ]
      ),
      n_donors_detected = as.integer(
        detection_n_matrix[
          ,
          cell_type_i
        ]
      ),
      detection_fraction = as.numeric(
        detection_fraction_matrix[
          ,
          cell_type_i
        ]
      ),
      mean_CPM = as.numeric(
        mean_cpm_matrix[
          ,
          cell_type_i
        ]
      ),
      median_CPM = as.numeric(
        median_cpm_matrix[
          ,
          cell_type_i
        ]
      ),
      specificity = as.numeric(
        specificity_matrix[
          ,
          cell_type_i
        ]
      ),
      marker_score = as.numeric(
        specificity_matrix[
          ,
          cell_type_i
        ] *
          log1p(
            mean_cpm_matrix[
              ,
              cell_type_i
            ]
          ) *
          detection_fraction_matrix[
            ,
            cell_type_i
          ]
      ),
      total_mean_CPM_across_types = as.numeric(
        total_mean_signal
      ),
      rare_reference_type = cell_type_i %in%
        c(
          "End",
          "Per"
        )
    )

    metrics_i <- merge(
      metrics_i,
      gene_annotation,
      by = "gene_id",
      all.x = TRUE,
      sort = FALSE
    )

    metrics_i <- metrics_i[
      match(
        gene_ids,
        gene_id
      )
    ]

    metric_list[[cell_type_i]] <- metrics_i
  }

  marker_metrics <- rbindlist(
    metric_list,
    use.names = TRUE,
    fill = TRUE
  )

  marker_metrics[
    ,
    cell_type_order :=
      match(
        cell_type,
        CELL_TYPE_ORDER
      )
  ]

  marker_metrics[
    ,
    gene_id_numeric :=
      as.numeric(
        gene_id
      )
  ]

  if (
    anyNA(
      marker_metrics$gene_id_numeric
    )
  ) {

    stop(
      "GeneID não numérico na tabela de métricas."
    )
  }

  setorder(
    marker_metrics,
    cell_type_order,
    gene_id_numeric
  )

  marker_metrics[
    ,
    c(
      "cell_type_order",
      "gene_id_numeric"
    ) := NULL
  ]

  fwrite(
    marker_metrics,
    marker_metrics_file,
    sep = "\t",
    na = "NA"
  )

  cat(
    "Linhas gene x tipo celular: ",
    nrow(
      marker_metrics
    ),
    "\n",
    sep = ""
  )

  expected_metric_rows <- EXPECTED_GENEIDS *
    length(
      EXPECTED_CELL_TYPES
    )

  if (
    nrow(
      marker_metrics
    ) !=
      expected_metric_rows
  ) {

    stop(
      "Número inesperado de linhas na tabela de métricas."
    )
  }


  # ============================================================================
  # 4.8 AUDITAR CANDIDATOS E marker_score
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "7. AUDITORIA DOS CANDIDATOS E marker_score\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  candidate_list <- vector(
    "list",
    length(CELL_TYPE_ORDER)
  )

  names(candidate_list) <- CELL_TYPE_ORDER

  for (cell_type_i in CELL_TYPE_ORDER) {

    x <- marker_metrics[
      cell_type == cell_type_i &
        mean_CPM > 0 &
        specificity > 0 &
        detection_fraction > MIN_POSITIVE_DETECTION_FRACTION &
        marker_score > 0
    ]

    candidate_list[[cell_type_i]] <- data.table(
      cell_type = cell_type_i,
      n_reference_donors = as.integer(unique(x$n_donors)),
      candidate_genes = as.integer(nrow(x)),
      median_mean_CPM = as.numeric(median(x$mean_CPM)),
      median_specificity = as.numeric(median(x$specificity)),
      median_detection_fraction = as.numeric(
        median(x$detection_fraction)
      ),
      median_marker_score = as.numeric(
        median(x$marker_score)
      ),
      enough_for_top500 = nrow(x) >= max(MARKER_SET_SIZES),
      rare_reference_type = cell_type_i %in% c("End", "Per")
    )
  }

  candidate_summary <- rbindlist(
    candidate_list,
    use.names = TRUE,
    fill = TRUE
  )

  fwrite(
    candidate_summary,
    file.path(
      result_dir,
      "04_marker_candidates_by_celltype.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(candidate_summary)

  if (any(!candidate_summary$enough_for_top500)) {
    stop(
      "Pelo menos um tipo celular possui menos de ",
      max(MARKER_SET_SIZES),
      " genes candidatos com marker_score > 0."
    )
  }

  marker_score_parameters <- data.table(
    parameter = c(
      "score_formula",
      "detection_denominator",
      "hard_detection_cutoff",
      "hard_mean_CPM_cutoff",
      "primary_marker_set",
      "sensitivity_marker_sets"
    ),
    value = c(
      "specificity * log1p(mean_CPM) * detection_fraction",
      "observed donors for that cell type",
      "none beyond >0 for candidate definition",
      "none beyond >0 for candidate definition",
      "top300 per cell type",
      "top100 and top500 per cell type"
    )
  )

  fwrite(
    marker_score_parameters,
    file.path(
      result_dir,
      "04_marker_score_parameters.tsv"
    ),
    sep = "\t",
    na = "NA"
  )


  # 4.9 DERIVAR TOP100 / TOP300 / TOP500
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "8. RANKING E CONJUNTOS DE MARCADORES\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  top_sets <- list()

  for (
    k in
    MARKER_SET_SIZES
  ) {

    top_k <- derive_top_k(
      marker_metrics,
      k
    )

    key <- paste0(
      "top",
      k
    )

    top_sets[[key]] <- top_k

    counts_by_type <- top_k[
      ,
      .(
        n_markers = .N
      ),
      by = cell_type
    ]

    if (
      any(
        counts_by_type$n_markers !=
          k
      )
    ) {

      stop(
        "Conjunto top",
        k,
        " não contém exatamente ",
        k,
        " marcadores para todos os tipos."
      )
    }

    cat(
      "top",
      k,
      ": ",
      nrow(
        top_k
      ),
      " atribuições; ",
      uniqueN(
        top_k$gene_id
      ),
      " GeneIDs únicos na união dos 8 tipos.\n",
      sep = ""
    )
  }


  # ============================================================================
  # 4.10 VALIDAR NESTING DOS CONJUNTOS
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "9. VALIDAÇÃO top100 ⊂ top300 ⊂ top500\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  nesting_rows <- list()

  counter <- 1L

  for (
    cell_type_i in
    CELL_TYPE_ORDER
  ) {

    g100 <- top_sets$top100[
      cell_type ==
        cell_type_i,
      gene_id
    ]

    g300 <- top_sets$top300[
      cell_type ==
        cell_type_i,
      gene_id
    ]

    g500 <- top_sets$top500[
      cell_type ==
        cell_type_i,
      gene_id
    ]

    nesting_rows[[counter]] <- data.table(
      cell_type = cell_type_i,
      top100_in_top300 = all(
        g100 %in%
          g300
      ),
      top300_in_top500 = all(
        g300 %in%
          g500
      )
    )

    counter <- counter + 1L
  }

  nesting_summary <- rbindlist(
    nesting_rows
  )

  fwrite(
    nesting_summary,
    file.path(
      result_dir,
      "04_marker_set_nesting.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(
    nesting_summary
  )

  if (
    any(
      !nesting_summary$top100_in_top300
    ) |
      any(
        !nesting_summary$top300_in_top500
      )
  ) {

    stop(
      "Os conjuntos de sensibilidade não são perfeitamente aninhados."
    )
  }


  # ============================================================================
  # 4.11 SALVAR LISTAS DE MARCADORES
  # ==============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "10. SALVAR LISTAS top100 / top300 / top500\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  fwrite(
    top_sets$top100,
    top100_file,
    sep = "\t",
    na = "NA"
  )

  fwrite(
    top_sets$top300,
    top300_file,
    sep = "\t",
    na = "NA"
  )

  fwrite(
    top_sets$top500,
    top500_file,
    sep = "\t",
    na = "NA"
  )


  # ============================================================================
  # 4.12 CONSTRUIR SIGNATURE MATRICES
  # ==============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "11. SIGNATURE MATRICES DE mean CPM\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  signature100 <- make_signature_table(
    top_sets$top100,
    mean_cpm_matrix,
    gene_ids
  )

  signature300 <- make_signature_table(
    top_sets$top300,
    mean_cpm_matrix,
    gene_ids
  )

  signature500 <- make_signature_table(
    top_sets$top500,
    mean_cpm_matrix,
    gene_ids
  )

  fwrite(
    signature100,
    signature100_file,
    sep = "\t",
    na = "NA"
  )

  fwrite(
    signature300,
    signature300_file,
    sep = "\t",
    na = "NA"
  )

  fwrite(
    signature500,
    signature500_file,
    sep = "\t",
    na = "NA"
  )

  cat(
    "Top100 signature: ",
    nrow(
      signature100
    ),
    " GeneIDs únicos\n",
    sep = ""
  )

  cat(
    "Top300 signature: ",
    nrow(
      signature300
    ),
    " GeneIDs únicos\n",
    sep = ""
  )

  cat(
    "Top500 signature: ",
    nrow(
      signature500
    ),
    " GeneIDs únicos\n",
    sep = ""
  )


  # ============================================================================
  # 4.13 AUDITAR OVERLAP ENTRE TIPOS CELULARES
  # ==============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "12. OVERLAP ENTRE LISTAS DE MARCADORES\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  overlap100 <- pairwise_marker_overlap(
    top_sets$top100,
    100L
  )

  overlap300 <- pairwise_marker_overlap(
    top_sets$top300,
    300L
  )

  overlap500 <- pairwise_marker_overlap(
    top_sets$top500,
    500L
  )

  overlap_all <- rbindlist(
    list(
      overlap100,
      overlap300,
      overlap500
    ),
    use.names = TRUE,
    fill = TRUE
  )

  fwrite(
    overlap_all,
    file.path(
      result_dir,
      "04_pairwise_marker_overlap.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  # Resumo somente dos pares diferentes para o conjunto principal.
  overlap300_between_types <- overlap300[
    cell_type_1 !=
      cell_type_2
  ]

  max_top300_overlap <- if (
    nrow(
      overlap300_between_types
    ) > 0L
  ) {
    max(
      overlap300_between_types$overlap_n
    )
  } else {
    NA_integer_
  }

  max_top300_jaccard <- if (
    nrow(
      overlap300_between_types
    ) > 0L
  ) {
    max(
      overlap300_between_types$jaccard
    )
  } else {
    NA_real_
  }

  cat(
    "Maior overlap entre dois tipos no top300: ",
    max_top300_overlap,
    " genes\n",
    sep = ""
  )

  cat(
    "Maior Jaccard entre dois tipos no top300: ",
    round(
      max_top300_jaccard,
      6
    ),
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.14 QC DE ROBUSTEZ DOS MARCADORES top300
  # ==============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "13. ROBUSTEZ DOS MARCADORES top300\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  top300_robustness_list <- vector(
    "list",
    length(
      CELL_TYPE_ORDER
    )
  )

  names(
    top300_robustness_list
  ) <- CELL_TYPE_ORDER

  for (
    cell_type_i in
    CELL_TYPE_ORDER
  ) {

    x <- top_sets$top300[
      cell_type ==
        cell_type_i
    ]

    top300_robustness_list[[cell_type_i]] <- data.table(
      cell_type = cell_type_i,
      n_markers = as.integer(
        nrow(
          x
        )
      ),
      n_reference_donors = as.integer(
        unique(
          x$n_donors
        )
      ),
      median_specificity = as.numeric(
        median(
          x$specificity
        )
      ),
      min_specificity = as.numeric(
        min(
          x$specificity
        )
      ),
      median_mean_CPM = as.numeric(
        median(
          x$mean_CPM
        )
      ),
      median_detection_fraction = as.numeric(
        median(
          x$detection_fraction
        )
      ),
      markers_detected_all_donors = as.integer(
        sum(
          x$n_donors_detected ==
            x$n_donors
        )
      ),
      markers_detected_at_least_half_donors = as.integer(
        sum(
          x$detection_fraction >= 0.5
        )
      ),
      markers_detected_less_than_half_donors = as.integer(
        sum(
          x$detection_fraction < 0.5
        )
      ),
      median_marker_score = as.numeric(
        median(
          x$marker_score
        )
      ),
      rare_reference_type = cell_type_i %in%
        c(
          "End",
          "Per"
        )
    )
  }

  top300_robustness <- rbindlist(
    top300_robustness_list,
    use.names = TRUE,
    fill = TRUE
  )

  fwrite(
    top300_robustness,
    file.path(
      result_dir,
      "04_top300_robustness_by_celltype.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(
    top300_robustness
  )


  # ============================================================================
  # 4.15 RESUMO DOS CONJUNTOS
  # ==============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "14. RESUMO DOS CONJUNTOS DE MARCADORES\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  marker_set_summary <- data.table(
    marker_set = c(
      "top100",
      "top300",
      "top500"
    ),
    markers_per_cell_type = c(
      100L,
      300L,
      500L
    ),
    total_assignments = c(
      nrow(
        top_sets$top100
      ),
      nrow(
        top_sets$top300
      ),
      nrow(
        top_sets$top500
      )
    ),
    unique_geneids_union = c(
      uniqueN(
        top_sets$top100$gene_id
      ),
      uniqueN(
        top_sets$top300$gene_id
      ),
      uniqueN(
        top_sets$top500$gene_id
      )
    ),
    role = c(
      "sensitivity",
      "primary",
      "sensitivity"
    )
  )

  fwrite(
    marker_set_summary,
    file.path(
      result_dir,
      "04_marker_set_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(
    marker_set_summary
  )


  # ============================================================================
  # 4.16 CHECKS FINAIS
  # ==============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "15. CHECKS FINAIS\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  top100_counts <- top_sets$top100[
    ,
    .N,
    by = cell_type
  ]$N

  top300_counts <- top_sets$top300[
    ,
    .N,
    by = cell_type
  ]$N

  top500_counts <- top_sets$top500[
    ,
    .N,
    by = cell_type
  ]$N

  checks <- rbindlist(
    list(

      make_check(
        "input_geneids",
        nrow(
          cpm
        ),
        EXPECTED_GENEIDS,
        if (
          nrow(
            cpm
          ) ==
            EXPECTED_GENEIDS
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "input_pseudobulks",
        ncol(
          cpm
        ),
        EXPECTED_PSEUDOBULKS,
        if (
          ncol(
            cpm
          ) ==
            EXPECTED_PSEUDOBULKS
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "eight_cell_types_present",
        setequal(
          observed_types,
          EXPECTED_CELL_TYPES
        ),
        TRUE,
        if (
          setequal(
            observed_types,
            EXPECTED_CELL_TYPES
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "specificity_sums_to_one",
        max_specificity_sum_deviation <=
          SPECIFICITY_SUM_TOLERANCE,
        TRUE,
        if (
          max_specificity_sum_deviation <=
            SPECIFICITY_SUM_TOLERANCE
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "all_celltypes_have_500_score_candidates",
        all(
          candidate_summary$candidate_genes >=
            max(MARKER_SET_SIZES)
        ),
        TRUE,
        if (
          all(
            candidate_summary$candidate_genes >=
              max(MARKER_SET_SIZES)
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        },
        "Candidates require mean_CPM > 0, specificity > 0, detection_fraction > 0 and marker_score > 0."
      ),

      make_check(
        "top300_marker_score_positive",
        all(
          top_sets$top300$marker_score > 0
        ),
        TRUE,
        if (
          all(
            top_sets$top300$marker_score > 0
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "top100_exact_per_celltype",
        all(
          top100_counts ==
            100L
        ),
        TRUE,
        if (
          all(
            top100_counts ==
              100L
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "top300_exact_per_celltype",
        all(
          top300_counts ==
            300L
        ),
        TRUE,
        if (
          all(
            top300_counts ==
              300L
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "top500_exact_per_celltype",
        all(
          top500_counts ==
            500L
        ),
        TRUE,
        if (
          all(
            top500_counts ==
              500L
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "marker_sets_nested",
        all(
          nesting_summary$top100_in_top300
        ) &
          all(
            nesting_summary$top300_in_top500
          ),
        TRUE,
        if (
          all(
            nesting_summary$top100_in_top300
          ) &
            all(
              nesting_summary$top300_in_top500
            )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "mean_matrix_has_no_negative_values",
        all(
          mean_cpm_matrix >=
            0
        ),
        TRUE,
        if (
          all(
            mean_cpm_matrix >=
              0
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "specificity_range_0_1",
        all(
          specificity_matrix >=
            -SPECIFICITY_SUM_TOLERANCE &
            specificity_matrix <=
              1 +
                SPECIFICITY_SUM_TOLERANCE
        ),
        TRUE,
        if (
          all(
            specificity_matrix >=
              -SPECIFICITY_SUM_TOLERANCE &
              specificity_matrix <=
                1 +
                  SPECIFICITY_SUM_TOLERANCE
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      )
    ),
    fill = TRUE
  )

  fwrite(
    checks,
    file.path(
      result_dir,
      "04_marker_derivation_status.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(
    checks
  )

  if (
    any(
      checks$status ==
        "FAIL",
      na.rm = TRUE
    )
  ) {

    stop(
      "Existe pelo menos um check FAIL após a derivação de marcadores."
    )
  }


  # ============================================================================
  # 4.17 MANIFESTO DOS OUTPUTS
  # ==============================================================================

  output_manifest <- data.table(
    object = c(
      "celltype_mean_cpm",
      "celltype_specificity",
      "marker_metrics_long",
      "marker_candidate_summary",
      "top100_markers",
      "top300_markers",
      "top500_markers",
      "signature_top100",
      "signature_top300",
      "signature_top500"
    ),
    path = c(
      mean_cpm_file,
      specificity_file,
      marker_metrics_file,
      file.path(
        result_dir,
        "04_marker_candidates_by_celltype.tsv"
      ),
      top100_file,
      top300_file,
      top500_file,
      signature100_file,
      signature300_file,
      signature500_file
    ),
    rows = c(
      nrow(
        mean_cpm_dt
      ),
      nrow(
        specificity_dt
      ),
      nrow(
        marker_metrics
      ),
      nrow(
        candidate_summary
      ),
      nrow(
        top_sets$top100
      ),
      nrow(
        top_sets$top300
      ),
      nrow(
        top_sets$top500
      ),
      nrow(
        signature100
      ),
      nrow(
        signature300
      ),
      nrow(
        signature500
      )
    ),
    columns = c(
      ncol(
        mean_cpm_dt
      ),
      ncol(
        specificity_dt
      ),
      ncol(
        marker_metrics
      ),
      ncol(
        candidate_summary
      ),
      ncol(
        top_sets$top100
      ),
      ncol(
        top_sets$top300
      ),
      ncol(
        top_sets$top500
      ),
      ncol(
        signature100
      ),
      ncol(
        signature300
      ),
      ncol(
        signature500
      )
    )
  )

  fwrite(
    output_manifest,
    file.path(
      result_dir,
      "04_output_manifest.tsv"
    ),
    sep = "\t",
    na = "NA"
  )


  # ============================================================================
  # 4.18 SESSION INFO
  # ==============================================================================

  sink(
    session_file
  )

  print(
    sessionInfo()
  )

  sink()


  # ============================================================================
  # 4.19 CONCLUSÃO
  # ==============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "16. CONCLUSÃO\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  cat(
    "GeneIDs avaliados: ",
    nrow(
      cpm
    ),
    "\n",
    sep = ""
  )

  cat(
    "Tipos celulares: ",
    length(
      CELL_TYPE_ORDER
    ),
    "\n",
    sep = ""
  )

  cat(
    "Pseudobulks: ",
    ncol(
      cpm
    ),
    "\n",
    sep = ""
  )

  cat(
    "Score de seleção: specificity * log1p(mean_CPM) * detection_fraction\n"
  )

  cat(
    "Conjunto principal: top",
    PRIMARY_MARKER_SET_SIZE,
    " por tipo celular pelo marker_score\n",
    sep = ""
  )

  cat(
    "GeneIDs únicos na união top300: ",
    uniqueN(
      top_sets$top300$gene_id
    ),
    "\n",
    sep = ""
  )

  cat(
    "Sensibilidades: top100 e top500 por tipo celular.\n"
  )

  cat(
    "End/Per permaneceram na referência e sua detecção entre doadores foi ",
    "mantida como métrica de robustez.\n"
  )

  cat(
    "\nPróxima etapa:\n"
  )

  cat(
    "  05_prepare_bulk_deconvolution.R\n"
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
    run_marker_derivation()
  },
  error = function(e) {

    execution_error <<- e

    cat(
      "\n============================================================\n"
    )

    cat(
      "ERRO NA DERIVAÇÃO DE MARCADORES\n"
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
  "04_derive_markers_Mathys2019.R CONCLUÍDO\n"
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
  "Resumo dos conjuntos:\n  ",
  file.path(
    result_dir,
    "04_marker_set_summary.tsv"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Robustez top300:\n  ",
  file.path(
    result_dir,
    "04_top300_robustness_by_celltype.tsv"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Status:\n  ",
  file.path(
    result_dir,
    "04_marker_derivation_status.tsv"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Marcadores principais top300:\n  ",
  top300_file,
  "\n\n",
  sep = ""
)

cat(
  "Signature matrix principal:\n  ",
  signature300_file,
  "\n",
  sep = ""
)
