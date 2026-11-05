#!/usr/bin/env Rscript

# ==============================================================================
# 06_run_MuSiC2.R
# ==============================================================================
#
# FASE 02 — DECONVOLUÇÃO CELULAR
# ETAPA 06 — MuSiC2
#
# Objetivo
# --------
# Estimar proporções celulares nos bulks AD / HD / PD usando MuSiC2 e a
# referência Mathys 2019 construída nas etapas 01–04.
#
# Datasets preparados no script 05:
#
#   ROSMAP       -> AD principal
#   GSE125583    -> AD validação independente
#   GSE64810     -> HD principal
#   GSE68719     -> PD principal
#
#
# IMPLEMENTAÇÃO DO MuSiC2
# -----------------------
# Este script usa a implementação atual de music2_prop() exposta pelo pacote
# MuSiC (repositório xuranw/MuSiC), cuja API espera:
#
#   bulk.control.mtx
#   bulk.case.mtx
#   sc.sce
#   clusters
#   samples
#   select.ct
#
# Essa API trabalha diretamente com uma matriz bulk e um SingleCellExperiment
# para a referência single-cell/single-nucleus.
#
# O script NÃO usa silenciosamente a API legada do pacote companheiro MuSiC2
# baseada em ExpressionSet. Antes de rodar, ele audita os argumentos formais da
# função instalada. Se a API não for a esperada, interrompe com mensagem clara.
#
#
# POR QUE raw counts?
# -------------------
# O tutorial e a implementação de MuSiC2 são construídos para matrizes de
# expressão não log-transformadas. Aqui usamos os raw counts preparados no 05.
#
# O script NÃO usa CPM nem log2(CPM+1) no MuSiC2.
#
#
# REFERÊNCIA
# ----------
# Referência principal:
#
#   24 doadores Mathys 2019 no-pathology
#
# A matriz nucleus-level harmonizada por GeneID é:
#
#   data/reference/mathys2019/phase02_prepared/reference/
#   mathys2019_pathology_free_geneid_counts_sparse.rds
#
# O metadata alinhado contém:
#
#   TAG
#   donor_id
#   cell_type
#   overlap_with_bulk_ROSMAP
#
#
# SENSIBILIDADE DE OVERLAP
# ------------------------
# O audit anterior encontrou 3 doadores Mathys compartilhados com o bulk ROSMAP.
#
# Por isso, além da referência principal de 24 doadores, este script pode rodar
# automaticamente uma análise de sensibilidade do ROSMAP removendo esses 3
# doadores:
#
#   main24                     -> 24 doadores
#   no_ROSMAP_overlap21        -> 21 doadores
#
# A referência de 21 doadores NÃO substitui a principal; ela serve apenas para
# verificar se as proporções do ROSMAP dependem materialmente do overlap.
#
#
# MARCADORES
# ----------
# Principal:
#
#   top300 por tipo celular
#
# O script usa a UNIÃO dos genes top300 realmente disponíveis naquele bulk.
#
# Importante:
#   - o bulk já foi intersectado com Mathys no script 05;
#   - genes ausentes NÃO são preenchidos com zero;
#   - a referência é subsetada para exatamente os GeneIDs presentes no bulk;
#   - genes all-zero no bulk ou na versão específica da referência são removidos
#     antes da chamada do MuSiC2 e auditados.
#
# O parâmetro markers=NULL é intencional:
#
#   como bulk e sc.sce já foram fisicamente restritos ao conjunto topK, todos os
#   genes disponíveis ao MuSiC2 pertencem ao conjunto de marcadores selecionado.
#
# Isso também permite que o algoritmo remova iterativamente genes com expressão
# diferencial célula-específica entre controle e caso sem tentar reintroduzir
# genes retirados da referência durante as iterações.
#
#
# CONFIGURAÇÃO PRINCIPAL
# ----------------------
# O padrão é rodar apenas top300, para não multiplicar o custo computacional:
#
#   MARKER_SETS_TO_RUN <- c("top300")
#
# Para a sensibilidade top100/top500, após a execução principal estar estável:
#
#   MARKER_SETS_TO_RUN <- c("top100", "top300", "top500")
#
#
# PARÂMETROS MuSiC2
# -----------------
# Mantemos explicitamente os defaults documentados, com maxiter = 200:
#
#   expr_low    = 20
#   prop_r      = 0.10
#   eps_c       = 0.05
#   eps_r       = 0.01
#   n_resample  = 20
#   sample_prop = 0.50
#   cutoff_expr = 0.05
#   cutoff_c    = 0.05
#   cutoff_r    = 0.01
#   maxiter     = 200
#   ct.cov      = FALSE
#   centered    = FALSE
#   normalize   = FALSE
#
#
# SAÍDAS
# ------
# results/deconvolution/MuSiC2/
#
# Para cada execução:
#
#   <dataset>/<reference_version>/<marker_set>/
#
#     proportions_wide.tsv
#     proportions_long.tsv
#     proportion_summary_by_group.tsv
#     reference_support.tsv
#     input_audit.tsv
#     convergence.tsv
#     removed_DE_genes.tsv
#     parameters.tsv
#     music2_result.rds
#
# Consolidados:
#
#   06_api_audit.tsv
#   06_run_summary.tsv
#   06_proportions_long_all.tsv
#   06_proportion_summary_all.tsv
#   06_status.tsv
#   06_output_manifest.tsv
#   06_run_MuSiC2_sessionInfo.txt
#
#
# INSTALAÇÃO
# ----------
# O script não instala pacotes automaticamente e não altera renv.
#
# Se MuSiC não estiver instalado, instale no ambiente R atual, por exemplo:
#
#   if (!requireNamespace("remotes", quietly = TRUE)) {
#     install.packages("remotes")
#   }
#
#   remotes::install_github(
#     "xuranw/MuSiC",
#     upgrade = "never"
#   )
#
# Depois execute novamente este script.
#
# ==============================================================================


# ==============================================================================
# 0. CONFIGURAÇÃO DO USUÁRIO
# ==============================================================================

# Principal: somente top300.
MARKER_SETS_TO_RUN <- c("top300")

# Sensibilidade específica do ROSMAP para os 3 doadores compartilhados.
RUN_ROSMAP_NO_OVERLAP_SENSITIVITY <- TRUE

# Parar imediatamente se uma chamada music2_prop() produzir erro.
STOP_ON_RUN_ERROR <- TRUE

# Seed base. Cada execução recebe um offset determinístico.
BASE_SEED <- 1234L


# ==============================================================================
# 1. PARÂMETROS MuSiC2
# ==============================================================================

MUSIC2_EXPR_LOW <- 20
MUSIC2_PROP_R <- 0.10
MUSIC2_EPS_C <- 0.05
MUSIC2_EPS_R <- 0.01
MUSIC2_N_RESAMPLE <- 20L
MUSIC2_SAMPLE_PROP <- 0.50
MUSIC2_CUTOFF_EXPR <- 0.05
MUSIC2_CUTOFF_C <- 0.05
MUSIC2_CUTOFF_R <- 0.01
MUSIC2_MAXITER <- 200L
MUSIC2_CT_COV <- FALSE
MUSIC2_CENTERED <- FALSE
MUSIC2_NORMALIZE <- FALSE

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

EXPECTED_MAIN_DONORS <- 24L
EXPECTED_NO_OVERLAP_DONORS <- 21L
EXPECTED_MAIN_NUCLEI <- 35140L

EXPECTED_MARKER_SETS <- c(
  "top100",
  "top300",
  "top500"
)

EXPECTED_MARKER_UNION <- c(
  top100 = 793L,
  top300 = 2242L,
  top500 = 3495L
)

PROP_NEGATIVE_TOLERANCE <- 1e-8
PROP_SUM_TOLERANCE <- 0.02


# ==============================================================================
# 2. PACOTES
# ==============================================================================

base_required <- c(
  "data.table",
  "Matrix",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "S4Vectors"
)

missing_base <- base_required[
  !vapply(
    base_required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_base) > 0L) {
  stop(
    "Pacotes necessários não instalados: ",
    paste(missing_base, collapse = ", "),
    "\nInstale apenas os pacotes ausentes no ambiente atual."
  )
}

if (!requireNamespace("MuSiC", quietly = TRUE)) {
  stop(
    paste0(
      "Pacote MuSiC não está instalado.\n\n",
      "Esta etapa usa a implementação atual de music2_prop() em xuranw/MuSiC.\n",
      "Instale no ambiente atual sem alterar o renv automaticamente, por exemplo:\n\n",
      "  if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes')\n",
      "  remotes::install_github('xuranw/MuSiC', upgrade = 'never')\n\n",
      "Depois execute novamente o script."
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(MuSiC)
})


# ==============================================================================
# 3. CAMINHOS
# ==============================================================================

reference_dir <- file.path(
  "data",
  "reference",
  "mathys2019",
  "phase02_prepared",
  "reference"
)

markers_dir <- file.path(
  "data",
  "reference",
  "mathys2019",
  "phase02_prepared",
  "markers"
)

bulk_root <- file.path(
  "data",
  "processed",
  "deconvolution"
)

result_root <- file.path(
  "results",
  "deconvolution",
  "MuSiC2"
)

dir.create(
  result_root,
  recursive = TRUE,
  showWarnings = FALSE
)

reference_counts_file <- file.path(
  reference_dir,
  "mathys2019_pathology_free_geneid_counts_sparse.rds"
)

reference_cell_metadata_file <- file.path(
  reference_dir,
  "mathys2019_pathology_free_cell_metadata_aligned.tsv"
)

step05_status_file <- file.path(
  "results",
  "deconvolution",
  "bulk_preparation",
  "05_status.tsv"
)

marker_files <- list(
  top100 = file.path(
    markers_dir,
    "mathys2019_top100_markers.tsv"
  ),
  top300 = file.path(
    markers_dir,
    "mathys2019_top300_markers.tsv"
  ),
  top500 = file.path(
    markers_dir,
    "mathys2019_top500_markers.tsv"
  )
)


# ==============================================================================
# 4. DATASETS
# ==============================================================================

dataset_registry <- data.table(
  dataset = c(
    "ROSMAP",
    "GSE125583",
    "GSE64810",
    "GSE68719"
  ),
  disease = c(
    "AD",
    "AD",
    "HD",
    "PD"
  ),
  cohort_role = c(
    "principal",
    "validation",
    "principal",
    "principal"
  )
)


# ==============================================================================
# 5. FUNÇÕES AUXILIARES
# ==============================================================================

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | toupper(x) == "NA"] <- NA_character_
  x
}


parse_logical <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  y <- toupper(clean_chr(x))

  out <- rep(
    NA,
    length(y)
  )

  out[
    y %in% c(
      "TRUE",
      "T",
      "1",
      "YES",
      "Y"
    )
  ] <- TRUE

  out[
    y %in% c(
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
  dataset,
  reference_version,
  marker_set,
  check,
  observed,
  expected,
  status,
  note = NA_character_
) {

  data.table(
    dataset = dataset,
    reference_version = reference_version,
    marker_set = marker_set,
    check = check,
    observed = as.character(observed),
    expected = as.character(expected),
    status = status,
    note = note
  )
}


validate_step05 <- function() {

  if (!file.exists(step05_status_file)) {
    stop(
      "Arquivo de status do script 05 não encontrado:\n",
      step05_status_file
    )
  }

  x <- fread(
    step05_status_file,
    na.strings = c(
      "NA",
      ""
    )
  )

  if (!"status" %in% names(x)) {
    stop(
      "05_status.tsv não contém coluna status."
    )
  }

  if (any(x$status == "FAIL", na.rm = TRUE)) {
    stop(
      "O script 05 possui pelo menos um check FAIL."
    )
  }

  invisible(TRUE)
}


validate_marker_sets <- function() {

  bad <- setdiff(
    MARKER_SETS_TO_RUN,
    EXPECTED_MARKER_SETS
  )

  if (length(bad) > 0L) {
    stop(
      "MARKER_SETS_TO_RUN contém conjunto inválido: ",
      paste(bad, collapse = ", ")
    )
  }

  invisible(TRUE)
}


safe_pkg_version <- function(pkg) {

  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(NA_character_)
  }

  as.character(
    utils::packageVersion(pkg)
  )
}


# ------------------------------------------------------------------------------
# Auditar API instalada do music2_prop.
# ------------------------------------------------------------------------------

audit_music2_api <- function() {

  if (
    !exists(
      "music2_prop",
      envir = asNamespace("MuSiC"),
      inherits = FALSE
    )
  ) {
    stop(
      paste0(
        "O pacote MuSiC instalado não contém music2_prop().\n",
        "A instalação parece incompatível com este pipeline.\n",
        "Use a implementação atual de xuranw/MuSiC."
      )
    )
  }

  fun <- get(
    "music2_prop",
    envir = asNamespace("MuSiC"),
    inherits = FALSE
  )

  actual_formals <- names(
    formals(fun)
  )

  required_formals <- c(
    "bulk.control.mtx",
    "bulk.case.mtx",
    "sc.sce",
    "clusters",
    "samples",
    "select.ct",
    "expr_low",
    "prop_r",
    "eps_c",
    "eps_r",
    "n_resample",
    "sample_prop",
    "cutoff_expr",
    "cutoff_c",
    "cutoff_r",
    "maxiter",
    "markers",
    "cell_size",
    "ct.cov",
    "centered",
    "normalize"
  )

  missing_formals <- setdiff(
    required_formals,
    actual_formals
  )

  if (length(missing_formals) > 0L) {
    stop(
      paste0(
        "A API instalada de MuSiC::music2_prop() não é a esperada.\n",
        "Argumentos ausentes: ",
        paste(missing_formals, collapse = ", "),
        "\n\nArgumentos encontrados:\n",
        paste(actual_formals, collapse = ", ")
      )
    )
  }

  api <- data.table(
    package = "MuSiC",
    package_version = safe_pkg_version("MuSiC"),
    function_name = "music2_prop",
    api_style = "bulk.control.mtx + bulk.case.mtx + SingleCellExperiment",
    formal_arguments = paste(
      actual_formals,
      collapse = ";"
    ),
    required_api_present = TRUE
  )

  fwrite(
    api,
    file.path(
      result_root,
      "06_api_audit.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  fun
}


# ------------------------------------------------------------------------------
# Validar matriz numérica de counts.
# ------------------------------------------------------------------------------

validate_counts <- function(
  x,
  label
) {

  if (
    is.null(
      rownames(x)
    ) ||
      is.null(
        colnames(x)
      )
  ) {
    stop(
      label,
      ": matriz sem rownames/colnames."
    )
  }

  if (anyDuplicated(rownames(x))) {
    stop(
      label,
      ": GeneID duplicado."
    )
  }

  if (anyDuplicated(colnames(x))) {
    stop(
      label,
      ": sample IDs duplicados."
    )
  }

  if (inherits(x, "sparseMatrix")) {
    values <- x@x
  } else {
    values <- as.numeric(x)
  }

  if (anyNA(values)) {
    stop(
      label,
      ": existem NA."
    )
  }

  if (any(!is.finite(values))) {
    stop(
      label,
      ": existem valores não finitos."
    )
  }

  if (any(values < 0)) {
    stop(
      label,
      ": existem valores negativos."
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# Validar que counts são inteiros.
# ------------------------------------------------------------------------------

validate_integer_counts <- function(
  x,
  label,
  tolerance = 1e-8
) {

  if (inherits(x, "sparseMatrix")) {
    values <- x@x
  } else {
    values <- as.numeric(x)
  }

  if (
    any(
      abs(
        values -
          round(values)
      ) >
        tolerance
    )
  ) {
    stop(
      label,
      ": valores fracionários detectados. MuSiC2 deve receber raw counts."
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# Carregar metadata bulk do script 05.
# ------------------------------------------------------------------------------

read_bulk_metadata <- function(dataset) {

  file <- file.path(
    bulk_root,
    dataset,
    "metadata_deconvolution.tsv"
  )

  if (!file.exists(file)) {
    stop(
      dataset,
      ": metadata do script 05 não encontrado:\n",
      file
    )
  }

  x <- fread(
    file,
    na.strings = c(
      "NA",
      ""
    )
  )

  required <- c(
    "sample_id",
    "group",
    "disease",
    "cohort_role"
  )

  missing <- setdiff(
    required,
    names(x)
  )

  if (length(missing) > 0L) {
    stop(
      dataset,
      ": metadata sem colunas: ",
      paste(missing, collapse = ", ")
    )
  }

  x[
    ,
    sample_id := clean_chr(sample_id)
  ]

  x[
    ,
    group := clean_chr(group)
  ]

  if (
    anyNA(x$sample_id) ||
      anyDuplicated(x$sample_id)
  ) {
    stop(
      dataset,
      ": sample_id inválido no metadata."
    )
  }

  x
}


# ------------------------------------------------------------------------------
# Carregar bulk topK.
# ------------------------------------------------------------------------------

read_bulk_marker_counts <- function(
  dataset,
  marker_set
) {

  file <- file.path(
    bulk_root,
    dataset,
    paste0(
      "bulk_counts_",
      marker_set,
      ".rds"
    )
  )

  if (!file.exists(file)) {
    stop(
      dataset,
      " / ",
      marker_set,
      ": bulk counts não encontrado:\n",
      file
    )
  }

  x <- readRDS(file)

  if (inherits(x, "sparseMatrix")) {
    x <- as(
      x,
      "dgCMatrix"
    )
  } else {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
  }

  validate_counts(
    x,
    paste(dataset, marker_set)
  )

  validate_integer_counts(
    x,
    paste(dataset, marker_set)
  )

  x
}


# ------------------------------------------------------------------------------
# Carregar união oficial de marcadores topK.
# ------------------------------------------------------------------------------

read_marker_union <- function(marker_set) {

  file <- marker_files[[marker_set]]

  if (!file.exists(file)) {
    stop(
      marker_set,
      ": arquivo de marcadores ausente:\n",
      file
    )
  }

  x <- fread(
    file,
    na.strings = c(
      "NA",
      ""
    )
  )

  required <- c(
    "gene_id",
    "cell_type",
    "marker_rank"
  )

  missing <- setdiff(
    required,
    names(x)
  )

  if (length(missing) > 0L) {
    stop(
      marker_set,
      ": colunas de marcador ausentes: ",
      paste(missing, collapse = ", ")
    )
  }

  x[
    ,
    gene_id := clean_chr(gene_id)
  ]

  genes <- unique(
    x$gene_id
  )

  expected <- EXPECTED_MARKER_UNION[
    marker_set
  ]

  if (
    length(genes) !=
      expected
  ) {
    warning(
      marker_set,
      ": união possui ",
      length(genes),
      " genes; esperado pelo QC anterior: ",
      expected,
      "."
    )
  }

  genes
}


# ------------------------------------------------------------------------------
# Alinhar metadata da referência às colunas da matriz.
# ------------------------------------------------------------------------------

align_reference_metadata <- function(
  counts,
  metadata
) {

  required <- c(
    "TAG",
    "donor_id",
    "cell_type",
    "overlap_with_bulk_ROSMAP"
  )

  missing <- setdiff(
    required,
    names(metadata)
  )

  if (length(missing) > 0L) {
    stop(
      "Reference metadata sem colunas: ",
      paste(missing, collapse = ", ")
    )
  }

  metadata[
    ,
    TAG := clean_chr(TAG)
  ]

  metadata[
    ,
    donor_id := clean_chr(donor_id)
  ]

  metadata[
    ,
    cell_type := clean_chr(cell_type)
  ]

  metadata[
    ,
    overlap_with_bulk_ROSMAP :=
      parse_logical(
        overlap_with_bulk_ROSMAP
      )
  ]

  if (
    anyNA(metadata$TAG) ||
      anyDuplicated(metadata$TAG)
  ) {
    stop(
      "TAG inválido/duplicado no reference metadata."
    )
  }

  idx <- match(
    colnames(counts),
    metadata$TAG
  )

  if (anyNA(idx)) {
    stop(
      "Existem colunas da referência sem TAG correspondente."
    )
  }

  metadata <- metadata[idx]

  if (
    !identical(
      metadata$TAG,
      colnames(counts)
    )
  ) {
    stop(
      "Falha ao alinhar reference metadata à matriz."
    )
  }

  metadata
}


# ------------------------------------------------------------------------------
# Construir uma versão da referência:
#
# main24:
#   todos os 24 no-pathology.
#
# no_ROSMAP_overlap21:
#   exclui os 3 doadores compartilhados com o bulk ROSMAP.
# ------------------------------------------------------------------------------

subset_reference_version <- function(
  counts,
  metadata,
  reference_version
) {

  if (reference_version == "main24") {

    keep <- rep(
      TRUE,
      nrow(metadata)
    )

    expected_donors <- EXPECTED_MAIN_DONORS

  } else if (
    reference_version ==
      "no_ROSMAP_overlap21"
  ) {

    overlap_donors <- unique(
      metadata[
        overlap_with_bulk_ROSMAP == TRUE,
        donor_id
      ]
    )

    if (
      length(overlap_donors) !=
        3L
    ) {
      stop(
        "Esperávamos 3 doadores Mathys com overlap ROSMAP; encontrados: ",
        length(overlap_donors),
        "."
      )
    }

    keep <- !metadata$donor_id %in%
      overlap_donors

    expected_donors <- EXPECTED_NO_OVERLAP_DONORS

  } else {

    stop(
      "reference_version desconhecida: ",
      reference_version
    )
  }

  meta_sub <- metadata[keep]
  counts_sub <- counts[
    ,
    keep,
    drop = FALSE
  ]

  n_donors <- uniqueN(
    meta_sub$donor_id
  )

  if (
    n_donors !=
      expected_donors
  ) {
    stop(
      reference_version,
      ": esperávamos ",
      expected_donors,
      " doadores; encontrados: ",
      n_donors,
      "."
    )
  }

  observed_types <- unique(
    meta_sub$cell_type
  )

  if (
    !setequal(
      observed_types,
      EXPECTED_CELL_TYPES
    )
  ) {
    stop(
      reference_version,
      ": os 8 tipos celulares não estão presentes."
    )
  }

  list(
    counts = counts_sub,
    metadata = meta_sub
  )
}


# ------------------------------------------------------------------------------
# Preparar inputs de uma execução.
# ------------------------------------------------------------------------------

prepare_run_inputs <- function(
  dataset,
  disease,
  marker_set,
  reference_version,
  reference_counts,
  reference_metadata
) {

  bulk <- read_bulk_marker_counts(
    dataset,
    marker_set
  )

  bulk_meta <- read_bulk_metadata(
    dataset
  )

  if (
    !identical(
      colnames(bulk),
      bulk_meta$sample_id
    )
  ) {

    idx <- match(
      colnames(bulk),
      bulk_meta$sample_id
    )

    if (anyNA(idx)) {
      stop(
        dataset,
        ": bulk e metadata não podem ser alinhados."
      )
    }

    bulk_meta <- bulk_meta[idx]

    if (
      !identical(
        colnames(bulk),
        bulk_meta$sample_id
      )
    ) {
      stop(
        dataset,
        ": ordem bulk/metadata permanece inconsistente."
      )
    }
  }

  if (
    !all(
      c(
        "Control",
        disease
      ) %in%
        unique(
          bulk_meta$group
        )
    )
  ) {
    stop(
      dataset,
      ": grupos Control e ",
      disease,
      " não estão ambos presentes."
    )
  }

  marker_union <- read_marker_union(
    marker_set
  )

  if (
    any(
      !rownames(bulk) %in%
        marker_union
    )
  ) {
    stop(
      dataset,
      " / ",
      marker_set,
      ": bulk topK contém GeneID fora da união oficial de marcadores."
    )
  }

  ref_version <- subset_reference_version(
    counts = reference_counts,
    metadata = reference_metadata,
    reference_version = reference_version
  )

  sc_counts <- ref_version$counts
  sc_meta <- ref_version$metadata

  common_genes <- rownames(bulk)[
    rownames(bulk) %in%
      rownames(sc_counts)
  ]

  if (length(common_genes) == 0L) {
    stop(
      dataset,
      " / ",
      marker_set,
      ": nenhuma interseção bulk x referência."
    )
  }

  bulk <- bulk[
    common_genes,
    ,
    drop = FALSE
  ]

  sc_counts <- sc_counts[
    common_genes,
    ,
    drop = FALSE
  ]

  # ---------------------------------------------------------------------------
  # Remover genes totalmente zerados no bulk inteiro ou na versão específica
  # da referência.
  # ---------------------------------------------------------------------------

  bulk_nonzero <- Matrix::rowSums(
    as(
      bulk,
      "dgCMatrix"
    )
  ) > 0

  sc_nonzero <- Matrix::rowSums(
    sc_counts
  ) > 0

  keep_gene <- bulk_nonzero &
    sc_nonzero

  removed_allzero_bulk <- common_genes[
    !bulk_nonzero
  ]

  removed_allzero_reference <- common_genes[
    !sc_nonzero
  ]

  genes_used <- common_genes[
    keep_gene
  ]

  if (
    length(genes_used) <
      100L
  ) {
    stop(
      dataset,
      " / ",
      marker_set,
      " / ",
      reference_version,
      ": menos de 100 genes permanecem após QC."
    )
  }

  bulk <- bulk[
    genes_used,
    ,
    drop = FALSE
  ]

  sc_counts <- sc_counts[
    genes_used,
    ,
    drop = FALSE
  ]

  # Converter bulk para matriz base: music2_prop documenta matrix para bulk.
  bulk_dense <- as.matrix(
    bulk
  )

  storage.mode(
    bulk_dense
  ) <- "double"

  # SingleCellExperiment mantém counts esparsos.
  sc_coldata <- S4Vectors::DataFrame(
    cellType = factor(
      sc_meta$cell_type,
      levels = EXPECTED_CELL_TYPES
    ),
    sampleID = factor(
      sc_meta$donor_id
    ),
    row.names = sc_meta$TAG
  )

  sc_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = sc_counts
    ),
    colData = sc_coldata
  )

  rownames(sc_sce) <- genes_used
  colnames(sc_sce) <- sc_meta$TAG

  control_ids <- bulk_meta[
    group == "Control",
    sample_id
  ]

  case_ids <- bulk_meta[
    group == disease,
    sample_id
  ]

  bulk_control <- bulk_dense[
    ,
    control_ids,
    drop = FALSE
  ]

  bulk_case <- bulk_dense[
    ,
    case_ids,
    drop = FALSE
  ]

  if (
    ncol(bulk_control) < 2L ||
      ncol(bulk_case) < 2L
  ) {
    stop(
      dataset,
      ": MuSiC2 requer múltiplas amostras nos dois grupos."
    )
  }

  reference_support <- sc_meta[
    ,
    .(
      n_nuclei = .N,
      n_donors = uniqueN(
        donor_id
      )
    ),
    by = cell_type
  ]

  reference_support[
    ,
    cell_type_order :=
      match(
        cell_type,
        EXPECTED_CELL_TYPES
      )
  ]

  setorder(
    reference_support,
    cell_type_order
  )

  reference_support[
    ,
    cell_type_order := NULL
  ]

  list(
    bulk_control = bulk_control,
    bulk_case = bulk_case,
    bulk_metadata = bulk_meta,
    sc_sce = sc_sce,
    reference_support = reference_support,
    genes_initial = common_genes,
    genes_used = genes_used,
    removed_allzero_bulk = unique(
      removed_allzero_bulk
    ),
    removed_allzero_reference = unique(
      removed_allzero_reference
    ),
    n_reference_donors = uniqueN(
      sc_meta$donor_id
    ),
    n_reference_nuclei = nrow(
      sc_meta
    )
  )
}


# ------------------------------------------------------------------------------
# Normalizar orientação da matriz de proporções SEM renormalizar os valores.
#
# Saída:
#   linhas = samples
#   colunas = cell types
# ------------------------------------------------------------------------------

standardize_proportions <- function(
  prop,
  expected_samples
) {

  if (
    is.data.frame(prop)
  ) {
    prop <- as.matrix(prop)
  }

  if (
    !is.matrix(prop)
  ) {
    stop(
      "Est.prop não é matrix/data.frame."
    )
  }

  storage.mode(prop) <- "double"

  row_match <- !is.null(
    rownames(prop)
  ) &&
    all(
      expected_samples %in%
        rownames(prop)
    )

  col_match <- !is.null(
    colnames(prop)
  ) &&
    all(
      expected_samples %in%
        colnames(prop)
    )

  if (row_match) {

    out <- prop[
      expected_samples,
      ,
      drop = FALSE
    ]

  } else if (col_match) {

    out <- t(
      prop[
        ,
        expected_samples,
        drop = FALSE
      ]
    )

  } else {

    stop(
      "Não foi possível identificar a dimensão de amostras em Est.prop."
    )
  }

  missing_types <- setdiff(
    EXPECTED_CELL_TYPES,
    colnames(out)
  )

  if (length(missing_types) > 0L) {
    stop(
      "Est.prop sem tipos celulares: ",
      paste(
        missing_types,
        collapse = ", "
      )
    )
  }

  out <- out[
    ,
    EXPECTED_CELL_TYPES,
    drop = FALSE
  ]

  if (
    anyNA(out) ||
      any(
        !is.finite(out)
      )
  ) {
    stop(
      "Est.prop contém NA ou valores não finitos."
    )
  }

  if (
    any(
      out <
        -PROP_NEGATIVE_TOLERANCE
    )
  ) {
    stop(
      "Est.prop contém proporções negativas abaixo da tolerância."
    )
  }

  out
}


# ------------------------------------------------------------------------------
# Converter proporções para long e anexar metadata.
# ------------------------------------------------------------------------------

make_long_proportions <- function(
  prop,
  metadata,
  dataset,
  disease,
  cohort_role,
  reference_version,
  marker_set
) {

  x <- as.data.table(
    prop,
    keep.rownames = "sample_id"
  )

  long <- melt(
    x,
    id.vars = "sample_id",
    variable.name = "cell_type",
    value.name = "proportion"
  )

  long[
    ,
    cell_type := as.character(cell_type)
  ]

  long <- merge(
    long,
    metadata[
      ,
      .(
        sample_id,
        group
      )
    ],
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  if (anyNA(long$group)) {
    stop(
      dataset,
      ": group ausente após merge das proporções."
    )
  }

  long[
    ,
    `:=`(
      dataset = dataset,
      disease = disease,
      cohort_role = cohort_role,
      reference_version = reference_version,
      marker_set = marker_set
    )
  ]

  setcolorder(
    long,
    c(
      "dataset",
      "disease",
      "cohort_role",
      "reference_version",
      "marker_set",
      "sample_id",
      "group",
      "cell_type",
      "proportion"
    )
  )

  long
}


# ------------------------------------------------------------------------------
# Resumo de proporções por grupo.
# ------------------------------------------------------------------------------

summarize_proportions <- function(long_dt) {

  long_dt[
    ,
    .(
      n_samples = uniqueN(
        sample_id
      ),
      mean_proportion = mean(
        proportion
      ),
      sd_proportion = stats::sd(
        proportion
      ),
      median_proportion = stats::median(
        proportion
      ),
      min_proportion = min(
        proportion
      ),
      max_proportion = max(
        proportion
      )
    ),
    by = .(
      dataset,
      disease,
      cohort_role,
      reference_version,
      marker_set,
      group,
      cell_type
    )
  ]
}


# ------------------------------------------------------------------------------
# Extrair genes removidos pelo MuSiC2.
# ------------------------------------------------------------------------------

extract_removed_genes <- function(result) {

  if (
    is.null(
      result$DE.genes
    )
  ) {
    return(
      data.table(
        gene_id = character()
      )
    )
  }

  genes <- unique(
    clean_chr(
      unlist(
        result$DE.genes,
        use.names = FALSE
      )
    )
  )

  genes <- genes[
    !is.na(genes)
  ]

  data.table(
    gene_id = genes
  )
}


# ------------------------------------------------------------------------------
# Interpretar convergência.
# ------------------------------------------------------------------------------

interpret_convergence <- function(result) {

  conv_raw <- result$convergence

  if (is.null(conv_raw)) {

    converged <- NA

  } else {

    conv_logical <- as.logical(
      conv_raw
    )

    converged <- if (
      all(
        is.na(
          conv_logical
        )
      )
    ) {
      NA
    } else {
      all(
        conv_logical,
        na.rm = TRUE
      )
    }
  }

  n_iter <- if (
    is.null(
      result$n.iter
    )
  ) {
    NA_integer_
  } else {
    suppressWarnings(
      as.integer(
        result$n.iter[1]
      )
    )
  }

  not_converged_ids <- if (
    is.null(
      result$id.not.converge
    )
  ) {
    character()
  } else {
    unique(
      clean_chr(
        result$id.not.converge
      )
    )
  }

  not_converged_ids <- not_converged_ids[
    !is.na(
      not_converged_ids
    )
  ]

  list(
    converged = converged,
    n_iter = n_iter,
    id_not_converged = not_converged_ids
  )
}


# ------------------------------------------------------------------------------
# Salvar parâmetros de uma execução.
# ------------------------------------------------------------------------------

make_parameter_table <- function(
  dataset,
  disease,
  reference_version,
  marker_set,
  seed
) {

  data.table(
    parameter = c(
      "dataset",
      "disease",
      "reference_version",
      "marker_set",
      "seed",
      "expr_low",
      "prop_r",
      "eps_c",
      "eps_r",
      "n_resample",
      "sample_prop",
      "cutoff_expr",
      "cutoff_c",
      "cutoff_r",
      "maxiter",
      "markers_argument",
      "ct.cov",
      "centered",
      "normalize",
      "selected_cell_types"
    ),
    value = c(
      dataset,
      disease,
      reference_version,
      marker_set,
      as.character(seed),
      as.character(MUSIC2_EXPR_LOW),
      as.character(MUSIC2_PROP_R),
      as.character(MUSIC2_EPS_C),
      as.character(MUSIC2_EPS_R),
      as.character(MUSIC2_N_RESAMPLE),
      as.character(MUSIC2_SAMPLE_PROP),
      as.character(MUSIC2_CUTOFF_EXPR),
      as.character(MUSIC2_CUTOFF_C),
      as.character(MUSIC2_CUTOFF_R),
      as.character(MUSIC2_MAXITER),
      "NULL because matrices are already restricted to topK",
      as.character(MUSIC2_CT_COV),
      as.character(MUSIC2_CENTERED),
      as.character(MUSIC2_NORMALIZE),
      paste(
        EXPECTED_CELL_TYPES,
        collapse = ";"
      )
    )
  )
}


# ==============================================================================
# 6. PRÉ-VALIDAÇÃO
# ==============================================================================

cat(
  "============================================================\n"
)

cat(
  "FASE 02 — 06_run_MuSiC2.R\n"
)

cat(
  "============================================================\n\n"
)

validate_step05()
validate_marker_sets()

music2_fun <- audit_music2_api()

cat(
  "MuSiC version: ",
  safe_pkg_version("MuSiC"),
  "\n",
  sep = ""
)

cat(
  "API music2_prop validada.\n"
)

cat(
  "Marker sets: ",
  paste(
    MARKER_SETS_TO_RUN,
    collapse = ", "
  ),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 7. CARREGAR REFERÊNCIA NUCLEUS-LEVEL
# ==============================================================================

if (!file.exists(reference_counts_file)) {
  stop(
    "Reference counts não encontrado:\n",
    reference_counts_file
  )
}

if (!file.exists(reference_cell_metadata_file)) {
  stop(
    "Reference cell metadata não encontrado:\n",
    reference_cell_metadata_file
  )
}

reference_counts <- readRDS(
  reference_counts_file
)

reference_counts <- as(
  reference_counts,
  "dgCMatrix"
)

validate_counts(
  reference_counts,
  "Mathys nucleus-level"
)

validate_integer_counts(
  reference_counts,
  "Mathys nucleus-level"
)

reference_metadata <- fread(
  reference_cell_metadata_file,
  na.strings = c(
    "NA",
    ""
  )
)

reference_metadata <- align_reference_metadata(
  counts = reference_counts,
  metadata = reference_metadata
)

if (
  ncol(reference_counts) !=
    EXPECTED_MAIN_NUCLEI
) {
  stop(
    "Esperávamos ",
    EXPECTED_MAIN_NUCLEI,
    " núcleos Mathys; encontrados: ",
    ncol(reference_counts),
    "."
  )
}

if (
  uniqueN(
    reference_metadata$donor_id
  ) !=
    EXPECTED_MAIN_DONORS
) {
  stop(
    "Esperávamos ",
    EXPECTED_MAIN_DONORS,
    " doadores Mathys."
  )
}

if (
  !setequal(
    unique(
      reference_metadata$cell_type
    ),
    EXPECTED_CELL_TYPES
  )
) {
  stop(
    "Os 8 tipos celulares Mathys não estão presentes."
  )
}

cat(
  "Referência Mathys: ",
  nrow(reference_counts),
  " GeneIDs x ",
  ncol(reference_counts),
  " núcleos; ",
  uniqueN(reference_metadata$donor_id),
  " doadores.\n\n",
  sep = ""
)


# ==============================================================================
# 8. PLANO DE EXECUÇÃO
# ==============================================================================

run_plan <- CJ(
  dataset = dataset_registry$dataset,
  marker_set = MARKER_SETS_TO_RUN,
  unique = TRUE
)

run_plan <- merge(
  run_plan,
  dataset_registry,
  by = "dataset",
  all.x = TRUE,
  sort = FALSE
)

run_plan[
  ,
  reference_version := "main24"
]

if (RUN_ROSMAP_NO_OVERLAP_SENSITIVITY) {

  rosmap_sensitivity <- run_plan[
    dataset == "ROSMAP"
  ]

  rosmap_sensitivity[
    ,
    reference_version :=
      "no_ROSMAP_overlap21"
  ]

  run_plan <- rbind(
    run_plan,
    rosmap_sensitivity,
    fill = TRUE
  )
}

run_plan[
  ,
  run_order := seq_len(.N)
]

setorder(
  run_plan,
  run_order
)

fwrite(
  run_plan,
  file.path(
    result_root,
    "06_run_plan.tsv"
  ),
  sep = "\t",
  na = "NA"
)

cat(
  "Execuções planejadas: ",
  nrow(run_plan),
  "\n",
  sep = ""
)

print(
  run_plan[
    ,
    .(
      run_order,
      dataset,
      disease,
      cohort_role,
      reference_version,
      marker_set
    )
  ]
)


# ==============================================================================
# 9. EXECUTAR MuSiC2
# ==============================================================================

run_summaries <- list()
all_long <- list()
all_group_summary <- list()
all_status <- list()
all_manifest <- list()

summary_i <- 1L
long_i <- 1L
group_i <- 1L
status_i <- 1L
manifest_i <- 1L


for (
  plan_i in seq_len(
    nrow(run_plan)
  )
) {

  row <- run_plan[
    plan_i
  ]

  dataset_i <- row$dataset
  disease_i <- row$disease
  cohort_role_i <- row$cohort_role
  reference_version_i <- row$reference_version
  marker_set_i <- row$marker_set
  seed_i <- BASE_SEED + row$run_order - 1L

  cat(
    "\n============================================================\n"
  )

  cat(
    "RUN ",
    row$run_order,
    " / ",
    nrow(run_plan),
    ": ",
    dataset_i,
    " | ",
    reference_version_i,
    " | ",
    marker_set_i,
    "\n",
    sep = ""
  )

  cat(
    "============================================================\n"
  )

  run_dir <- file.path(
    result_root,
    dataset_i,
    reference_version_i,
    marker_set_i
  )

  dir.create(
    run_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # ---------------------------------------------------------------------------
  # Preparar inputs.
  # ---------------------------------------------------------------------------

  inputs <- prepare_run_inputs(
    dataset = dataset_i,
    disease = disease_i,
    marker_set = marker_set_i,
    reference_version = reference_version_i,
    reference_counts = reference_counts,
    reference_metadata = reference_metadata
  )

  expected_samples <- inputs$bulk_metadata$sample_id

  input_audit <- data.table(
    metric = c(
      "bulk_control_samples",
      "bulk_case_samples",
      "reference_donors",
      "reference_nuclei",
      "genes_before_allzero_QC",
      "genes_used",
      "genes_removed_allzero_bulk",
      "genes_removed_allzero_reference"
    ),
    value = c(
      ncol(
        inputs$bulk_control
      ),
      ncol(
        inputs$bulk_case
      ),
      inputs$n_reference_donors,
      inputs$n_reference_nuclei,
      length(
        inputs$genes_initial
      ),
      length(
        inputs$genes_used
      ),
      length(
        inputs$removed_allzero_bulk
      ),
      length(
        inputs$removed_allzero_reference
      )
    )
  )

  fwrite(
    input_audit,
    file.path(
      run_dir,
      "input_audit.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  fwrite(
    inputs$reference_support,
    file.path(
      run_dir,
      "reference_support.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  fwrite(
    data.table(
      gene_id = inputs$removed_allzero_bulk
    ),
    file.path(
      run_dir,
      "genes_removed_allzero_bulk.tsv"
    ),
    sep = "\t"
  )

  fwrite(
    data.table(
      gene_id = inputs$removed_allzero_reference
    ),
    file.path(
      run_dir,
      "genes_removed_allzero_reference.tsv"
    ),
    sep = "\t"
  )

  parameters <- make_parameter_table(
    dataset = dataset_i,
    disease = disease_i,
    reference_version = reference_version_i,
    marker_set = marker_set_i,
    seed = seed_i
  )

  fwrite(
    parameters,
    file.path(
      run_dir,
      "parameters.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  cat(
    "Controls: ",
    ncol(inputs$bulk_control),
    " | Cases: ",
    ncol(inputs$bulk_case),
    "\n",
    sep = ""
  )

  cat(
    "Genes usados: ",
    length(inputs$genes_used),
    "\n",
    sep = ""
  )

  cat(
    "Reference: ",
    inputs$n_reference_donors,
    " donors / ",
    inputs$n_reference_nuclei,
    " nuclei\n",
    sep = ""
  )

  # ---------------------------------------------------------------------------
  # Rodar MuSiC2.
  # ---------------------------------------------------------------------------

  set.seed(
    seed_i
  )

  start_time <- Sys.time()

  run_error <- NULL

  result <- tryCatch(
    {
      music2_fun(
        bulk.control.mtx = inputs$bulk_control,
        bulk.case.mtx = inputs$bulk_case,
        sc.sce = inputs$sc_sce,
        clusters = "cellType",
        samples = "sampleID",
        select.ct = EXPECTED_CELL_TYPES,
        expr_low = MUSIC2_EXPR_LOW,
        prop_r = MUSIC2_PROP_R,
        eps_c = MUSIC2_EPS_C,
        eps_r = MUSIC2_EPS_R,
        n_resample = MUSIC2_N_RESAMPLE,
        sample_prop = MUSIC2_SAMPLE_PROP,
        cutoff_expr = MUSIC2_CUTOFF_EXPR,
        cutoff_c = MUSIC2_CUTOFF_C,
        cutoff_r = MUSIC2_CUTOFF_R,
        maxiter = MUSIC2_MAXITER,
        markers = NULL,
        cell_size = NULL,
        ct.cov = MUSIC2_CT_COV,
        centered = MUSIC2_CENTERED,
        normalize = MUSIC2_NORMALIZE
      )
    },
    error = function(e) {
      run_error <<- conditionMessage(e)
      NULL
    }
  )

  end_time <- Sys.time()

  elapsed_seconds <- as.numeric(
    difftime(
      end_time,
      start_time,
      units = "secs"
    )
  )

  if (!is.null(run_error)) {

    error_dt <- data.table(
      dataset = dataset_i,
      reference_version = reference_version_i,
      marker_set = marker_set_i,
      error = run_error
    )

    fwrite(
      error_dt,
      file.path(
        run_dir,
        "ERROR.tsv"
      ),
      sep = "\t",
      na = "NA"
    )

    all_status[[status_i]] <- make_check(
      dataset = dataset_i,
      reference_version = reference_version_i,
      marker_set = marker_set_i,
      check = "music2_execution",
      observed = run_error,
      expected = "successful execution",
      status = "FAIL",
      note = "Veja ERROR.tsv no diretório da execução."
    )

    status_i <- status_i + 1L

    if (STOP_ON_RUN_ERROR) {
      stop(
        dataset_i,
        " / ",
        reference_version_i,
        " / ",
        marker_set_i,
        ": MuSiC2 falhou:\n",
        run_error
      )
    } else {
      next
    }
  }

  # ---------------------------------------------------------------------------
  # Salvar objeto bruto.
  # ---------------------------------------------------------------------------

  result_rds <- file.path(
    run_dir,
    "music2_result.rds"
  )

  saveRDS(
    result,
    result_rds,
    compress = FALSE
  )

  # ---------------------------------------------------------------------------
  # Proporções.
  # ---------------------------------------------------------------------------

  if (
    is.null(
      result$Est.prop
    )
  ) {
    stop(
      dataset_i,
      ": music2_prop retornou objeto sem Est.prop."
    )
  }

  prop <- standardize_proportions(
    result$Est.prop,
    expected_samples = expected_samples
  )

  row_sums <- rowSums(
    prop
  )

  prop_wide <- as.data.table(
    prop,
    keep.rownames = "sample_id"
  )

  prop_wide <- merge(
    inputs$bulk_metadata[
      ,
      .(
        sample_id,
        group
      )
    ],
    prop_wide,
    by = "sample_id",
    all.y = TRUE,
    sort = FALSE
  )

  prop_wide <- prop_wide[
    match(
      expected_samples,
      sample_id
    )
  ]

  fwrite(
    prop_wide,
    file.path(
      run_dir,
      "proportions_wide.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  prop_long <- make_long_proportions(
    prop = prop,
    metadata = inputs$bulk_metadata,
    dataset = dataset_i,
    disease = disease_i,
    cohort_role = cohort_role_i,
    reference_version = reference_version_i,
    marker_set = marker_set_i
  )

  group_summary <- summarize_proportions(
    prop_long
  )

  fwrite(
    prop_long,
    file.path(
      run_dir,
      "proportions_long.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  fwrite(
    group_summary,
    file.path(
      run_dir,
      "proportion_summary_by_group.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  # ---------------------------------------------------------------------------
  # Convergência.
  # ---------------------------------------------------------------------------

  convergence_info <- interpret_convergence(
    result
  )

  convergence_dt <- data.table(
    dataset = dataset_i,
    disease = disease_i,
    cohort_role = cohort_role_i,
    reference_version = reference_version_i,
    marker_set = marker_set_i,
    converged = convergence_info$converged,
    n_iter = convergence_info$n_iter,
    n_id_not_converged = length(
      convergence_info$id_not_converged
    ),
    id_not_converged = if (
      length(
        convergence_info$id_not_converged
      ) == 0L
    ) {
      NA_character_
    } else {
      paste(
        convergence_info$id_not_converged,
        collapse = ";"
      )
    }
  )

  fwrite(
    convergence_dt,
    file.path(
      run_dir,
      "convergence.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  # ---------------------------------------------------------------------------
  # Genes removidos iterativamente.
  # ---------------------------------------------------------------------------

  removed_genes <- extract_removed_genes(
    result
  )

  fwrite(
    removed_genes,
    file.path(
      run_dir,
      "removed_DE_genes.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  # ---------------------------------------------------------------------------
  # Resumo da execução.
  # ---------------------------------------------------------------------------

  run_summary <- data.table(
    dataset = dataset_i,
    disease = disease_i,
    cohort_role = cohort_role_i,
    reference_version = reference_version_i,
    marker_set = marker_set_i,
    seed = seed_i,
    controls = ncol(
      inputs$bulk_control
    ),
    cases = ncol(
      inputs$bulk_case
    ),
    samples_total = length(
      expected_samples
    ),
    reference_donors = inputs$n_reference_donors,
    reference_nuclei = inputs$n_reference_nuclei,
    genes_before_allzero_QC = length(
      inputs$genes_initial
    ),
    genes_used = length(
      inputs$genes_used
    ),
    genes_removed_allzero_bulk = length(
      inputs$removed_allzero_bulk
    ),
    genes_removed_allzero_reference = length(
      inputs$removed_allzero_reference
    ),
    music2_DE_genes_removed = nrow(
      removed_genes
    ),
    converged = convergence_info$converged,
    n_iter = convergence_info$n_iter,
    n_id_not_converged = length(
      convergence_info$id_not_converged
    ),
    proportion_row_sum_min = min(
      row_sums
    ),
    proportion_row_sum_median = median(
      row_sums
    ),
    proportion_row_sum_max = max(
      row_sums
    ),
    proportion_min = min(
      prop
    ),
    proportion_max = max(
      prop
    ),
    elapsed_seconds = elapsed_seconds,
    MuSiC_version = safe_pkg_version(
      "MuSiC"
    )
  )

  fwrite(
    run_summary,
    file.path(
      run_dir,
      "run_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  run_summaries[[summary_i]] <- run_summary
  summary_i <- summary_i + 1L

  all_long[[long_i]] <- prop_long
  long_i <- long_i + 1L

  all_group_summary[[group_i]] <- group_summary
  group_i <- group_i + 1L

  # ---------------------------------------------------------------------------
  # Status.
  # ---------------------------------------------------------------------------

  convergence_status <- if (
    isTRUE(
      convergence_info$converged
    ) &&
      length(
        convergence_info$id_not_converged
      ) == 0L
  ) {
    "PASS"
  } else {
    "WARN"
  }

  convergence_note <- if (
    convergence_status ==
      "PASS"
  ) {
    NA_character_
  } else {
    paste0(
      "MuSiC2 não reportou convergência completa. ",
      "Preservar estimativas, mas revisar convergence.tsv."
    )
  }

  all_status[[status_i]] <- make_check(
    dataset = dataset_i,
    reference_version = reference_version_i,
    marker_set = marker_set_i,
    check = "music2_execution",
    observed = "completed",
    expected = "completed",
    status = "PASS"
  )
  status_i <- status_i + 1L

  all_status[[status_i]] <- make_check(
    dataset = dataset_i,
    reference_version = reference_version_i,
    marker_set = marker_set_i,
    check = "all_samples_returned",
    observed = nrow(prop),
    expected = length(expected_samples),
    status = if (
      nrow(prop) ==
        length(expected_samples)
    ) {
      "PASS"
    } else {
      "FAIL"
    }
  )
  status_i <- status_i + 1L

  all_status[[status_i]] <- make_check(
    dataset = dataset_i,
    reference_version = reference_version_i,
    marker_set = marker_set_i,
    check = "eight_cell_types_returned",
    observed = ncol(prop),
    expected = length(EXPECTED_CELL_TYPES),
    status = if (
      ncol(prop) ==
        length(EXPECTED_CELL_TYPES)
    ) {
      "PASS"
    } else {
      "FAIL"
    }
  )
  status_i <- status_i + 1L

  all_status[[status_i]] <- make_check(
    dataset = dataset_i,
    reference_version = reference_version_i,
    marker_set = marker_set_i,
    check = "proportions_nonnegative",
    observed = min(prop),
    expected = paste0(
      ">=",
      -PROP_NEGATIVE_TOLERANCE
    ),
    status = if (
      min(prop) >=
        -PROP_NEGATIVE_TOLERANCE
    ) {
      "PASS"
    } else {
      "FAIL"
    }
  )
  status_i <- status_i + 1L

  max_sum_deviation <- max(
    abs(
      row_sums -
        1
    )
  )

  all_status[[status_i]] <- make_check(
    dataset = dataset_i,
    reference_version = reference_version_i,
    marker_set = marker_set_i,
    check = "proportion_rows_sum_to_one",
    observed = max_sum_deviation,
    expected = paste0(
      "<=",
      PROP_SUM_TOLERANCE
    ),
    status = if (
      max_sum_deviation <=
        PROP_SUM_TOLERANCE
    ) {
      "PASS"
    } else {
      "WARN"
    },
    note = if (
      max_sum_deviation <=
        PROP_SUM_TOLERANCE
    ) {
      NA_character_
    } else {
      "Os valores foram preservados exatamente como retornados; não houve renormalização pós-hoc."
    }
  )
  status_i <- status_i + 1L

  all_status[[status_i]] <- make_check(
    dataset = dataset_i,
    reference_version = reference_version_i,
    marker_set = marker_set_i,
    check = "music2_convergence",
    observed = convergence_info$converged,
    expected = TRUE,
    status = convergence_status,
    note = convergence_note
  )
  status_i <- status_i + 1L

  # ---------------------------------------------------------------------------
  # Manifesto.
  # ---------------------------------------------------------------------------

  output_files <- c(
    file.path(
      run_dir,
      "proportions_wide.tsv"
    ),
    file.path(
      run_dir,
      "proportions_long.tsv"
    ),
    file.path(
      run_dir,
      "proportion_summary_by_group.tsv"
    ),
    file.path(
      run_dir,
      "reference_support.tsv"
    ),
    file.path(
      run_dir,
      "input_audit.tsv"
    ),
    file.path(
      run_dir,
      "convergence.tsv"
    ),
    file.path(
      run_dir,
      "removed_DE_genes.tsv"
    ),
    file.path(
      run_dir,
      "parameters.tsv"
    ),
    file.path(
      run_dir,
      "run_summary.tsv"
    ),
    result_rds
  )

  all_manifest[[manifest_i]] <- data.table(
    dataset = dataset_i,
    reference_version = reference_version_i,
    marker_set = marker_set_i,
    path = output_files,
    exists = file.exists(
      output_files
    )
  )

  manifest_i <- manifest_i + 1L

  cat(
    "Convergência: ",
    as.character(
      convergence_info$converged
    ),
    " | iterações: ",
    as.character(
      convergence_info$n_iter
    ),
    " | genes DE removidos: ",
    nrow(
      removed_genes
    ),
    " | tempo: ",
    round(
      elapsed_seconds,
      1
    ),
    " s\n",
    sep = ""
  )

  # Liberar objetos grandes da execução atual.
  rm(
    inputs,
    result,
    prop,
    prop_wide,
    prop_long,
    group_summary,
    removed_genes
  )

  invisible(
    gc()
  )
}


# ==============================================================================
# 10. CONSOLIDAR RESULTADOS
# ==============================================================================

run_summary_all <- rbindlist(
  run_summaries,
  use.names = TRUE,
  fill = TRUE
)

proportions_long_all <- rbindlist(
  all_long,
  use.names = TRUE,
  fill = TRUE
)

proportion_summary_all <- rbindlist(
  all_group_summary,
  use.names = TRUE,
  fill = TRUE
)

status_all <- rbindlist(
  all_status,
  use.names = TRUE,
  fill = TRUE
)

manifest_all <- rbindlist(
  all_manifest,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  run_summary_all,
  file.path(
    result_root,
    "06_run_summary.tsv"
  ),
  sep = "\t",
  na = "NA"
)

fwrite(
  proportions_long_all,
  file.path(
    result_root,
    "06_proportions_long_all.tsv"
  ),
  sep = "\t",
  na = "NA"
)

fwrite(
  proportion_summary_all,
  file.path(
    result_root,
    "06_proportion_summary_all.tsv"
  ),
  sep = "\t",
  na = "NA"
)

fwrite(
  status_all,
  file.path(
    result_root,
    "06_status.tsv"
  ),
  sep = "\t",
  na = "NA"
)

fwrite(
  manifest_all,
  file.path(
    result_root,
    "06_output_manifest.tsv"
  ),
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 11. COMPARAÇÃO ROSMAP main24 vs no-overlap21
# ==============================================================================

rosmap_sensitivity_file <- file.path(
  result_root,
  "06_ROSMAP_reference_overlap_sensitivity.tsv"
)

if (
  RUN_ROSMAP_NO_OVERLAP_SENSITIVITY &&
    "top300" %in%
      MARKER_SETS_TO_RUN
) {

  x_main <- proportions_long_all[
    dataset == "ROSMAP" &
      reference_version == "main24" &
      marker_set == "top300"
  ]

  x_sens <- proportions_long_all[
    dataset == "ROSMAP" &
      reference_version == "no_ROSMAP_overlap21" &
      marker_set == "top300"
  ]

  if (
    nrow(x_main) > 0L &&
      nrow(x_sens) > 0L
  ) {

    comp <- merge(
      x_main[
        ,
        .(
          sample_id,
          cell_type,
          proportion_main24 = proportion
        )
      ],
      x_sens[
        ,
        .(
          sample_id,
          cell_type,
          proportion_no_overlap21 = proportion
        )
      ],
      by = c(
        "sample_id",
        "cell_type"
      ),
      all = FALSE
    )

    rosmap_sensitivity <- comp[
      ,
      .(
        n_samples = .N,
        spearman = suppressWarnings(
          stats::cor(
            proportion_main24,
            proportion_no_overlap21,
            method = "spearman"
          )
        ),
        pearson = suppressWarnings(
          stats::cor(
            proportion_main24,
            proportion_no_overlap21,
            method = "pearson"
          )
        ),
        mean_absolute_difference = mean(
          abs(
            proportion_main24 -
              proportion_no_overlap21
          )
        ),
        max_absolute_difference = max(
          abs(
            proportion_main24 -
              proportion_no_overlap21
          )
        )
      ),
      by = cell_type
    ]

    rosmap_sensitivity[
      ,
      cell_type_order :=
        match(
          cell_type,
          EXPECTED_CELL_TYPES
        )
    ]

    setorder(
      rosmap_sensitivity,
      cell_type_order
    )

    rosmap_sensitivity[
      ,
      cell_type_order := NULL
    ]

    fwrite(
      rosmap_sensitivity,
      rosmap_sensitivity_file,
      sep = "\t",
      na = "NA"
    )
  }
}


# ==============================================================================
# 12. SESSION INFO
# ==============================================================================

session_file <- file.path(
  result_root,
  "06_run_MuSiC2_sessionInfo.txt"
)

sink(
  session_file
)

print(
  sessionInfo()
)

sink()


# ==============================================================================
# 13. VEREDITO FINAL
# ==============================================================================

if (
  any(
    status_all$status ==
      "FAIL",
    na.rm = TRUE
  )
) {
  stop(
    "Existe pelo menos um check FAIL em results/deconvolution/MuSiC2/06_status.tsv."
  )
}


# ==============================================================================
# 14. RESUMO NO TERMINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "06_run_MuSiC2.R CONCLUÍDO\n"
)

cat(
  "============================================================\n\n"
)

print(
  run_summary_all[
    ,
    .(
      dataset,
      reference_version,
      marker_set,
      controls,
      cases,
      reference_donors,
      reference_nuclei,
      genes_used,
      music2_DE_genes_removed,
      converged,
      n_iter,
      n_id_not_converged,
      proportion_row_sum_min,
      proportion_row_sum_max,
      elapsed_seconds
    )
  ]
)

if (
  any(
    status_all$status ==
      "WARN",
    na.rm = TRUE
  )
) {

  cat(
    "\nExistem WARNs; revise:\n  ",
    file.path(
      result_root,
      "06_status.tsv"
    ),
    "\n",
    sep = ""
  )

} else {

  cat(
    "\nTodas as execuções concluíram sem WARN/FAIL.\n"
  )
}

if (
  file.exists(
    rosmap_sensitivity_file
  )
) {

  cat(
    "\nSensibilidade ROSMAP 24 vs 21 doadores:\n  ",
    rosmap_sensitivity_file,
    "\n",
    sep = ""
  )
}

cat(
  "\nOutputs consolidados:\n  ",
  result_root,
  "\n",
  sep = ""
)

cat(
  "\nPróxima etapa após revisão dos resultados:\n"
)

cat(
  "  07_run_BisqueRNA.R\n"
)
