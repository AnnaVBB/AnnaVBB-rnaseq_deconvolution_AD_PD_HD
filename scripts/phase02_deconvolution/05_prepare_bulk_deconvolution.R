#!/usr/bin/env Rscript

# ==============================================================================
# 05_prepare_bulk_deconvolution.R
# ==============================================================================
#
# FASE 02 — DECONVOLUÇÃO CELULAR
# ETAPA 05 — PREPARAÇÃO DOS BULKS PARA DECONVOLUÇÃO
#
# Objetivo
# --------
# Preparar, de forma padronizada e auditável, os datasets bulk RNA-seq que
# serão usados nos algoritmos de deconvolução.
#
# Datasets:
#
#   AD principal:
#     ROSMAP — DLPFC
#
#   AD validação independente:
#     GSE125583 — fusiform gyrus
#
#   HD:
#     GSE64810
#
#   PD:
#     GSE68719
#
#
# PRINCÍPIOS METODOLÓGICOS
# ------------------------
#
# 1. Cada bulk é harmonizado SEPARADAMENTE com Mathys 2019 por NCBI GeneID.
#
# 2. NÃO é criado um universo universal:
#
#       ROSMAP ∩ GSE125583 ∩ HD ∩ PD ∩ Mathys
#
#    Isso reduziria desnecessariamente a informação disponível.
#
#    Em vez disso:
#
#       ROSMAP    ∩ Mathys
#       GSE125583 ∩ Mathys
#       GSE64810  ∩ Mathys
#       GSE68719  ∩ Mathys
#
#    são calculados separadamente.
#
# 3. Genes ausentes em um bulk NÃO são preenchidos com zero.
#
# 4. Para os GEOs, a preferência é pela matriz raw counts NCBI, pois:
#
#       - preserva counts inteiros;
#       - preserva o universo gênico antes do filtro de DE;
#       - evita usar DESeq2-normalized counts como se fossem raw counts.
#
# 5. Para ROSMAP, o script procura uma matriz já harmonizada para GeneID
#    produzida na Fase 01c. Ele NÃO tenta converter silenciosamente a matriz
#    Ensembl raw do ROSMAP nesta etapa.
#
#    Isso é intencional: a harmonização Ensembl -> GeneID já foi auditada na
#    Fase 01c e não deve ser refeita de maneira diferente aqui.
#
# 6. CPM é calculado usando a library size da matriz bulk recebida ANTES de
#    restringir aos genes Mathys ou aos marcadores.
#
#    Portanto:
#
#        CPM(g,s) = count(g,s) / library_size_full(s) * 1e6
#
#    e NÃO:
#
#        count(g,s) / sum(counts somente nos genes intersectados) * 1e6
#
# 7. São preparados três níveis:
#
#       A. reference intersection
#          todos os GeneIDs em comum Mathys ∩ bulk
#
#       B. top100 marker intersection
#
#       C. top300 marker intersection   <- PRINCIPAL
#
#       D. top500 marker intersection
#
# 8. São salvas três escalas:
#
#       raw counts
#       CPM
#       log2(CPM + 1)
#
#    O uso final depende do algoritmo:
#
#       MuSiC2     -> raw counts / objeto single-cell apropriado
#       BisqueRNA  -> escala linear conforme API/protocolo a confirmar no 07
#       dtangle    -> log2(CPM + 1)
#
#    Esta etapa NÃO chama nenhum algoritmo de deconvolução.
#
#
# REFERÊNCIA MATHYS
# -----------------
#
# Usa os objetos finais construídos nas etapas 03 e 04:
#
#   - 17.394 GeneIDs expressos nos 24 no-pathology;
#   - 8 tipos celulares;
#   - top300 principal derivado por marker_score;
#   - top100/top500 para sensibilidade.
#
#
# PRIVACIDADE / ROSMAP
# --------------------
#
# Os identificadores de amostra são necessários para alinhar counts e metadata,
# mas esta etapa não exporta campos clínicos desnecessários do ROSMAP.
#
# O metadata de saída contém apenas:
#
#   sample_id
#   group_original
#   group
#   disease
#   cohort_role
#
#
# OUTPUT GERAL
# ------------
#
# data/processed/deconvolution/
#   ROSMAP/
#   GSE125583/
#   GSE64810/
#   GSE68719/
#
# Em cada pasta:
#
#   bulk_counts_reference_intersection.rds
#   bulk_cpm_reference_intersection.rds
#   bulk_log2cpm_reference_intersection.rds
#
#   bulk_counts_top100.rds
#   bulk_cpm_top100.rds
#   bulk_log2cpm_top100.rds
#
#   bulk_counts_top300.rds
#   bulk_cpm_top300.rds
#   bulk_log2cpm_top300.rds
#
#   bulk_counts_top500.rds
#   bulk_cpm_top500.rds
#   bulk_log2cpm_top500.rds
#
#   reference_signature_top100_matched.tsv
#   reference_signature_top300_matched.tsv
#   reference_signature_top500_matched.tsv
#
#   metadata_deconvolution.tsv
#   library_sizes.tsv
#   reference_intersection_geneids.tsv
#
#
# QC
# --
#
# results/deconvolution/bulk_preparation/
#
#   05_dataset_summary.tsv
#   05_marker_coverage_by_dataset_celltype.tsv
#   05_marker_set_coverage.tsv
#   05_input_resolution.tsv
#   05_status.tsv
#   05_output_manifest.tsv
#   05_prepare_bulk_deconvolution.txt
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
    paste(missing_packages, collapse = ", "),
    "\nInstale somente os pacotes ausentes no ambiente R atual. ",
    "Não altere o renv automaticamente."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})


# ==============================================================================
# 1. PARÂMETROS
# ==============================================================================

EXPECTED_REFERENCE_GENEIDS <- 17394L
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

EXPECTED_MARKERS_PER_TYPE <- c(
  top100 = 100L,
  top300 = 300L,
  top500 = 500L
)

PRIMARY_MARKER_SET <- "top300"

COUNT_TOLERANCE <- 1e-8

# Cobertura não é usada para excluir dataset automaticamente.
# Serve somente como QC.
MARKER_COVERAGE_WARN_THRESHOLD <- 0.80

# ROSMAP principal validado anteriormente.
EXPECTED_ROSMAP_MAIN_SAMPLES <- 144L

# GSE125583 completo validado anteriormente.
EXPECTED_GSE125583_SAMPLES <- 289L


# ==============================================================================
# 2. CAMINHOS
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

bulk_output_root <- file.path(
  "data",
  "processed",
  "deconvolution"
)

result_dir <- file.path(
  "results",
  "deconvolution",
  "bulk_preparation"
)

dir.create(
  bulk_output_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Referência
# ------------------------------------------------------------------------------

reference_counts_file <- file.path(
  reference_dir,
  "mathys2019_pathology_free_geneid_counts_sparse.rds"
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

signature_files <- list(
  top100 = file.path(
    markers_dir,
    "mathys2019_signature_top100_mean_cpm.tsv"
  ),
  top300 = file.path(
    markers_dir,
    "mathys2019_signature_top300_mean_cpm.tsv"
  ),
  top500 = file.path(
    markers_dir,
    "mathys2019_signature_top500_mean_cpm.tsv"
  )
)

step04_status_file <- file.path(
  "results",
  "deconvolution",
  "mathys2019",
  "marker_derivation",
  "04_marker_derivation_status.tsv"
)


# ------------------------------------------------------------------------------
# Metadata — candidatos explícitos
# ------------------------------------------------------------------------------

metadata_candidates <- list(

  ROSMAP = c(
    file.path(
      "results",
      "candidate_ROSMAP",
      "qc",
      "ROSMAP_metadata_main.tsv"
    ),
    file.path(
      "results",
      "candidate_ROSMAP",
      "qc",
      "metadata_ROSMAP_main.tsv"
    ),
    file.path(
      "results",
      "candidate_ROSMAP",
      "metadata",
      "ROSMAP_metadata_main.tsv"
    ),
    file.path(
      "results",
      "candidate_ROSMAP",
      "metadata",
      "ROSMAP_metadata.tsv"
    ),
    file.path(
      "metadata",
      "ROSMAP_metadata.tsv"
    ),
    file.path(
      "metadata",
      "ROSMAP_metadata.csv"
    )
  ),

  GSE125583 = c(
    file.path(
      "metadata",
      "GSE125583_metadata.csv"
    ),
    file.path(
      "results",
      "candidate_GSE125583",
      "metadata",
      "GSE125583_metadata.tsv"
    ),
    file.path(
      "results",
      "candidate_GSE125583",
      "metadata",
      "GSE125583_metadata.csv"
    )
  ),

  GSE64810 = c(
    file.path(
      "metadata",
      "GSE64810_metadata.csv"
    )
  ),

  GSE68719 = c(
    file.path(
      "metadata",
      "GSE68719_metadata.csv"
    )
  )
)


# ------------------------------------------------------------------------------
# Counts — candidatos explícitos
#
# Para GEOs, raw counts NCBI são prioridade.
# ------------------------------------------------------------------------------

count_candidates <- list(

  ROSMAP = c(
    file.path(
      "results",
      "candidate_ROSMAP",
      "qc",
      "counts_ROSMAP_GeneID.tsv"
    ),
    file.path(
      "results",
      "candidate_ROSMAP",
      "qc",
      "ROSMAP_counts_GeneID.tsv"
    ),
    file.path(
      "results",
      "candidate_ROSMAP",
      "qc",
      "counts_ROSMAP_harmonized.tsv"
    ),
    file.path(
      "data",
      "processed",
      "ROSMAP",
      "ROSMAP_counts_GeneID.tsv"
    ),
    file.path(
      "data",
      "processed",
      "ROSMAP_counts_GeneID.tsv"
    )
  ),

  GSE125583 = c(
    file.path(
      "data",
      "raw",
      "GSE125583",
      "GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz"
    ),
    file.path(
      "results",
      "candidate_GSE125583",
      "qc",
      "counts_GSE125583_DE_filtered.tsv"
    )
  ),

  GSE64810 = c(
    file.path(
      "data",
      "raw",
      "GSE64810",
      "GSE64810_raw_counts_GRCh38.p13_NCBI.tsv.gz"
    ),
    file.path(
      "results",
      "qc",
      "counts_HD_DE_filtered.tsv"
    )
  ),

  GSE68719 = c(
    file.path(
      "data",
      "raw",
      "GSE68719",
      "GSE68719_raw_counts_GRCh38.p13_NCBI.tsv.gz"
    ),
    file.path(
      "results",
      "qc",
      "counts_PD_DE_filtered.tsv"
    )
  )
)


# ------------------------------------------------------------------------------
# Diretórios de busca controlada.
# ------------------------------------------------------------------------------

metadata_search_dirs <- list(
  ROSMAP = c(
    file.path("results", "candidate_ROSMAP"),
    "metadata"
  ),
  GSE125583 = c(
    file.path("results", "candidate_GSE125583"),
    "metadata"
  ),
  GSE64810 = c(
    "metadata",
    file.path("results", "qc")
  ),
  GSE68719 = c(
    "metadata",
    file.path("results", "qc")
  )
)

count_search_dirs <- list(
  ROSMAP = c(
    file.path("results", "candidate_ROSMAP"),
    file.path("data", "processed", "ROSMAP"),
    file.path("data", "processed")
  ),
  GSE125583 = c(
    file.path("data", "raw", "GSE125583"),
    file.path("results", "candidate_GSE125583")
  ),
  GSE64810 = c(
    file.path("data", "raw", "GSE64810"),
    file.path("results", "qc")
  ),
  GSE68719 = c(
    file.path("data", "raw", "GSE68719"),
    file.path("results", "qc")
  )
)


# ==============================================================================
# 3. REGISTRO DOS DATASETS
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
# 4. FUNÇÕES AUXILIARES
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

  x <- fread(
    file,
    na.strings = c(
      "NA",
      ""
    )
  )

  if (!"status" %in% names(x)) {
    stop(
      label,
      ": coluna status ausente."
    )
  }

  if (any(x$status == "FAIL", na.rm = TRUE)) {
    stop(
      label,
      ": existe pelo menos um check FAIL."
    )
  }

  invisible(TRUE)
}


safe_file_size_mb <- function(file) {

  if (
    is.na(file) ||
      !file.exists(file)
  ) {
    return(NA_real_)
  }

  as.numeric(
    file.info(file)$size
  ) /
    1024^2
}


# ------------------------------------------------------------------------------
# Ler somente cabeçalho de TSV/CSV.
# ------------------------------------------------------------------------------

read_header_names <- function(file) {

  if (
    grepl(
      "\\.rds$",
      file,
      ignore.case = TRUE
    )
  ) {
    return(character())
  }

  x <- tryCatch(
    fread(
      file,
      nrows = 2L,
      na.strings = c(
        "NA",
        ""
      ),
      showProgress = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(x)) {
    return(character())
  }

  names(x)
}


# ------------------------------------------------------------------------------
# Localizar metadata de maneira controlada.
# ------------------------------------------------------------------------------

resolve_metadata_file <- function(dataset) {

  explicit <- metadata_candidates[[dataset]]

  explicit_existing <- explicit[
    file.exists(explicit)
  ]

  if (length(explicit_existing) > 0L) {
    return(
      list(
        file = explicit_existing[1],
        resolution = "explicit_candidate"
      )
    )
  }

  dirs <- metadata_search_dirs[[dataset]]

  dirs <- dirs[
    dir.exists(dirs)
  ]

  if (length(dirs) == 0L) {
    stop(
      dataset,
      ": nenhum diretório disponível para procurar metadata."
    )
  }

  files <- unique(
    unlist(
      lapply(
        dirs,
        function(d) {
          list.files(
            d,
            pattern = "\\.(tsv|csv|txt)$",
            full.names = TRUE,
            recursive = TRUE,
            ignore.case = TRUE
          )
        }
      )
    )
  )

  files <- files[
    grepl(
      paste0(
        dataset,
        "|metadata|sample"
      ),
      basename(files),
      ignore.case = TRUE
    )
  ]

  if (length(files) == 0L) {
    stop(
      dataset,
      ": metadata não encontrado automaticamente."
    )
  }

  valid <- logical(
    length(files)
  )

  scores <- numeric(
    length(files)
  )

  for (i in seq_along(files)) {

    nm <- read_header_names(
      files[i]
    )

    has_sample <- any(
      c(
        "sample_id",
        "geo_accession",
        "sample",
        "specimenID"
      ) %in%
        nm
    )

    has_group <- any(
      c(
        "group",
        "diagnosis",
        "raw_diagnosis",
        "cogdx"
      ) %in%
        nm
    )

    valid[i] <- has_sample &&
      has_group

    scores[i] <- sum(
      c(
        grepl(
          dataset,
          basename(files[i]),
          ignore.case = TRUE
        ),
        grepl(
          "metadata",
          basename(files[i]),
          ignore.case = TRUE
        ),
        "sample_id" %in% nm,
        "group" %in% nm
      )
    )
  }

  candidates <- files[
    valid
  ]

  candidate_scores <- scores[
    valid
  ]

  if (length(candidates) == 0L) {
    stop(
      dataset,
      ": arquivos foram encontrados, mas nenhum parece metadata válido ",
      "(sample ID + grupo/diagnóstico)."
    )
  }

  ord <- order(
    -candidate_scores,
    nchar(candidates),
    candidates
  )

  list(
    file = candidates[ord][1],
    resolution = "automatic_search"
  )
}


# ------------------------------------------------------------------------------
# Localizar counts.
#
# Para evitar escolher resultado DESeq2 por engano:
#
#   - arquivo precisa ter muitas colunas OU ser RDS;
#   - nomes "deseq", "deg", "result", "summary" recebem exclusão;
#   - preferência por nomes com count/matrix/harmon/geneid.
# ------------------------------------------------------------------------------

resolve_count_file <- function(
  dataset,
  metadata
) {

  explicit <- count_candidates[[dataset]]

  explicit_existing <- explicit[
    file.exists(explicit)
  ]

  if (length(explicit_existing) > 0L) {
    return(
      list(
        file = explicit_existing[1],
        resolution = "explicit_candidate"
      )
    )
  }

  dirs <- count_search_dirs[[dataset]]

  dirs <- dirs[
    dir.exists(dirs)
  ]

  if (length(dirs) == 0L) {
    stop(
      dataset,
      ": nenhum diretório disponível para procurar counts."
    )
  }

  files <- unique(
    unlist(
      lapply(
        dirs,
        function(d) {
          list.files(
            d,
            pattern = "\\.(tsv|csv|txt|rds)(\\.gz)?$",
            full.names = TRUE,
            recursive = TRUE,
            ignore.case = TRUE
          )
        }
      )
    )
  )

  if (length(files) == 0L) {
    stop(
      dataset,
      ": nenhum arquivo tabular/RDS encontrado para counts."
    )
  }

  bad_name <- grepl(
    "deseq|deg_|deg\\.|result|summary|enrich|marker|metadata|clinical|biospecimen|assay|session|status|report",
    basename(files),
    ignore.case = TRUE
  )

  files <- files[
    !bad_name
  ]

  if (length(files) == 0L) {
    stop(
      dataset,
      ": todos os arquivos encontrados foram excluídos como resultados/metadata."
    )
  }

  sample_candidates <- unique(
    na.omit(
      c(
        if ("sample_id" %in% names(metadata)) {
          metadata$sample_id
        } else {
          NULL
        },
        if ("geo_accession" %in% names(metadata)) {
          metadata$geo_accession
        } else {
          NULL
        },
        if ("specimenID" %in% names(metadata)) {
          metadata$specimenID
        } else {
          NULL
        }
      )
    )
  )

  audit <- vector(
    "list",
    length(files)
  )

  for (i in seq_along(files)) {

    f <- files[i]

    if (
      grepl(
        "\\.rds$",
        f,
        ignore.case = TRUE
      )
    ) {

      audit[[i]] <- data.table(
        file = f,
        n_header_columns = NA_integer_,
        sample_overlap = NA_integer_,
        score = as.numeric(
          2 +
            grepl(
              "count|matrix|harmon|geneid",
              basename(f),
              ignore.case = TRUE
            )
        )
      )

      next
    }

    nm <- read_header_names(
      f
    )

    if (length(nm) < 3L) {

      audit[[i]] <- data.table(
        file = f,
        n_header_columns = length(nm),
        sample_overlap = 0L,
        score = -Inf
      )

      next
    }

    overlap_n <- length(
      intersect(
        nm,
        sample_candidates
      )
    )

    name_score <- sum(
      c(
        grepl(
          "count",
          basename(f),
          ignore.case = TRUE
        ),
        grepl(
          "matrix",
          basename(f),
          ignore.case = TRUE
        ),
        grepl(
          "harmon",
          basename(f),
          ignore.case = TRUE
        ),
        grepl(
          "geneid",
          basename(f),
          ignore.case = TRUE
        ),
        grepl(
          dataset,
          basename(f),
          ignore.case = TRUE
        )
      )
    )

    audit[[i]] <- data.table(
      file = f,
      n_header_columns = as.integer(
        length(nm)
      ),
      sample_overlap = as.integer(
        overlap_n
      ),
      score = as.numeric(
        10 * overlap_n +
          name_score +
          log10(
            max(
              length(nm),
              1
            )
          )
      )
    )
  }

  audit_dt <- rbindlist(
    audit,
    use.names = TRUE,
    fill = TRUE
  )

  audit_dt <- audit_dt[
    is.finite(score)
  ]

  if (nrow(audit_dt) == 0L) {
    stop(
      dataset,
      ": nenhum candidato plausível de counts foi encontrado."
    )
  }

  setorder(
    audit_dt,
    -score,
    file
  )

  fwrite(
    audit_dt,
    file.path(
      result_dir,
      paste0(
        "05_",
        dataset,
        "_count_discovery_candidates.tsv"
      )
    ),
    sep = "\t",
    na = "NA"
  )

  list(
    file = audit_dt$file[1],
    resolution = "automatic_search"
  )
}


# ------------------------------------------------------------------------------
# Detectar coluna GeneID.
# ------------------------------------------------------------------------------

detect_geneid_column <- function(dt) {

  candidates <- c(
    "GeneID",
    "gene_id",
    "geneid",
    "ENTREZID",
    "entrez_id"
  )

  hit <- intersect(
    candidates,
    names(dt)
  )

  if (length(hit) > 0L) {
    return(hit[1])
  }

  # Algumas matrizes NCBI podem ter GeneID na primeira coluna com outro nome.
  first_col <- names(dt)[1]

  x <- clean_chr(
    dt[[first_col]]
  )

  valid_fraction <- mean(
    grepl(
      "^[0-9]+$",
      x
    ),
    na.rm = TRUE
  )

  if (
    is.finite(valid_fraction) &&
      valid_fraction > 0.90
  ) {
    return(first_col)
  }

  stop(
    "Não foi possível identificar uma coluna NCBI GeneID."
  )
}


# ------------------------------------------------------------------------------
# Ler matriz de counts.
# ------------------------------------------------------------------------------

read_counts_table <- function(file) {

  if (
    grepl(
      "\\.rds$",
      file,
      ignore.case = TRUE
    )
  ) {

    obj <- readRDS(
      file
    )

    if (
      is.matrix(obj) ||
        inherits(
          obj,
          "Matrix"
        )
    ) {

      if (
        is.null(
          rownames(obj)
        ) ||
          is.null(
            colnames(obj)
          )
      ) {
        stop(
          "RDS de counts sem rownames/colnames: ",
          file
        )
      }

      gene_ids <- clean_chr(
        rownames(obj)
      )

      if (
        anyNA(gene_ids) ||
          anyDuplicated(gene_ids)
      ) {
        stop(
          "RDS possui GeneIDs ausentes/duplicados: ",
          file
        )
      }

      counts <- as.matrix(
        obj
      )

      return(
        list(
          gene_ids = gene_ids,
          counts = counts,
          gene_id_source = "rownames"
        )
      )
    }

    if (is.data.frame(obj)) {
      dt <- as.data.table(obj)
    } else {
      stop(
        "RDS não contém matriz/data.frame reconhecido: ",
        file
      )
    }

  } else {

    dt <- fread(
      file,
      na.strings = c(
        "NA",
        ""
      ),
      check.names = FALSE,
      showProgress = TRUE
    )
  }

  gene_col <- detect_geneid_column(
    dt
  )

  gene_ids <- clean_chr(
    dt[[gene_col]]
  )

  if (anyNA(gene_ids)) {
    stop(
      "Existem GeneIDs ausentes em ",
      file
    )
  }

  if (
    any(
      !grepl(
        "^[0-9]+$",
        gene_ids
      )
    )
  ) {
    stop(
      "Existem identificadores não numéricos na coluna GeneID de ",
      file
    )
  }

  if (anyDuplicated(gene_ids)) {
    stop(
      "Existem GeneIDs duplicados em ",
      file
    )
  }

  sample_cols <- setdiff(
    names(dt),
    gene_col
  )

  if (length(sample_cols) < 2L) {
    stop(
      "Matriz possui menos de duas colunas de amostra: ",
      file
    )
  }

  counts <- as.matrix(
    dt[
      ,
      ..sample_cols
    ]
  )

  storage.mode(counts) <- "double"

  rownames(counts) <- gene_ids
  colnames(counts) <- sample_cols

  rm(dt)
  invisible(gc())

  list(
    gene_ids = gene_ids,
    counts = counts,
    gene_id_source = gene_col
  )
}


# ------------------------------------------------------------------------------
# Validar raw counts.
# ------------------------------------------------------------------------------

validate_raw_counts <- function(
  counts,
  dataset
) {

  if (anyNA(counts)) {
    stop(
      dataset,
      ": counts contêm NA."
    )
  }

  if (
    any(
      !is.finite(counts)
    )
  ) {
    stop(
      dataset,
      ": counts contêm valores não finitos."
    )
  }

  if (
    any(
      counts < 0
    )
  ) {
    stop(
      dataset,
      ": counts contêm valores negativos."
    )
  }

  fractional <- abs(
    counts -
      round(counts)
  ) >
    COUNT_TOLERANCE

  if (any(fractional)) {
    stop(
      dataset,
      ": a matriz selecionada possui valores fracionários. ",
      "Isso sugere counts normalizados/transformados, não raw counts. ",
      "Não prosseguir com MuSiC2 usando este arquivo."
    )
  }

  if (
    any(
      colSums(counts) <= 0
    )
  ) {
    stop(
      dataset,
      ": existe amostra com library size <= 0."
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# Ler metadata.
# ------------------------------------------------------------------------------

read_metadata <- function(file) {

  fread(
    file,
    na.strings = c(
      "NA",
      ""
    ),
    check.names = FALSE
  )
}


# ------------------------------------------------------------------------------
# Detectar coluna de amostra principal no metadata.
# ------------------------------------------------------------------------------

detect_metadata_sample_column <- function(metadata) {

  candidates <- c(
    "sample_id",
    "sample",
    "specimenID",
    "specimen_id",
    "geo_accession"
  )

  hit <- intersect(
    candidates,
    names(metadata)
  )

  if (length(hit) == 0L) {
    stop(
      "Nenhuma coluna de sample ID reconhecida no metadata."
    )
  }

  hit[1]
}


# ------------------------------------------------------------------------------
# Detectar coluna original de grupo/diagnóstico.
# ------------------------------------------------------------------------------

detect_group_column <- function(metadata) {

  candidates <- c(
    "group",
    "diagnosis",
    "raw_diagnosis",
    "disease_status",
    "disease",
    "cogdx"
  )

  hit <- intersect(
    candidates,
    names(metadata)
  )

  if (length(hit) == 0L) {
    stop(
      "Nenhuma coluna de grupo/diagnóstico reconhecida no metadata."
    )
  }

  hit[1]
}


# ------------------------------------------------------------------------------
# Normalizar grupo para Control / AD / HD / PD.
# ------------------------------------------------------------------------------

normalize_group <- function(
  x,
  disease
) {

  y <- tolower(
    clean_chr(x)
  )

  out <- rep(
    NA_character_,
    length(y)
  )

  control_patterns <- c(
    "^control$",
    "^ctrl$",
    "control",
    "normal",
    "unaffected",
    "^1$"
  )

  disease_patterns <- switch(
    disease,
    AD = c(
      "^ad$",
      "alzheimer",
      "^4$"
    ),
    HD = c(
      "^hd$",
      "huntington"
    ),
    PD = c(
      "^pd$",
      "parkinson"
    ),
    stop(
      "Doença não reconhecida: ",
      disease
    )
  )

  for (p in control_patterns) {
    out[
      is.na(out) &
        !is.na(y) &
        grepl(
          p,
          y,
          ignore.case = TRUE
        )
    ] <- "Control"
  }

  for (p in disease_patterns) {
    out[
      is.na(out) &
        !is.na(y) &
        grepl(
          p,
          y,
          ignore.case = TRUE
        )
    ] <- disease
  }

  out
}


# ------------------------------------------------------------------------------
# Alinhar matriz e metadata.
#
# Estratégias:
#
#   1) counts colnames == sample_id
#   2) counts colnames == geo_accession, renomeando para sample_id
#   3) counts colnames == specimenID, renomeando para sample_id
#
# O resultado sempre usa sample_id nas colunas.
# ------------------------------------------------------------------------------

align_counts_metadata <- function(
  counts,
  metadata,
  dataset,
  disease,
  cohort_role
) {

  metadata <- copy(
    metadata
  )

  primary_sample_col <- detect_metadata_sample_column(
    metadata
  )

  metadata[
    ,
    sample_id_internal :=
      clean_chr(
        get(
          primary_sample_col
        )
      )
  ]

  if (
    anyNA(
      metadata$sample_id_internal
    )
  ) {
    stop(
      dataset,
      ": sample ID ausente no metadata."
    )
  }

  if (
    anyDuplicated(
      metadata$sample_id_internal
    )
  ) {
    stop(
      dataset,
      ": sample IDs duplicados no metadata."
    )
  }

  group_col <- detect_group_column(
    metadata
  )

  metadata[
    ,
    group_original :=
      clean_chr(
        get(
          group_col
        )
      )
  ]

  metadata[
    ,
    group :=
      normalize_group(
        group_original,
        disease
      )
  ]

  # Em metadata ROSMAP que ainda contenha outros cogdx, manter somente
  # Control/AD definidos no pipeline principal.
  metadata <- metadata[
    !is.na(group)
  ]

  if (nrow(metadata) < 2L) {
    stop(
      dataset,
      ": menos de duas amostras após normalização dos grupos."
    )
  }

  count_ids <- colnames(
    counts
  )

  # ---------------------------------------------------------------------------
  # Estratégia 1: sample ID principal.
  # ---------------------------------------------------------------------------

  direct_common <- intersect(
    count_ids,
    metadata$sample_id_internal
  )

  best_mode <- "sample_id"
  best_common <- direct_common
  mapping_col <- "sample_id_internal"

  # ---------------------------------------------------------------------------
  # Estratégia 2: geo_accession.
  # ---------------------------------------------------------------------------

  if (
    "geo_accession" %in%
      names(metadata)
  ) {

    geo <- clean_chr(
      metadata$geo_accession
    )

    geo_common <- intersect(
      count_ids,
      geo
    )

    if (
      length(geo_common) >
        length(best_common)
    ) {

      best_mode <- "geo_accession"
      best_common <- geo_common
      mapping_col <- "geo_accession"
    }
  }

  # ---------------------------------------------------------------------------
  # Estratégia 3: specimenID.
  # ---------------------------------------------------------------------------

  if (
    "specimenID" %in%
      names(metadata)
  ) {

    specimen <- clean_chr(
      metadata$specimenID
    )

    specimen_common <- intersect(
      count_ids,
      specimen
    )

    if (
      length(specimen_common) >
        length(best_common)
    ) {

      best_mode <- "specimenID"
      best_common <- specimen_common
      mapping_col <- "specimenID"
    }
  }

  if (length(best_common) < 2L) {
    stop(
      dataset,
      ": não foi possível alinhar counts e metadata. ",
      "Maior overlap encontrado: ",
      length(best_common),
      "."
    )
  }

  metadata[
    ,
    mapping_value :=
      clean_chr(
        get(
          mapping_col
        )
      )
  ]

  metadata <- metadata[
    mapping_value %in%
      best_common
  ]

  if (
    anyDuplicated(
      metadata$mapping_value
    )
  ) {
    stop(
      dataset,
      ": chave de alinhamento duplicada no metadata: ",
      mapping_col,
      "."
    )
  }

  # Manter a ordem do metadata e buscar no counts.
  count_idx <- match(
    metadata$mapping_value,
    colnames(counts)
  )

  if (anyNA(count_idx)) {
    stop(
      dataset,
      ": falha interna ao alinhar counts."
    )
  }

  counts <- counts[
    ,
    count_idx,
    drop = FALSE
  ]

  # Renomear counts para o ID interno escolhido no metadata.
  colnames(counts) <- metadata$sample_id_internal

  # Metadata mínimo para deconvolução.
  metadata_out <- data.table(
    sample_id = metadata$sample_id_internal,
    group_original = metadata$group_original,
    group = metadata$group,
    disease = disease,
    cohort_role = cohort_role
  )

  if (
    anyDuplicated(
      metadata_out$sample_id
    )
  ) {
    stop(
      dataset,
      ": sample_id duplicado após alinhamento."
    )
  }

  if (
    !identical(
      colnames(counts),
      metadata_out$sample_id
    )
  ) {
    stop(
      dataset,
      ": ordem counts/metadata não coincide."
    )
  }

  group_levels <- unique(
    metadata_out$group
  )

  if (
    !all(
      c(
        "Control",
        disease
      ) %in%
        group_levels
    )
  ) {
    stop(
      dataset,
      ": grupos esperados Control e ",
      disease,
      " não estão ambos presentes. Encontrados: ",
      paste(
        group_levels,
        collapse = ", "
      )
    )
  }

  list(
    counts = counts,
    metadata = metadata_out,
    alignment_mode = best_mode,
    n_metadata_before_group_filter = nrow(
      metadata
    )
  )
}


# ------------------------------------------------------------------------------
# CPM usando library sizes da matriz COMPLETA recebida.
# ------------------------------------------------------------------------------

counts_to_cpm <- function(
  counts,
  library_sizes
) {

  if (
    length(library_sizes) !=
      ncol(counts)
  ) {
    stop(
      "counts_to_cpm: tamanho de library_sizes incompatível."
    )
  }

  if (
    any(
      library_sizes <= 0
    )
  ) {
    stop(
      "counts_to_cpm: library size <= 0."
    )
  }

  sweep(
    counts,
    2,
    library_sizes,
    "/"
  ) *
    1e6
}


# ------------------------------------------------------------------------------
# Salvar matriz como RDS sem compressão.
# ------------------------------------------------------------------------------

save_matrix_rds <- function(
  x,
  file
) {

  saveRDS(
    x,
    file,
    compress = FALSE
  )
}


# ------------------------------------------------------------------------------
# Preparar signature matrix da referência para genes disponíveis no bulk.
# ------------------------------------------------------------------------------

prepare_matched_signature <- function(
  signature_dt,
  available_geneids,
  output_file
) {

  sig <- copy(
    signature_dt
  )

  if (
    !"gene_id" %in%
      names(sig)
  ) {
    stop(
      "Signature matrix sem gene_id."
    )
  }

  sig[
    ,
    gene_id :=
      clean_chr(
        gene_id
      )
  ]

  sig <- sig[
    gene_id %in%
      available_geneids
  ]

  # Preservar ordem do bulk, não a ordem original da signature.
  sig <- sig[
    match(
      available_geneids,
      gene_id
    )
  ]

  if (
    anyNA(
      sig$gene_id
    )
  ) {
    stop(
      "Falha ao alinhar signature aos genes disponíveis."
    )
  }

  fwrite(
    sig,
    output_file,
    sep = "\t",
    na = "NA"
  )

  sig
}


# ------------------------------------------------------------------------------
# Cobertura por tipo celular.
#
# Um mesmo GeneID pode pertencer a mais de uma lista de tipos.
# Portanto, coverage por cell type usa as ATRIBUIÇÕES topK, não somente a união.
# ------------------------------------------------------------------------------

marker_coverage_by_type <- function(
  marker_dt,
  bulk_geneids,
  dataset,
  marker_set
) {

  x <- copy(
    marker_dt
  )

  x[
    ,
    available_in_bulk :=
      gene_id %in%
        bulk_geneids
  ]

  out <- x[
    ,
    .(
      expected_assignments = .N,
      available_assignments = sum(
        available_in_bulk
      ),
      missing_assignments = sum(
        !available_in_bulk
      ),
      coverage = mean(
        available_in_bulk
      )
    ),
    by = cell_type
  ]

  out[
    ,
    `:=`(
      dataset = dataset,
      marker_set = marker_set
    )
  ]

  setcolorder(
    out,
    c(
      "dataset",
      "marker_set",
      "cell_type",
      "expected_assignments",
      "available_assignments",
      "missing_assignments",
      "coverage"
    )
  )

  out
}


# ==============================================================================
# 5. CARREGAR E VALIDAR REFERÊNCIA
# ==============================================================================

cat(
  "============================================================\n"
)

cat(
  "FASE 02 — 05_prepare_bulk_deconvolution.R\n"
)

cat(
  "PREPARAÇÃO DOS BULKS PARA DECONVOLUÇÃO\n"
)

cat(
  "============================================================\n\n"
)

validate_status_file(
  step04_status_file,
  "Etapa 04"
)

reference_required_files <- c(
  reference_counts_file,
  unlist(
    marker_files,
    use.names = FALSE
  ),
  unlist(
    signature_files,
    use.names = FALSE
  )
)

missing_reference_files <- reference_required_files[
  !file.exists(
    reference_required_files
  )
]

if (
  length(
    missing_reference_files
  ) > 0L
) {

  stop(
    "Arquivos da referência ausentes:\n",
    paste(
      missing_reference_files,
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------------------------
# Universo Mathys por GeneID
# ------------------------------------------------------------------------------

reference_counts <- readRDS(
  reference_counts_file
)

if (
  is.null(
    rownames(
      reference_counts
    )
  )
) {
  stop(
    "Matriz Mathys nucleus-level sem rownames GeneID."
  )
}

reference_geneids <- clean_chr(
  rownames(
    reference_counts
  )
)

if (
  length(
    reference_geneids
  ) !=
    EXPECTED_REFERENCE_GENEIDS
) {
  stop(
    "Esperávamos ",
    EXPECTED_REFERENCE_GENEIDS,
    " GeneIDs na referência; encontrados: ",
    length(
      reference_geneids
    ),
    "."
  )
}

if (
  anyNA(reference_geneids) ||
    anyDuplicated(reference_geneids)
) {
  stop(
    "GeneIDs da referência ausentes ou duplicados."
  )
}

rm(reference_counts)
invisible(gc())


# ------------------------------------------------------------------------------
# Marcadores e signatures
# ------------------------------------------------------------------------------

markers <- list()
signatures <- list()

for (
  marker_set in
  names(marker_files)
) {

  marker_dt <- fread(
    marker_files[[marker_set]],
    na.strings = c(
      "NA",
      ""
    )
  )

  required_marker_cols <- c(
    "gene_id",
    "cell_type",
    "marker_rank"
  )

  missing_cols <- setdiff(
    required_marker_cols,
    names(marker_dt)
  )

  if (
    length(
      missing_cols
    ) > 0L
  ) {
    stop(
      marker_set,
      ": colunas ausentes no arquivo de marcadores: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }

  marker_dt[
    ,
    gene_id :=
      clean_chr(
        gene_id
      )
  ]

  marker_dt[
    ,
    cell_type :=
      clean_chr(
        cell_type
      )
  ]

  if (
    !setequal(
      unique(
        marker_dt$cell_type
      ),
      EXPECTED_CELL_TYPES
    )
  ) {
    stop(
      marker_set,
      ": tipos celulares inesperados."
    )
  }

  expected_k <- EXPECTED_MARKERS_PER_TYPE[
    marker_set
  ]

  marker_counts <- marker_dt[
    ,
    .N,
    by = cell_type
  ]

  if (
    any(
      marker_counts$N !=
        expected_k
    )
  ) {
    stop(
      marker_set,
      ": número de marcadores por tipo diferente do esperado."
    )
  }

  if (
    any(
      !marker_dt$gene_id %in%
        reference_geneids
    )
  ) {
    stop(
      marker_set,
      ": existem marcadores fora do universo Mathys."
    )
  }

  markers[[marker_set]] <- marker_dt

  sig <- fread(
    signature_files[[marker_set]],
    na.strings = c(
      "NA",
      ""
    )
  )

  if (
    !"gene_id" %in%
      names(sig)
  ) {
    stop(
      marker_set,
      ": signature matrix sem gene_id."
    )
  }

  sig[
    ,
    gene_id :=
      clean_chr(
        gene_id
      )
  ]

  signatures[[marker_set]] <- sig
}

cat(
  "Referência carregada: ",
  length(reference_geneids),
  " GeneIDs.\n",
  sep = ""
)

cat(
  "União top300: ",
  uniqueN(
    markers$top300$gene_id
  ),
  " GeneIDs.\n\n",
  sep = ""
)


# ==============================================================================
# 6. PREPARAR CADA DATASET
# ==============================================================================

dataset_summaries <- list()
marker_coverage_rows <- list()
marker_set_summary_rows <- list()
input_resolution_rows <- list()
status_rows <- list()
output_manifest_rows <- list()

summary_counter <- 1L
coverage_counter <- 1L
marker_set_counter <- 1L
input_counter <- 1L
status_counter <- 1L
manifest_counter <- 1L


for (
  dataset_i in
  dataset_registry$dataset
) {

  disease_i <- dataset_registry[
    dataset == dataset_i,
    disease
  ]

  cohort_role_i <- dataset_registry[
    dataset == dataset_i,
    cohort_role
  ]

  cat(
    "\n============================================================\n"
  )

  cat(
    dataset_i,
    " — ",
    disease_i,
    " — ",
    cohort_role_i,
    "\n",
    sep = ""
  )

  cat(
    "============================================================\n"
  )


  # ============================================================================
  # 6.1 Resolver metadata
  # ============================================================================

  metadata_resolution <- resolve_metadata_file(
    dataset_i
  )

  metadata_file <- metadata_resolution$file

  cat(
    "Metadata:\n  ",
    metadata_file,
    "\n",
    sep = ""
  )

  metadata_raw <- read_metadata(
    metadata_file
  )


  # ============================================================================
  # 6.2 Resolver counts
  # ============================================================================

  count_resolution <- resolve_count_file(
    dataset_i,
    metadata_raw
  )

  count_file <- count_resolution$file

  cat(
    "Counts:\n  ",
    count_file,
    "\n",
    sep = ""
  )

  input_resolution_rows[[input_counter]] <- data.table(
    dataset = dataset_i,
    input_type = c(
      "metadata",
      "counts"
    ),
    file = c(
      metadata_file,
      count_file
    ),
    resolution = c(
      metadata_resolution$resolution,
      count_resolution$resolution
    ),
    size_MB = c(
      safe_file_size_mb(
        metadata_file
      ),
      safe_file_size_mb(
        count_file
      )
    )
  )

  input_counter <- input_counter + 1L


  # ============================================================================
  # 6.3 Ler counts
  # ============================================================================

  counts_obj <- read_counts_table(
    count_file
  )

  counts_full <- counts_obj$counts

  validate_raw_counts(
    counts_full,
    dataset_i
  )

  cat(
    "Matriz recebida: ",
    nrow(
      counts_full
    ),
    " GeneIDs x ",
    ncol(
      counts_full
    ),
    " colunas de amostra.\n",
    sep = ""
  )


  # ============================================================================
  # 6.4 Alinhar counts e metadata
  # ============================================================================

  aligned <- align_counts_metadata(
    counts = counts_full,
    metadata = metadata_raw,
    dataset = dataset_i,
    disease = disease_i,
    cohort_role = cohort_role_i
  )

  counts_full <- aligned$counts
  metadata <- aligned$metadata
  alignment_mode_i <- aligned$alignment_mode

  rm(
    metadata_raw,
    aligned
  )

  invisible(
    gc()
  )

  cat(
    "Amostras após alinhamento/grupo: ",
    ncol(
      counts_full
    ),
    "\n",
    sep = ""
  )

  print(
    table(
      metadata$group
    )
  )


  # ============================================================================
  # 6.5 Sanity checks de tamanho conhecidos
  # ============================================================================

  if (
    dataset_i == "ROSMAP" &&
      ncol(
        counts_full
      ) !=
        EXPECTED_ROSMAP_MAIN_SAMPLES
  ) {

    warning(
      "ROSMAP: esperávamos ",
      EXPECTED_ROSMAP_MAIN_SAMPLES,
      " amostras no fenótipo principal, mas foram alinhadas ",
      ncol(
        counts_full
      ),
      ". Verifique se o metadata selecionado é o metadata principal da Fase 01c."
    )
  }

  if (
    dataset_i == "GSE125583" &&
      ncol(
        counts_full
      ) !=
        EXPECTED_GSE125583_SAMPLES
  ) {

    warning(
      "GSE125583: esperávamos ",
      EXPECTED_GSE125583_SAMPLES,
      " amostras, mas foram alinhadas ",
      ncol(
        counts_full
      ),
      "."
    )
  }


  # ============================================================================
  # 6.6 Library sizes ANTES da interseção
  # ============================================================================

  library_sizes <- colSums(
    counts_full
  )

  if (
    any(
      library_sizes <= 0
    )
  ) {
    stop(
      dataset_i,
      ": library size <= 0."
    )
  }

  library_dt <- data.table(
    sample_id = colnames(
      counts_full
    ),
    library_size_full_input = as.numeric(
      library_sizes
    )
  )


  # ============================================================================
  # 6.7 Interseção Mathys ∩ bulk
  # ============================================================================

  bulk_geneids <- rownames(
    counts_full
  )

  common_reference_geneids <- reference_geneids[
    reference_geneids %in%
      bulk_geneids
  ]

  if (
    length(
      common_reference_geneids
    ) == 0L
  ) {
    stop(
      dataset_i,
      ": nenhuma interseção GeneID com Mathys."
    )
  }

  reference_coverage <- length(
    common_reference_geneids
  ) /
    length(
      reference_geneids
    )

  counts_reference <- counts_full[
    common_reference_geneids,
    ,
    drop = FALSE
  ]

  cpm_reference <- counts_to_cpm(
    counts_reference,
    library_sizes
  )

  log2cpm_reference <- log2(
    cpm_reference +
      1
  )


  # ============================================================================
  # 6.8 Pasta do dataset
  # ============================================================================

  dataset_out_dir <- file.path(
    bulk_output_root,
    dataset_i
  )

  dir.create(
    dataset_out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  # ============================================================================
  # 6.9 Salvar metadata e library sizes
  # ============================================================================

  metadata_file_out <- file.path(
    dataset_out_dir,
    "metadata_deconvolution.tsv"
  )

  library_file_out <- file.path(
    dataset_out_dir,
    "library_sizes.tsv"
  )

  intersection_file_out <- file.path(
    dataset_out_dir,
    "reference_intersection_geneids.tsv"
  )

  fwrite(
    metadata,
    metadata_file_out,
    sep = "\t",
    na = "NA"
  )

  fwrite(
    library_dt,
    library_file_out,
    sep = "\t",
    na = "NA"
  )

  fwrite(
    data.table(
      gene_id = common_reference_geneids
    ),
    intersection_file_out,
    sep = "\t"
  )


  # ============================================================================
  # 6.10 Salvar interseção completa da referência
  # ============================================================================

  ref_counts_out <- file.path(
    dataset_out_dir,
    "bulk_counts_reference_intersection.rds"
  )

  ref_cpm_out <- file.path(
    dataset_out_dir,
    "bulk_cpm_reference_intersection.rds"
  )

  ref_log2cpm_out <- file.path(
    dataset_out_dir,
    "bulk_log2cpm_reference_intersection.rds"
  )

  save_matrix_rds(
    counts_reference,
    ref_counts_out
  )

  save_matrix_rds(
    cpm_reference,
    ref_cpm_out
  )

  save_matrix_rds(
    log2cpm_reference,
    ref_log2cpm_out
  )


  # ============================================================================
  # 6.11 Preparar top100 / top300 / top500
  # ============================================================================

  marker_outputs <- list()

  for (
    marker_set_i in
    c(
      "top100",
      "top300",
      "top500"
    )
  ) {

    marker_dt <- markers[[marker_set_i]]

    marker_union <- unique(
      marker_dt$gene_id
    )

    available_markers <- marker_union[
      marker_union %in%
        bulk_geneids
    ]

    missing_markers <- setdiff(
      marker_union,
      available_markers
    )

    # Preservar ordem da referência/signature.
    signature_gene_order <- signatures[[marker_set_i]]$gene_id

    available_markers <- signature_gene_order[
      signature_gene_order %in%
        available_markers
    ]

    if (
      length(
        available_markers
      ) == 0L
    ) {
      stop(
        dataset_i,
        ": nenhum marcador disponível para ",
        marker_set_i,
        "."
      )
    }

    counts_marker <- counts_full[
      available_markers,
      ,
      drop = FALSE
    ]

    cpm_marker <- counts_to_cpm(
      counts_marker,
      library_sizes
    )

    log2cpm_marker <- log2(
      cpm_marker +
        1
    )

    counts_marker_out <- file.path(
      dataset_out_dir,
      paste0(
        "bulk_counts_",
        marker_set_i,
        ".rds"
      )
    )

    cpm_marker_out <- file.path(
      dataset_out_dir,
      paste0(
        "bulk_cpm_",
        marker_set_i,
        ".rds"
      )
    )

    log2cpm_marker_out <- file.path(
      dataset_out_dir,
      paste0(
        "bulk_log2cpm_",
        marker_set_i,
        ".rds"
      )
    )

    signature_matched_out <- file.path(
      dataset_out_dir,
      paste0(
        "reference_signature_",
        marker_set_i,
        "_matched.tsv"
      )
    )

    missing_marker_out <- file.path(
      dataset_out_dir,
      paste0(
        "missing_markers_",
        marker_set_i,
        ".tsv"
      )
    )

    save_matrix_rds(
      counts_marker,
      counts_marker_out
    )

    save_matrix_rds(
      cpm_marker,
      cpm_marker_out
    )

    save_matrix_rds(
      log2cpm_marker,
      log2cpm_marker_out
    )

    matched_signature <- prepare_matched_signature(
      signature_dt = signatures[[marker_set_i]],
      available_geneids = available_markers,
      output_file = signature_matched_out
    )

    fwrite(
      data.table(
        gene_id = missing_markers
      ),
      missing_marker_out,
      sep = "\t"
    )

    # --------------------------------------------------------------------------
    # Cobertura por tipo celular.
    # --------------------------------------------------------------------------

    coverage_by_type <- marker_coverage_by_type(
      marker_dt = marker_dt,
      bulk_geneids = bulk_geneids,
      dataset = dataset_i,
      marker_set = marker_set_i
    )

    marker_coverage_rows[[coverage_counter]] <- coverage_by_type
    coverage_counter <- coverage_counter + 1L

    marker_set_summary_rows[[marker_set_counter]] <- data.table(
      dataset = dataset_i,
      disease = disease_i,
      cohort_role = cohort_role_i,
      marker_set = marker_set_i,
      marker_assignments_total = nrow(
        marker_dt
      ),
      marker_unique_geneids_total = uniqueN(
        marker_dt$gene_id
      ),
      marker_unique_geneids_available = length(
        available_markers
      ),
      marker_unique_geneids_missing = length(
        missing_markers
      ),
      unique_geneid_coverage = length(
        available_markers
      ) /
        uniqueN(
          marker_dt$gene_id
        ),
      min_celltype_assignment_coverage = min(
        coverage_by_type$coverage
      ),
      median_celltype_assignment_coverage = median(
        coverage_by_type$coverage
      ),
      all_celltypes_at_least_80pct = all(
        coverage_by_type$coverage >=
          MARKER_COVERAGE_WARN_THRESHOLD
      )
    )

    marker_set_counter <- marker_set_counter + 1L

    marker_outputs[[marker_set_i]] <- list(
      counts = counts_marker_out,
      cpm = cpm_marker_out,
      log2cpm = log2cpm_marker_out,
      signature = signature_matched_out,
      n_genes = length(
        available_markers
      )
    )

    # Liberar intermediários.
    rm(
      counts_marker,
      cpm_marker,
      log2cpm_marker,
      matched_signature
    )

    invisible(
      gc()
    )
  }


  # ============================================================================
  # 6.12 Resumo do dataset
  # ============================================================================

  input_scope <- if (
    nrow(
      counts_full
    ) >=
      30000L
  ) {
    "broad_gene_universe"
  } else {
    "possibly_prefiltered_gene_universe"
  }

  dataset_summaries[[summary_counter]] <- data.table(
    dataset = dataset_i,
    disease = disease_i,
    cohort_role = cohort_role_i,
    count_input = count_file,
    metadata_input = metadata_file,
    count_input_geneid_source = counts_obj$gene_id_source,
    alignment_mode = alignment_mode_i,
    bulk_geneids_input = nrow(
      counts_full
    ),
    bulk_samples = ncol(
      counts_full
    ),
    controls = sum(
      metadata$group ==
        "Control"
    ),
    cases = sum(
      metadata$group ==
        disease_i
    ),
    reference_geneids = length(
      reference_geneids
    ),
    bulk_reference_intersection = length(
      common_reference_geneids
    ),
    reference_geneid_coverage = reference_coverage,
    library_size_min = min(
      library_sizes
    ),
    library_size_median = median(
      library_sizes
    ),
    library_size_max = max(
      library_sizes
    ),
    input_scope = input_scope
  )

  summary_counter <- summary_counter + 1L


  # ============================================================================
  # 6.13 Status do dataset
  # ============================================================================

  principal_coverage <- marker_set_summary_rows[[marker_set_counter - 2L]]

  # O índice acima corresponde ao top300 porque a ordem é top100, top300, top500.
  # Conferência defensiva:
  if (
    principal_coverage$marker_set !=
      PRIMARY_MARKER_SET
  ) {
    stop(
      dataset_i,
      ": erro interno ao recuperar resumo top300."
    )
  }

  dataset_status <- "PASS"
  dataset_note <- NA_character_

  if (
    principal_coverage$min_celltype_assignment_coverage <
      MARKER_COVERAGE_WARN_THRESHOLD
  ) {

    dataset_status <- "WARN"
    dataset_note <- paste0(
      "Top300 possui pelo menos um tipo celular com cobertura < ",
      MARKER_COVERAGE_WARN_THRESHOLD,
      ". Não houve preenchimento de genes ausentes com zero."
    )
  }

  if (
    input_scope ==
      "possibly_prefiltered_gene_universe"
  ) {

    dataset_status <- "WARN"

    extra_note <- paste0(
      "Matriz bulk possui <30000 GeneIDs e pode ter sido pré-filtrada. ",
      "Preferir universo raw/harmonizado amplo se disponível."
    )

    if (
      is.na(dataset_note)
    ) {
      dataset_note <- extra_note
    } else {
      dataset_note <- paste(
        dataset_note,
        extra_note
      )
    }
  }

  status_rows[[status_counter]] <- data.table(
    dataset = dataset_i,
    check = "dataset_prepared",
    observed = paste0(
      nrow(
        counts_full
      ),
      " genes; ",
      ncol(
        counts_full
      ),
      " samples"
    ),
    expected = "raw integer counts; Control + disease; GeneID intersection with Mathys",
    status = dataset_status,
    note = dataset_note
  )

  status_counter <- status_counter + 1L

  status_rows[[status_counter]] <- data.table(
    dataset = dataset_i,
    check = "counts_are_integer_nonnegative",
    observed = TRUE,
    expected = TRUE,
    status = "PASS",
    note = NA_character_
  )

  status_counter <- status_counter + 1L

  status_rows[[status_counter]] <- data.table(
    dataset = dataset_i,
    check = "reference_intersection_nonempty",
    observed = length(
      common_reference_geneids
    ),
    expected = ">0",
    status = if (
      length(
        common_reference_geneids
      ) > 0L
    ) {
      "PASS"
    } else {
      "FAIL"
    },
    note = NA_character_
  )

  status_counter <- status_counter + 1L

  status_rows[[status_counter]] <- data.table(
    dataset = dataset_i,
    check = "no_zero_filling",
    observed = TRUE,
    expected = TRUE,
    status = "PASS",
    note = "Somente GeneIDs realmente presentes no bulk foram mantidos."
  )

  status_counter <- status_counter + 1L


  # ============================================================================
  # 6.14 Manifesto de outputs
  # ============================================================================

  manifest_files <- c(
    metadata_file_out,
    library_file_out,
    intersection_file_out,
    ref_counts_out,
    ref_cpm_out,
    ref_log2cpm_out,
    marker_outputs$top100$counts,
    marker_outputs$top100$cpm,
    marker_outputs$top100$log2cpm,
    marker_outputs$top100$signature,
    marker_outputs$top300$counts,
    marker_outputs$top300$cpm,
    marker_outputs$top300$log2cpm,
    marker_outputs$top300$signature,
    marker_outputs$top500$counts,
    marker_outputs$top500$cpm,
    marker_outputs$top500$log2cpm,
    marker_outputs$top500$signature
  )

  manifest_objects <- c(
    "metadata",
    "library_sizes",
    "reference_intersection_geneids",
    "counts_reference_intersection",
    "cpm_reference_intersection",
    "log2cpm_reference_intersection",
    "counts_top100",
    "cpm_top100",
    "log2cpm_top100",
    "signature_top100",
    "counts_top300",
    "cpm_top300",
    "log2cpm_top300",
    "signature_top300",
    "counts_top500",
    "cpm_top500",
    "log2cpm_top500",
    "signature_top500"
  )

  output_manifest_rows[[manifest_counter]] <- data.table(
    dataset = dataset_i,
    object = manifest_objects,
    path = manifest_files,
    exists = file.exists(
      manifest_files
    ),
    size_MB = vapply(
      manifest_files,
      safe_file_size_mb,
      numeric(1)
    )
  )

  manifest_counter <- manifest_counter + 1L


  # ============================================================================
  # 6.15 Liberar memória antes do próximo dataset
  # ============================================================================

  rm(
    counts_obj,
    counts_full,
    counts_reference,
    cpm_reference,
    log2cpm_reference,
    metadata,
    library_sizes,
    library_dt
  )

  invisible(
    gc()
  )

  cat(
    dataset_i,
    " concluído.\n",
    sep = ""
  )
}


# ==============================================================================
# 7. CONSOLIDAR RESULTADOS
# ==============================================================================

dataset_summary <- rbindlist(
  dataset_summaries,
  use.names = TRUE,
  fill = TRUE
)

marker_coverage <- rbindlist(
  marker_coverage_rows,
  use.names = TRUE,
  fill = TRUE
)

marker_set_summary <- rbindlist(
  marker_set_summary_rows,
  use.names = TRUE,
  fill = TRUE
)

input_resolution <- rbindlist(
  input_resolution_rows,
  use.names = TRUE,
  fill = TRUE
)

status_table <- rbindlist(
  status_rows,
  use.names = TRUE,
  fill = TRUE
)

output_manifest <- rbindlist(
  output_manifest_rows,
  use.names = TRUE,
  fill = TRUE
)


# ==============================================================================
# 8. CHECAGENS GLOBAIS
# ==============================================================================

global_checks <- rbindlist(
  list(

    make_check(
      "four_datasets_prepared",
      uniqueN(
        dataset_summary$dataset
      ),
      4L,
      if (
        uniqueN(
          dataset_summary$dataset
        ) ==
          4L
      ) {
        "PASS"
      } else {
        "FAIL"
      }
    ),

    make_check(
      "ROSMAP_is_principal_AD",
      dataset_registry[
        dataset == "ROSMAP",
        cohort_role
      ],
      "principal",
      if (
        dataset_registry[
          dataset == "ROSMAP",
          cohort_role
        ] ==
          "principal"
      ) {
        "PASS"
      } else {
        "FAIL"
      }
    ),

    make_check(
      "GSE125583_is_validation_AD",
      dataset_registry[
        dataset == "GSE125583",
        cohort_role
      ],
      "validation",
      if (
        dataset_registry[
          dataset == "GSE125583",
          cohort_role
        ] ==
          "validation"
      ) {
        "PASS"
      } else {
        "FAIL"
      }
    ),

    make_check(
      "all_datasets_have_reference_intersection",
      all(
        dataset_summary$bulk_reference_intersection >
          0L
      ),
      TRUE,
      if (
        all(
          dataset_summary$bulk_reference_intersection >
            0L
        )
      ) {
        "PASS"
      } else {
        "FAIL"
      }
    ),

    make_check(
      "all_top300_have_available_markers",
      all(
        marker_set_summary[
          marker_set == "top300",
          marker_unique_geneids_available
        ] >
          0L
      ),
      TRUE,
      if (
        all(
          marker_set_summary[
            marker_set == "top300",
            marker_unique_geneids_available
          ] >
            0L
        )
      ) {
        "PASS"
      } else {
        "FAIL"
      }
    ),

    make_check(
      "all_outputs_exist",
      all(
        output_manifest$exists
      ),
      TRUE,
      if (
        all(
          output_manifest$exists
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

global_checks[
  ,
  dataset := "GLOBAL"
]

setcolorder(
  global_checks,
  c(
    "dataset",
    "check",
    "observed",
    "expected",
    "status",
    "note"
  )
)

status_table <- rbind(
  status_table,
  global_checks,
  fill = TRUE
)


# ==============================================================================
# 9. SALVAR QC CONSOLIDADO
# ==============================================================================

dataset_summary_file <- file.path(
  result_dir,
  "05_dataset_summary.tsv"
)

marker_coverage_file <- file.path(
  result_dir,
  "05_marker_coverage_by_dataset_celltype.tsv"
)

marker_set_summary_file <- file.path(
  result_dir,
  "05_marker_set_coverage.tsv"
)

input_resolution_file <- file.path(
  result_dir,
  "05_input_resolution.tsv"
)

status_file <- file.path(
  result_dir,
  "05_status.tsv"
)

manifest_file <- file.path(
  result_dir,
  "05_output_manifest.tsv"
)

fwrite(
  dataset_summary,
  dataset_summary_file,
  sep = "\t",
  na = "NA"
)

fwrite(
  marker_coverage,
  marker_coverage_file,
  sep = "\t",
  na = "NA"
)

fwrite(
  marker_set_summary,
  marker_set_summary_file,
  sep = "\t",
  na = "NA"
)

fwrite(
  input_resolution,
  input_resolution_file,
  sep = "\t",
  na = "NA"
)

fwrite(
  status_table,
  status_file,
  sep = "\t",
  na = "NA"
)

fwrite(
  output_manifest,
  manifest_file,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 10. LOG
# ==============================================================================

log_file <- file.path(
  result_dir,
  "05_prepare_bulk_deconvolution.txt"
)

sink(
  log_file
)

cat(
  "============================================================\n"
)

cat(
  "05_prepare_bulk_deconvolution.R\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Princípio de harmonização:\n"
)

cat(
  "  Cada bulk foi intersectado separadamente com Mathys por GeneID.\n"
)

cat(
  "  Nenhum universo universal AD-HD-PD-Mathys foi imposto.\n"
)

cat(
  "  Genes ausentes não foram preenchidos com zero.\n\n"
)

cat(
  "Escalas produzidas:\n"
)

cat(
  "  raw counts\n"
)

cat(
  "  CPM calculado com library size da matriz bulk de entrada\n"
)

cat(
  "  log2(CPM + 1)\n\n"
)

cat(
  "Datasets:\n"
)

print(
  dataset_summary
)

cat(
  "\nCobertura dos conjuntos de marcadores:\n"
)

print(
  marker_set_summary
)

cat(
  "\nStatus:\n"
)

print(
  status_table
)

sink()


# ==============================================================================
# 11. SESSION INFO
# ==============================================================================

session_file <- file.path(
  result_dir,
  "05_prepare_bulk_deconvolution_sessionInfo.txt"
)

sink(
  session_file
)

print(
  sessionInfo()
)

sink()


# ==============================================================================
# 12. VEREDITO
# ==============================================================================

if (
  any(
    status_table$status ==
      "FAIL",
    na.rm = TRUE
  )
) {
  stop(
    "Existe pelo menos um check FAIL em 05_status.tsv."
  )
}


# ==============================================================================
# 13. RESUMO NO TERMINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "05_prepare_bulk_deconvolution.R CONCLUÍDO\n"
)

cat(
  "============================================================\n\n"
)

print(
  dataset_summary[
    ,
    .(
      dataset,
      disease,
      cohort_role,
      bulk_geneids_input,
      bulk_samples,
      controls,
      cases,
      bulk_reference_intersection,
      reference_geneid_coverage,
      input_scope
    )
  ]
)

cat(
  "\nCobertura top300:\n"
)

print(
  marker_set_summary[
    marker_set ==
      "top300"
  ]
)

if (
  any(
    status_table$status ==
      "WARN",
    na.rm = TRUE
  )
) {

  cat(
    "\nATENÇÃO: existem WARNs. Revise:\n  ",
    status_file,
    "\n",
    sep = ""
  )

} else {

  cat(
    "\nTodos os datasets foram preparados sem WARN/FAIL.\n"
  )
}

cat(
  "\nOutputs:\n  ",
  bulk_output_root,
  "\n",
  sep = ""
)

cat(
  "QC:\n  ",
  result_dir,
  "\n",
  sep = ""
)

cat(
  "\nPróximas etapas:\n"
)

cat(
  "  06_run_MuSiC2.R\n"
)

cat(
  "  07_run_BisqueRNA.R\n"
)

cat(
  "  08_run_dtangle.R\n"
)
