#!/usr/bin/env Rscript

# ==============================================================================
# 00_audit_Mathys2019.R
# ==============================================================================
#
# FASE 02 — DECONVOLUÇÃO CELULAR
# ETAPA 00 — AUDITORIA DA REFERÊNCIA Mathys et al. 2019
#
# Objetivo:
#   Auditar a referência snRNA-seq antes de qualquer filtro, normalização,
#   harmonização, pseudobulk ou deconvolução.
#
# Este script responde:
#   1. matriz × genes × metadata estão estruturalmente compatíveis?
#   2. quantos genes, núcleos e doadores existem?
#   3. os counts são inteiros, não negativos e finitos?
#   4. os 8 tipos celulares esperados estão presentes?
#   5. qual a cobertura donor × cell type?
#   6. quais campos podem definir o status neuropatológico?
#   7. existe algum candidato inequívoco para os 24 pathology-free?
#   8. é possível avaliar overlap entre doadores Mathys e bulk ROSMAP?
#
# IMPORTANTE:
#   Este script NÃO seleciona os 24 pathology-free.
#   Ele NÃO remove genes/células, NÃO normaliza, NÃO cria pseudobulk,
#   NÃO harmoniza para GeneID e NÃO roda deconvolução.
#
# A seleção pathology-free será feita somente no:
#   01_prepare_metadata_Mathys2019.R
#
# PRIVACIDADE:
#   Outputs agregados não exportam individualID/projid originais.
# ==============================================================================


# ==============================================================================
# 0. PACOTES
# ==============================================================================

required_packages <- c("Matrix", "data.table")

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
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})


# ==============================================================================
# 1. PARÂMETROS DA AUDITORIA
# ==============================================================================

EXPECTED_DONORS <- 48L
EXPECTED_PATHOLOGY_FREE_DONORS <- 24L

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

FRACTION_TOLERANCE <- 1e-8


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

ref_dir <- file.path(
  "data",
  "reference",
  "mathys2019",
  "raw"
)

matrix_file <- file.path(
  ref_dir,
  "filtered_count_matrix.mtx"
)

genes_file <- file.path(
  ref_dir,
  "filtered_gene_row_names.txt"
)

cells_file <- file.path(
  ref_dir,
  "filtered_column_metadata.txt"
)

# Metadados auxiliares: opcionais no 00.
# Se presentes, ajudam a resolver doadores e status neuropatológico.
optional_metadata_files <- c(
  sample_key = file.path(
    ref_dir,
    "snRNAseqPFC_BA10_Sample_key.csv"
  ),
  assay = file.path(
    ref_dir,
    "snRNAseqPFC_BA10_assay_scRNAseq_metadata.csv"
  ),
  biospecimen = file.path(
    ref_dir,
    "snRNAseqPFC_BA10_biospecimen_metadata.csv"
  ),
  id_mapping = file.path(
    ref_dir,
    "snRNAseqPFC_BA10_id_mapping.csv"
  )
)

# Metadata bulk ROSMAP da Fase 01c.
# É usado somente para tentar estimar overlap de doadores.
rosmap_bulk_metadata_file <- file.path(
  "results",
  "candidate_ROSMAP",
  "metadata",
  "metadata_ROSMAP.tsv"
)

# Manifesto privado/local opcional.
# Se contiver IDs/versões Synapse, eles não serão copiados para o relatório
# agregado deste script.
private_synapse_manifest <- file.path(
  ref_dir,
  "synapse_manifest.tsv"
)

out_dir <- file.path(
  "results",
  "deconvolution",
  "mathys2019",
  "audit"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

log_file <- file.path(
  out_dir,
  "00_audit_Mathys2019.txt"
)

session_file <- file.path(
  out_dir,
  "00_audit_Mathys2019_sessionInfo.txt"
)


# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================

fmt_int <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }

  format(
    as.integer(round(as.numeric(x))),
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


fmt_num <- function(x, digits = 4L) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }

  format(
    round(as.numeric(x), digits),
    nsmall = digits,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


normalize_character <- function(x) {
  x <- trimws(as.character(x))

  x[
    is.na(x) |
      x == "" |
      toupper(x) == "NA"
  ] <- NA_character_

  x
}


detect_donor_column <- function(dt) {
  candidates <- c(
    "individualID",
    "individual_id",
    "projid",
    "donor_id",
    "donor",
    "subject_id",
    "participant_id",
    "subject",
    "participant"
  )

  hit <- candidates[
    candidates %in% names(dt)
  ]

  if (length(hit) == 0L) {
    return(NA_character_)
  }

  hit[1]
}


find_pathology_candidate_columns <- function(dt) {
  pattern <- paste(
    c(
      "patholog",
      "diagn",
      "cogdx",
      "cogn",
      "braak",
      "cerad",
      "disease",
      "dement",
      "amyloid",
      "plaque",
      "case",
      "control",
      "clinical",
      "status"
    ),
    collapse = "|"
  )

  names(dt)[
    grepl(
      pattern,
      names(dt),
      ignore.case = TRUE,
      perl = TRUE
    )
  ]
}


looks_pathology_free <- function(x) {
  x <- normalize_character(x)

  ifelse(
    is.na(x),
    FALSE,
    grepl(
      paste(
        c(
          "no[ _-]?path",
          "pathology[ _-]?free",
          "without[ _-]?path",
          "no[ _-]?ad",
          "^control$",
          "^normal$",
          "unaffected"
        ),
        collapse = "|"
      ),
      x,
      ignore.case = TRUE,
      perl = TRUE
    )
  )
}


classify_gene_id <- function(x) {
  x <- normalize_character(x)

  fifelse(
    is.na(x),
    "missing",
    fifelse(
      grepl(
        "^ENSG[0-9]+(\\.[0-9]+)?(_PAR_Y)?$",
        x
      ),
      "Ensembl_gene",
      fifelse(
        grepl(
          "^[0-9]+$",
          x
        ),
        "numeric_possible_Entrez",
        "other_or_symbol"
      )
    )
  )
}


metadata_column_inventory <- function(dt, source_name) {
  rbindlist(
    lapply(
      names(dt),
      function(col) {
        x <- dt[[col]]

        data.table(
          source = source_name,
          column = col,
          class = paste(class(x), collapse = ";"),
          n_rows = length(x),
          n_missing = sum(is.na(x)),
          n_unique_nonmissing = uniqueN(
            x[!is.na(x)]
          )
        )
      }
    ),
    fill = TRUE
  )
}


pathology_candidate_summary <- function(dt, source_name) {
  candidate_cols <- find_pathology_candidate_columns(dt)

  if (length(candidate_cols) == 0L) {
    return(
      list(
        columns = data.table(),
        values = data.table()
      )
    )
  }

  donor_col <- detect_donor_column(dt)

  column_summary <- rbindlist(
    lapply(
      candidate_cols,
      function(col) {
        x <- dt[[col]]

        data.table(
          source = source_name,
          column = col,
          class = paste(class(x), collapse = ";"),
          donor_column_used = donor_col,
          n_rows = length(x),
          n_missing = sum(is.na(x)),
          n_unique_nonmissing = uniqueN(
            x[!is.na(x)]
          )
        )
      }
    ),
    fill = TRUE
  )

  value_tables <- list()

  for (col in candidate_cols) {
    x <- normalize_character(dt[[col]])

    n_unique <- uniqueN(
      x[!is.na(x)]
    )

    # Evitar exportar variáveis contínuas ou campos excessivamente granulares.
    if (n_unique == 0L || n_unique > 50L) {
      next
    }

    tmp <- data.table(
      value = x
    )[
      !is.na(value),
      .(
        n_rows = .N
      ),
      by = value
    ]

    if (!is.na(donor_col)) {
      donor_values <- normalize_character(
        dt[[donor_col]]
      )

      donor_tmp <- data.table(
        value = x,
        donor = donor_values
      )[
        !is.na(value) &
          !is.na(donor),
        .(
          n_donors = uniqueN(donor)
        ),
        by = value
      ]

      tmp <- merge(
        tmp,
        donor_tmp,
        by = "value",
        all.x = TRUE,
        sort = FALSE
      )

    } else {
      tmp[
        ,
        n_donors := NA_integer_
      ]
    }

    tmp[
      ,
      `:=`(
        source = source_name,
        column = col,
        donor_column_used = donor_col,
        looks_pathology_free = looks_pathology_free(value),
        exactly_24_donors = !is.na(n_donors) &
          n_donors == EXPECTED_PATHOLOGY_FREE_DONORS
      )
    ]

    setcolorder(
      tmp,
      c(
        "source",
        "column",
        "donor_column_used",
        "value",
        "n_rows",
        "n_donors",
        "looks_pathology_free",
        "exactly_24_donors"
      )
    )

    value_tables[[length(value_tables) + 1L]] <- tmp
  }

  values <- if (length(value_tables) > 0L) {
    rbindlist(
      value_tables,
      fill = TRUE
    )
  } else {
    data.table()
  }

  list(
    columns = column_summary,
    values = values
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


# ==============================================================================
# 4. AUDITORIA PRINCIPAL
# ==============================================================================

run_audit <- function() {

  cat(
    "============================================================\n"
  )
  cat(
    "FASE 02 — 00_audit_Mathys2019.R\n"
  )
  cat(
    "AUDITORIA DA REFERÊNCIA snRNA-seq MATHYS 2019\n"
  )
  cat(
    "============================================================\n\n"
  )


  # ============================================================================
  # 4.1 ARQUIVOS
  # ============================================================================

  cat(
    "------------------------------------------------------------\n"
  )
  cat(
    "1. ARQUIVOS\n"
  )
  cat(
    "------------------------------------------------------------\n"
  )

  required_files <- c(
    matrix = matrix_file,
    genes = genes_file,
    cells = cells_file
  )

  missing_required <- required_files[
    !file.exists(required_files)
  ]

  if (length(missing_required) > 0L) {
    cat(
      "\nArquivos obrigatórios ausentes:\n"
    )
    print(missing_required)

    stop(
      "Auditoria interrompida: arquivos principais ausentes."
    )
  }

  all_files <- c(
    required_files,
    optional_metadata_files
  )

  file_audit <- data.table(
    role = names(all_files),
    file_name = basename(all_files),
    required = names(all_files) %in% names(required_files),
    exists = file.exists(all_files),
    size_MB = fifelse(
      file.exists(all_files),
      round(
        file.info(all_files)$size / 1024^2,
        3
      ),
      NA_real_
    )
  )

  fwrite(
    file_audit,
    file.path(
      out_dir,
      "00_file_audit.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(file_audit)

  cat(
    "\nManifesto Synapse privado/local presente: ",
    file.exists(private_synapse_manifest),
    "\n",
    sep = ""
  )

  cat(
    "IDs/versões Synapse não são copiados para os outputs agregados.\n\n"
  )


  # ============================================================================
  # 4.2 MATRIZ DE CONTAGENS
  # ============================================================================

  cat(
    "------------------------------------------------------------\n"
  )
  cat(
    "2. MATRIZ DE CONTAGENS\n"
  )
  cat(
    "------------------------------------------------------------\n"
  )

  cat(
    "Lendo matriz Matrix Market...\n"
  )

  counts <- Matrix::readMM(
    matrix_file
  )

  counts <- as(
    counts,
    "CsparseMatrix"
  )

  nonzero_values <- counts@x

  matrix_na <- sum(
    is.na(nonzero_values)
  )

  matrix_nonfinite <- sum(
    !is.finite(nonzero_values)
  )

  matrix_negative <- sum(
    nonzero_values < 0,
    na.rm = TRUE
  )

  matrix_fractional <- sum(
    abs(
      nonzero_values -
        round(nonzero_values)
    ) > FRACTION_TOLERANCE,
    na.rm = TRUE
  )

  gene_sums <- Matrix::rowSums(counts)
  cell_sums <- Matrix::colSums(counts)

  all_zero_genes <- sum(
    gene_sums == 0
  )

  all_zero_cells <- sum(
    cell_sums == 0
  )

  total_entries <- as.double(
    nrow(counts)
  ) *
    as.double(
      ncol(counts)
    )

  nnz <- length(
    nonzero_values
  )

  zero_fraction <- 1 -
    (
      nnz /
        total_entries
    )

  matrix_audit <- data.table(
    metric = c(
      "class",
      "genes_rows",
      "nuclei_columns",
      "nonzero_entries",
      "zero_fraction",
      "NA_nonzero_values",
      "nonfinite_nonzero_values",
      "negative_nonzero_values",
      "fractional_nonzero_values",
      "all_zero_gene_rows",
      "all_zero_nucleus_columns",
      "minimum_nonzero_value",
      "maximum_value"
    ),
    value = c(
      class(counts)[1],
      nrow(counts),
      ncol(counts),
      nnz,
      zero_fraction,
      matrix_na,
      matrix_nonfinite,
      matrix_negative,
      matrix_fractional,
      all_zero_genes,
      all_zero_cells,
      if (length(nonzero_values) > 0L) {
        min(
          nonzero_values,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
      if (length(nonzero_values) > 0L) {
        max(
          nonzero_values,
          na.rm = TRUE
        )
      } else {
        NA_real_
      }
    )
  )

  fwrite(
    matrix_audit,
    file.path(
      out_dir,
      "00_matrix_audit.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  cat(
    "Genes na matriz: ",
    fmt_int(nrow(counts)),
    "\n",
    sep = ""
  )

  cat(
    "Núcleos na matriz: ",
    fmt_int(ncol(counts)),
    "\n",
    sep = ""
  )

  cat(
    "Zeros: ",
    fmt_num(
      100 * zero_fraction,
      2
    ),
    "%\n",
    sep = ""
  )

  cat(
    "NA: ",
    matrix_na,
    "\n",
    sep = ""
  )

  cat(
    "Não finitos: ",
    matrix_nonfinite,
    "\n",
    sep = ""
  )

  cat(
    "Negativos: ",
    matrix_negative,
    "\n",
    sep = ""
  )

  cat(
    "Fracionários: ",
    matrix_fractional,
    "\n",
    sep = ""
  )

  cat(
    "Genes all-zero: ",
    fmt_int(all_zero_genes),
    "\n",
    sep = ""
  )

  cat(
    "Núcleos all-zero: ",
    fmt_int(all_zero_cells),
    "\n\n",
    sep = ""
  )


  # ============================================================================
  # 4.3 ARQUIVO DE GENES
  # ============================================================================

  cat(
    "------------------------------------------------------------\n"
  )
  cat(
    "3. ARQUIVO DE GENES\n"
  )
  cat(
    "------------------------------------------------------------\n"
  )

  genes <- fread(
    genes_file,
    header = FALSE,
    fill = TRUE,
    na.strings = c(
      "NA",
      ""
    )
  )

  if (ncol(genes) == 0L) {
    stop(
      "Arquivo de genes sem colunas."
    )
  }

  setnames(
    genes,
    paste0(
      "gene_col_",
      seq_len(ncol(genes))
    )
  )

  primary_gene_id <- normalize_character(
    genes$gene_col_1
  )

  gene_format <- classify_gene_id(
    primary_gene_id
  )

  gene_identifier_audit <- data.table(
    identifier_format = gene_format
  )[
    ,
    .N,
    by = identifier_format
  ][
    order(-N)
  ]

  fwrite(
    gene_identifier_audit,
    file.path(
      out_dir,
      "00_gene_identifier_audit.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  gene_summary <- data.table(
    metric = c(
      "rows",
      "columns",
      "missing_primary_id",
      "unique_primary_id",
      "duplicated_primary_id_rows",
      "matrix_rows_equal_gene_rows"
    ),
    value = c(
      nrow(genes),
      ncol(genes),
      sum(
        is.na(primary_gene_id)
      ),
      uniqueN(
        primary_gene_id[
          !is.na(primary_gene_id)
        ]
      ),
      sum(
        duplicated(
          primary_gene_id[
            !is.na(primary_gene_id)
          ]
        )
      ),
      nrow(counts) == nrow(genes)
    )
  )

  fwrite(
    gene_summary,
    file.path(
      out_dir,
      "00_gene_file_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  cat(
    "Linhas no arquivo de genes: ",
    fmt_int(nrow(genes)),
    "\n",
    sep = ""
  )

  cat(
    "Colunas no arquivo de genes: ",
    ncol(genes),
    "\n",
    sep = ""
  )

  cat(
    "Linhas matriz == linhas genes: ",
    nrow(counts) == nrow(genes),
    "\n\n",
    sep = ""
  )

  cat(
    "Formatos observados na primeira coluna:\n"
  )

  print(
    gene_identifier_audit
  )

  if (
    nrow(gene_identifier_audit) == 1L &&
      gene_identifier_audit$identifier_format[1] == "other_or_symbol"
  ) {
    cat(
      "\nObservação: 'other_or_symbol' significa apenas que os IDs não têm ",
      "formato Ensembl nem são puramente numéricos. ",
      "A natureza exata do identificador será confirmada antes da harmonização ",
      "para GeneID.\n"
    )
  }


  # ============================================================================
  # 4.4 METADATA DOS NÚCLEOS
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )
  cat(
    "4. METADATA DOS NÚCLEOS\n"
  )
  cat(
    "------------------------------------------------------------\n"
  )

  cell_meta <- fread(
    cells_file,
    na.strings = c(
      "NA",
      ""
    )
  )

  cell_inventory <- metadata_column_inventory(
    cell_meta,
    "cell_metadata"
  )

  fwrite(
    cell_inventory,
    file.path(
      out_dir,
      "00_cell_metadata_column_inventory.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  cells_match <- ncol(counts) == nrow(cell_meta)

  cat(
    "Linhas do metadata celular: ",
    fmt_int(nrow(cell_meta)),
    "\n",
    sep = ""
  )

  cat(
    "Colunas do metadata celular: ",
    ncol(cell_meta),
    "\n",
    sep = ""
  )

  cat(
    "Colunas matriz == linhas metadata: ",
    cells_match,
    "\n",
    sep = ""
  )

  cat(
    "Colunas disponíveis:\n"
  )

  print(
    names(cell_meta)
  )

  required_cell_columns <- c(
    "TAG",
    "projid",
    "broad.cell.type"
  )

  missing_cell_columns <- setdiff(
    required_cell_columns,
    names(cell_meta)
  )

  if (length(missing_cell_columns) > 0L) {
    stop(
      paste0(
        "Metadata celular sem colunas obrigatórias: ",
        paste(
          missing_cell_columns,
          collapse = ", "
        )
      )
    )
  }

  tag <- normalize_character(
    cell_meta$TAG
  )

  projid <- normalize_character(
    cell_meta$projid
  )

  broad_type <- normalize_character(
    cell_meta$broad.cell.type
  )

  n_tag_missing <- sum(
    is.na(tag)
  )

  n_tag_duplicated <- sum(
    duplicated(
      tag[
        !is.na(tag)
      ]
    )
  )

  n_projid_missing <- sum(
    is.na(projid)
  )

  n_donors_projid <- uniqueN(
    projid[
      !is.na(projid)
    ]
  )

  n_celltype_missing <- sum(
    is.na(broad_type)
  )

  observed_cell_types <- sort(
    unique(
      broad_type[
        !is.na(broad_type)
      ]
    )
  )

  missing_expected_cell_types <- setdiff(
    EXPECTED_CELL_TYPES,
    observed_cell_types
  )

  unexpected_cell_types <- setdiff(
    observed_cell_types,
    EXPECTED_CELL_TYPES
  )

  cell_metadata_summary <- data.table(
    metric = c(
      "rows",
      "columns",
      "matrix_columns_equal_metadata_rows",
      "TAG_missing",
      "TAG_duplicated_rows",
      "projid_missing",
      "unique_projid",
      "broad_cell_type_missing",
      "observed_cell_types",
      "expected_cell_types_missing",
      "unexpected_cell_types"
    ),
    value = c(
      nrow(cell_meta),
      ncol(cell_meta),
      cells_match,
      n_tag_missing,
      n_tag_duplicated,
      n_projid_missing,
      n_donors_projid,
      n_celltype_missing,
      paste(
        observed_cell_types,
        collapse = ";"
      ),
      if (length(missing_expected_cell_types) == 0L) {
        "none"
      } else {
        paste(
          missing_expected_cell_types,
          collapse = ";"
        )
      },
      if (length(unexpected_cell_types) == 0L) {
        "none"
      } else {
        paste(
          unexpected_cell_types,
          collapse = ";"
        )
      }
    )
  )

  fwrite(
    cell_metadata_summary,
    file.path(
      out_dir,
      "00_cell_metadata_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  cat(
    "\nDoadores por projid: ",
    n_donors_projid,
    "\n",
    sep = ""
  )

  cat(
    "Esperado pelo desenho: ",
    EXPECTED_DONORS,
    "\n",
    sep = ""
  )

  cat(
    "Tipos celulares observados: ",
    paste(
      observed_cell_types,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )

  if (length(missing_expected_cell_types) > 0L) {
    cat(
      "ATENÇÃO — tipos esperados ausentes: ",
      paste(
        missing_expected_cell_types,
        collapse = ", "
      ),
      "\n",
      sep = ""
    )
  }

  if (length(unexpected_cell_types) > 0L) {
    cat(
      "ATENÇÃO — tipos não previstos: ",
      paste(
        unexpected_cell_types,
        collapse = ", "
      ),
      "\n",
      sep = ""
    )
  }


  # ============================================================================
  # 4.5 COBERTURA DONOR × CELL TYPE
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )
  cat(
    "5. COBERTURA DONOR × CELL TYPE\n"
  )
  cat(
    "------------------------------------------------------------\n"
  )

  donor_cell_observed <- data.table(
    projid = projid,
    cell_type = broad_type
  )[
    !is.na(projid) &
      !is.na(cell_type),
    .(
      n_nuclei = .N
    ),
    by = .(
      projid,
      cell_type
    )
  ]

  all_combinations <- CJ(
    projid = sort(
      unique(
        projid[
          !is.na(projid)
        ]
      )
    ),
    cell_type = sort(
      unique(
        broad_type[
          !is.na(broad_type)
        ]
      )
    ),
    unique = TRUE
  )

  donor_cell_complete <- merge(
    all_combinations,
    donor_cell_observed,
    by = c(
      "projid",
      "cell_type"
    ),
    all.x = TRUE,
    sort = TRUE
  )

  donor_cell_complete[
    is.na(n_nuclei),
    n_nuclei := 0L
  ]

  # Não exportar projid individuais; salvar somente resumo por tipo.
  # IMPORTANTE:
  # data.table exige que uma mesma coluna produzida por grupos tenha exatamente
  # o mesmo tipo em todos os grupos. Dependendo do número de observações,
  # median() sobre um vetor integer pode retornar internamente integer em um
  # grupo e double em outro. Isso pode interromper a agregação com a mensagem:
  #
  #   "Column ... is type 'integer' but expecting type 'double'"
  #
  # Por isso, as estatísticas quantitativas abaixo são convertidas
  # explicitamente para numeric/double.
  celltype_summary <- donor_cell_complete[
    ,
    .(
      donors_total = as.integer(
        uniqueN(projid)
      ),
      donors_present = as.integer(
        sum(
          n_nuclei > 0
        )
      ),
      donors_absent = as.integer(
        sum(
          n_nuclei == 0
        )
      ),
      total_nuclei = as.numeric(
        sum(n_nuclei)
      ),
      min_nuclei_including_zero = as.numeric(
        min(n_nuclei)
      ),
      median_nuclei_including_zero = as.numeric(
        median(n_nuclei)
      ),
      max_nuclei = as.numeric(
        max(n_nuclei)
      ),
      min_nuclei_when_present = as.numeric(
        if (any(n_nuclei > 0)) {
          min(
            n_nuclei[
              n_nuclei > 0
            ]
          )
        } else {
          0
        }
      ),
      median_nuclei_when_present = as.numeric(
        if (any(n_nuclei > 0)) {
          median(
            n_nuclei[
              n_nuclei > 0
            ]
          )
        } else {
          0
        }
      )
    ),
    by = cell_type
  ][
    order(-total_nuclei)
  ]

  fwrite(
    celltype_summary,
    file.path(
      out_dir,
      "00_celltype_donor_coverage.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  cat(
    "Cobertura agregada por tipo celular:\n"
  )

  print(
    celltype_summary
  )

  zero_donor_cell_combinations <- sum(
    donor_cell_complete$n_nuclei == 0
  )

  cat(
    "\nCombinações donor × cell type sem núcleos: ",
    fmt_int(zero_donor_cell_combinations),
    " / ",
    fmt_int(nrow(donor_cell_complete)),
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.6 METADADOS AUXILIARES
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )
  cat(
    "6. METADADOS AUXILIARES\n"
  )
  cat(
    "------------------------------------------------------------\n"
  )

  metadata_objects <- list(
    cell_metadata = cell_meta
  )

  aux_inventory_list <- list()

  all_column_inventory <- list(
    metadata_column_inventory(
      cell_meta,
      "cell_metadata"
    )
  )

  for (nm in names(optional_metadata_files)) {

    f <- optional_metadata_files[[nm]]

    if (!file.exists(f)) {

      aux_inventory_list[[length(aux_inventory_list) + 1L]] <- data.table(
        source = nm,
        file_name = basename(f),
        exists = FALSE,
        n_rows = NA_integer_,
        n_columns = NA_integer_,
        donor_column_detected = NA_character_
      )

      cat(
        nm,
        ": AUSENTE — ",
        basename(f),
        "\n",
        sep = ""
      )

      next
    }

    dt <- fread(
      f,
      na.strings = c(
        "NA",
        ""
      )
    )

    metadata_objects[[nm]] <- dt

    donor_col <- detect_donor_column(dt)

    aux_inventory_list[[length(aux_inventory_list) + 1L]] <- data.table(
      source = nm,
      file_name = basename(f),
      exists = TRUE,
      n_rows = nrow(dt),
      n_columns = ncol(dt),
      donor_column_detected = donor_col
    )

    all_column_inventory[[length(all_column_inventory) + 1L]] <- metadata_column_inventory(
      dt,
      nm
    )

    cat(
      nm,
      ": ",
      nrow(dt),
      " linhas × ",
      ncol(dt),
      " colunas\n",
      sep = ""
    )

    cat(
      "  colunas: ",
      paste(
        names(dt),
        collapse = ", "
      ),
      "\n",
      sep = ""
    )
  }

  aux_inventory <- rbindlist(
    aux_inventory_list,
    fill = TRUE
  )

  fwrite(
    aux_inventory,
    file.path(
      out_dir,
      "00_aux_metadata_inventory.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  full_column_inventory <- rbindlist(
    all_column_inventory,
    fill = TRUE
  )

  fwrite(
    full_column_inventory,
    file.path(
      out_dir,
      "00_all_metadata_column_inventory.tsv"
    ),
    sep = "\t",
    na = "NA"
  )


  # ============================================================================
  # 4.7 CANDIDATOS A PATHOLOGY-FREE
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )
  cat(
    "7. CANDIDATOS PARA DEFINIR PATHOLOGY-FREE\n"
  )
  cat(
    "------------------------------------------------------------\n"
  )

  pathology_columns_list <- list()
  pathology_values_list <- list()

  for (nm in names(metadata_objects)) {

    res <- pathology_candidate_summary(
      metadata_objects[[nm]],
      nm
    )

    if (nrow(res$columns) > 0L) {
      pathology_columns_list[[length(pathology_columns_list) + 1L]] <- res$columns
    }

    if (nrow(res$values) > 0L) {
      pathology_values_list[[length(pathology_values_list) + 1L]] <- res$values
    }
  }

  pathology_columns <- if (
    length(pathology_columns_list) > 0L
  ) {
    rbindlist(
      pathology_columns_list,
      fill = TRUE
    )
  } else {
    data.table()
  }

  pathology_values <- if (
    length(pathology_values_list) > 0L
  ) {
    rbindlist(
      pathology_values_list,
      fill = TRUE
    )
  } else {
    data.table()
  }

  if (nrow(pathology_columns) > 0L) {
    fwrite(
      pathology_columns,
      file.path(
        out_dir,
        "00_pathology_candidate_columns.tsv"
      ),
      sep = "\t",
      na = "NA"
    )
  } else {
    fwrite(
      data.table(
        note = "Nenhuma coluna candidata identificada automaticamente."
      ),
      file.path(
        out_dir,
        "00_pathology_candidate_columns.tsv"
      ),
      sep = "\t"
    )
  }

  if (nrow(pathology_values) > 0L) {
    fwrite(
      pathology_values,
      file.path(
        out_dir,
        "00_pathology_candidate_value_counts.tsv"
      ),
      sep = "\t",
      na = "NA"
    )
  } else {
    fwrite(
      data.table(
        note = "Nenhuma categoria de patologia pôde ser resumida automaticamente."
      ),
      file.path(
        out_dir,
        "00_pathology_candidate_value_counts.tsv"
      ),
      sep = "\t"
    )
  }

  if (nrow(pathology_columns) == 0L) {

    cat(
      "Nenhuma coluna candidata a patologia foi localizada automaticamente.\n"
    )

  } else {

    cat(
      "Colunas candidatas encontradas:\n"
    )

    print(
      pathology_columns
    )
  }

  exact_24_candidates <- if (
    nrow(pathology_values) > 0L
  ) {
    pathology_values[
      exactly_24_donors == TRUE |
        (
          looks_pathology_free == TRUE &
            !is.na(n_donors)
        )
    ]
  } else {
    data.table()
  }

  if (nrow(exact_24_candidates) > 0L) {

    cat(
      "\nCategorias potencialmente informativas para os 24 pathology-free:\n"
    )

    print(
      exact_24_candidates
    )

  } else {

    cat(
      "\nNenhuma categoria foi automaticamente confirmada como os 24 pathology-free.\n"
    )

    cat(
      "Isso NÃO é erro; a seleção será definida após inspeção do metadata.\n"
    )
  }


  # ============================================================================
  # 4.8 MAPA projid -> individualID
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )
  cat(
    "8. MAPA DE DOADORES\n"
  )
  cat(
    "------------------------------------------------------------\n"
  )

  donor_map <- NULL
  donor_map_source <- NA_character_
  donor_map_ambiguous <- NA_integer_

  for (nm in names(metadata_objects)) {

    dt <- metadata_objects[[nm]]

    if (
      all(
        c(
          "projid",
          "individualID"
        ) %in% names(dt)
      )
    ) {

      candidate_map <- unique(
        dt[
          !is.na(projid) &
            !is.na(individualID),
          .(
            projid = normalize_character(projid),
            individualID = normalize_character(individualID)
          )
        ]
      )

      ambiguity <- candidate_map[
        ,
        .(
          n_individualID = uniqueN(individualID)
        ),
        by = projid
      ][
        n_individualID > 1L
      ]

      if (nrow(ambiguity) == 0L) {
        donor_map <- candidate_map
        donor_map_source <- nm
        donor_map_ambiguous <- 0L
        break
      }

      donor_map_ambiguous <- nrow(ambiguity)
    }
  }

  if (!is.null(donor_map)) {

    mathys_projids <- sort(
      unique(
        projid[
          !is.na(projid)
        ]
      )
    )

    idx <- match(
      mathys_projids,
      donor_map$projid
    )

    mathys_individual_ids <- donor_map$individualID[
      idx
    ]

    n_missing_individual_map <- sum(
      is.na(mathys_individual_ids)
    )

    cat(
      "Mapa projid -> individualID encontrado em: ",
      donor_map_source,
      "\n",
      sep = ""
    )

    cat(
      "Doadores da matriz sem individualID após o mapa: ",
      n_missing_individual_map,
      "\n",
      sep = ""
    )

  } else {

    mathys_individual_ids <- NULL
    n_missing_individual_map <- NA_integer_

    cat(
      "Mapa projid -> individualID não pôde ser estabelecido nesta execução.\n"
    )
  }

  donor_map_summary <- data.table(
    metric = c(
      "map_available",
      "source",
      "ambiguous_projid_detected",
      "matrix_unique_projid",
      "matrix_projid_without_individualID"
    ),
    value = c(
      !is.null(donor_map),
      donor_map_source,
      donor_map_ambiguous,
      n_donors_projid,
      n_missing_individual_map
    )
  )

  fwrite(
    donor_map_summary,
    file.path(
      out_dir,
      "00_donor_mapping_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )


  # ============================================================================
  # 4.9 OVERLAP MATHYS × BULK ROSMAP
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )
  cat(
    "9. OVERLAP DE DOADORES MATHYS × BULK ROSMAP\n"
  )
  cat(
    "------------------------------------------------------------\n"
  )

  overlap_summary <- data.table(
    overlap_assessed = FALSE,
    mathys_identifier_type = NA_character_,
    bulk_identifier_column = NA_character_,
    mathys_unique_donors = NA_integer_,
    bulk_unique_donors = NA_integer_,
    overlap_n = NA_integer_,
    note = NA_character_
  )

  if (!file.exists(rosmap_bulk_metadata_file)) {

    overlap_summary[
      ,
      note :=
        "Metadata bulk ROSMAP não encontrado; overlap não avaliado."
    ]

    cat(
      "Metadata bulk ROSMAP não encontrado. Overlap não avaliado.\n"
    )

  } else {

    bulk_meta <- fread(
      rosmap_bulk_metadata_file,
      na.strings = c(
        "NA",
        ""
      )
    )

    bulk_donor_col <- detect_donor_column(
      bulk_meta
    )

    if (is.na(bulk_donor_col)) {

      overlap_summary[
        ,
        note :=
          "Nenhuma coluna de doador reconhecida no metadata bulk ROSMAP."
      ]

      cat(
        "Nenhuma coluna de doador reconhecida no metadata bulk ROSMAP.\n"
      )

    } else {

      bulk_donors <- normalize_character(
        bulk_meta[[bulk_donor_col]]
      )

      bulk_donors <- sort(
        unique(
          bulk_donors[
            !is.na(bulk_donors)
          ]
        )
      )

      # Comparar individualID quando ambos os lados usam esse identificador.
      if (
        !is.null(mathys_individual_ids) &&
          bulk_donor_col %in% c(
            "individualID",
            "individual_id"
          )
      ) {

        mathys_ids <- sort(
          unique(
            normalize_character(
              mathys_individual_ids
            )
          )
        )

        mathys_ids <- mathys_ids[
          !is.na(mathys_ids)
        ]

        overlap_ids <- intersect(
          mathys_ids,
          bulk_donors
        )

        overlap_summary <- data.table(
          overlap_assessed = TRUE,
          mathys_identifier_type = "individualID",
          bulk_identifier_column = bulk_donor_col,
          mathys_unique_donors = length(mathys_ids),
          bulk_unique_donors = length(bulk_donors),
          overlap_n = length(overlap_ids),
          note = paste0(
            "Somente a contagem de overlap foi exportada; ",
            "IDs individuais não foram salvos."
          )
        )

      } else if (
        bulk_donor_col == "projid"
      ) {

        mathys_ids <- sort(
          unique(
            projid[
              !is.na(projid)
            ]
          )
        )

        overlap_ids <- intersect(
          mathys_ids,
          bulk_donors
        )

        overlap_summary <- data.table(
          overlap_assessed = TRUE,
          mathys_identifier_type = "projid",
          bulk_identifier_column = bulk_donor_col,
          mathys_unique_donors = length(mathys_ids),
          bulk_unique_donors = length(bulk_donors),
          overlap_n = length(overlap_ids),
          note = paste0(
            "Somente a contagem de overlap foi exportada; ",
            "IDs individuais não foram salvos."
          )
        )

      } else {

        overlap_summary[
          ,
          `:=`(
            mathys_identifier_type = if (
              !is.null(mathys_individual_ids)
            ) {
              "individualID available"
            } else {
              "projid only"
            },
            bulk_identifier_column = bulk_donor_col,
            note = paste0(
              "Tipos de identificador não puderam ser comparados com segurança; ",
              "nenhum overlap foi inferido."
            )
          )
        ]
      }
    }
  }

  fwrite(
    overlap_summary,
    file.path(
      out_dir,
      "00_donor_overlap_ROSMAP_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(
    overlap_summary
  )


  # ============================================================================
  # 4.10 STATUS GLOBAL
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )
  cat(
    "10. STATUS GLOBAL\n"
  )
  cat(
    "------------------------------------------------------------\n"
  )

  pathology_24_found <- if (
    nrow(pathology_values) > 0L
  ) {
    any(
      pathology_values$exactly_24_donors == TRUE &
        pathology_values$looks_pathology_free == TRUE,
      na.rm = TRUE
    )
  } else {
    FALSE
  }

  checks <- rbindlist(
    list(
      make_check(
        "matrix_rows_equal_gene_rows",
        nrow(counts) == nrow(genes),
        TRUE,
        if (nrow(counts) == nrow(genes)) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "matrix_columns_equal_cell_metadata_rows",
        cells_match,
        TRUE,
        if (cells_match) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "matrix_negative_values",
        matrix_negative,
        0,
        if (matrix_negative == 0L) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "matrix_fractional_values",
        matrix_fractional,
        0,
        if (matrix_fractional == 0L) {
          "PASS"
        } else {
          "WARN"
        },
        "A referência deverá representar raw counts."
      ),

      make_check(
        "TAG_missing",
        n_tag_missing,
        0,
        if (n_tag_missing == 0L) {
          "PASS"
        } else {
          "WARN"
        }
      ),

      make_check(
        "TAG_duplicated_rows",
        n_tag_duplicated,
        0,
        if (n_tag_duplicated == 0L) {
          "PASS"
        } else {
          "WARN"
        }
      ),

      make_check(
        "unique_projid",
        n_donors_projid,
        EXPECTED_DONORS,
        if (
          n_donors_projid ==
            EXPECTED_DONORS
        ) {
          "PASS"
        } else {
          "WARN"
        },
        "Esperamos 48 doadores; confirmar com o metadata."
      ),

      make_check(
        "eight_expected_cell_types_present",
        length(
          setdiff(
            EXPECTED_CELL_TYPES,
            observed_cell_types
          )
        ) == 0L,
        TRUE,
        if (
          length(
            setdiff(
              EXPECTED_CELL_TYPES,
              observed_cell_types
            )
          ) == 0L
        ) {
          "PASS"
        } else {
          "WARN"
        }
      ),

      make_check(
        "pathology_free_24_candidate_identified",
        pathology_24_found,
        TRUE,
        if (pathology_24_found) {
          "PASS"
        } else {
          "WARN"
        },
        paste0(
          "O 00 não filtra indivíduos. ",
          "Se WARN, revisar 00_pathology_candidate_value_counts.tsv ",
          "e/ou obter metadata auxiliar com status neuropatológico."
        )
      ),

      make_check(
        "donor_overlap_ROSMAP_assessed",
        overlap_summary$overlap_assessed[1],
        TRUE,
        if (
          isTRUE(
            overlap_summary$overlap_assessed[1]
          )
        ) {
          "PASS"
        } else {
          "WARN"
        },
        overlap_summary$note[1]
      )
    ),
    fill = TRUE
  )

  fwrite(
    checks,
    file.path(
      out_dir,
      "00_audit_status.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(
    checks
  )


  # ============================================================================
  # 4.11 CONCLUSÃO
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )
  cat(
    "11. CONCLUSÃO\n"
  )
  cat(
    "------------------------------------------------------------\n"
  )

  cat(
    "- Nenhum gene foi removido.\n"
  )
  cat(
    "- Nenhum núcleo foi removido.\n"
  )
  cat(
    "- Nenhum doador foi classificado automaticamente como pathology-free.\n"
  )
  cat(
    "- Nenhuma normalização foi aplicada.\n"
  )
  cat(
    "- Nenhum identificador gênico foi convertido.\n"
  )
  cat(
    "- Nenhum identificador individual foi exportado nos outputs agregados.\n"
  )

  cat(
    "\nPróxima etapa:\n"
  )
  cat(
    "  01_prepare_metadata_Mathys2019.R\n"
  )
  cat(
    "\nEssa etapa só será escrita após revisarmos qual campo/categoria\n"
  )
  cat(
    "define de maneira inequívoca os 24 doadores pathology-free.\n"
  )


  # ============================================================================
  # 4.12 SESSION INFO
  # ============================================================================

  sink(
    session_file
  )

  print(
    sessionInfo()
  )

  sink()

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
    run_audit()
  },
  error = function(e) {

    execution_error <<- e

    cat(
      "\n============================================================\n"
    )
    cat(
      "ERRO NA AUDITORIA\n"
    )
    cat(
      "============================================================\n"
    )
    cat(
      conditionMessage(e),
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
    conditionMessage(execution_error),
    call. = FALSE
  )
}

cat(
  "\n============================================================\n"
)
cat(
  "00_audit_Mathys2019.R CONCLUÍDO\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Relatório:\n  ",
  log_file,
  "\n\n",
  sep = ""
)

cat(
  "Status dos checks:\n  ",
  file.path(
    out_dir,
    "00_audit_status.tsv"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Candidatos a pathology-free:\n  ",
  file.path(
    out_dir,
    "00_pathology_candidate_value_counts.tsv"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Overlap Mathys × ROSMAP:\n  ",
  file.path(
    out_dir,
    "00_donor_overlap_ROSMAP_summary.tsv"
  ),
  "\n",
  sep = ""
)
