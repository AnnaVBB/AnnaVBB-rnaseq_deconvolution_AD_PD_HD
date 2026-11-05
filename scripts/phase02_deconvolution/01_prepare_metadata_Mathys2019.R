#!/usr/bin/env Rscript

# ==============================================================================
# 01_prepare_metadata_Mathys2019.R
# Fase 02 — Deconvolução celular
# Etapa 01 — Seleção dos 24 doadores "no-pathology"
# ==============================================================================
#
# O script usa a classificação ORIGINAL do estudo Mathys et al. 2019.
#
# Inputs:
#   - filtered_column_metadata.txt
#   - snRNAseqPFC_BA10_id_mapping.csv
#   - snRNAseqPFC_BA10_biospecimen_metadata.csv
#   - 41586_2019_1195_MOESM3_ESM.xlsx  (Supplementary Table 1; sheet S1)
#   - 41586_2019_1195_MOESM5_ESM.xlsx  (Supplementary Table 3; sheet 1)
#
# Fluxo:
#   Supplement 1 + Supplement 3
#       -> Subject
#       -> projid
#       -> individualID
#       -> pathology.group == "no-pathology"
#       -> exatamente 24 doadores
#       -> núcleos desses 24 doadores
#       -> QC donor x cell type
#       -> overlap com o bulk ROSMAP
#
# IMPORTANTE:
#   * NÃO reconstrói pathology-free usando Braak/CERAD/cogdx.
#   * NÃO modifica a matriz de counts.
#   * NÃO normaliza expressão.
#   * NÃO harmoniza genes.
#   * NÃO cria pseudobulk.
#
# PRIVACIDADE:
#   projid e individualID são usados apenas em memória.
#   Outputs analíticos usam donor_id anônimo (RefDonor_01 ... RefDonor_24).
# ==============================================================================


# ==============================================================================
# 0. PACOTES
# ==============================================================================

required_packages <- c("data.table", "openxlsx")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Pacotes ausentes: ",
    paste(missing_packages, collapse = ", "),
    "\nInstale-os antes de executar o script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})


# ==============================================================================
# 1. PARÂMETROS
# ==============================================================================

EXPECTED_TOTAL_DONORS <- 48L
EXPECTED_PATHOLOGY_FREE_DONORS <- 24L
PATHOLOGY_FREE_LABEL <- "no-pathology"

EXPECTED_CELL_TYPES <- c(
  "Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per"
)


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

raw_dir <- file.path("data", "reference", "mathys2019", "raw")

# Diretório novo, separado dos objetos "prepared" do pipeline antigo.
prepared_dir <- file.path(
  "data", "reference", "mathys2019", "phase02_prepared"
)

result_dir <- file.path(
  "results", "deconvolution", "mathys2019", "metadata"
)

dir.create(prepared_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

cell_file <- file.path(raw_dir, "filtered_column_metadata.txt")
idmap_file <- file.path(raw_dir, "snRNAseqPFC_BA10_id_mapping.csv")
biospecimen_file <- file.path(
  raw_dir, "snRNAseqPFC_BA10_biospecimen_metadata.csv"
)
supp1_file <- file.path(raw_dir, "41586_2019_1195_MOESM3_ESM.xlsx")
supp3_file <- file.path(raw_dir, "41586_2019_1195_MOESM5_ESM.xlsx")

bulk_meta_file <- file.path(
  "results", "candidate_ROSMAP", "metadata", "metadata_ROSMAP.tsv"
)

log_file <- file.path(result_dir, "01_prepare_metadata_Mathys2019.txt")
session_file <- file.path(
  result_dir, "01_prepare_metadata_Mathys2019_sessionInfo.txt"
)


# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | toupper(x) == "NA"] <- NA_character_
  x
}

normalize_label <- function(x) {
  x <- tolower(clean_chr(x))
  x <- gsub("_", "-", x, fixed = TRUE)
  x <- gsub(" +", "-", x)
  x <- gsub("-+", "-", x)
  x
}

normalize_colname <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

find_col <- function(dt, candidates, context, required = TRUE) {
  idx <- match(
    normalize_colname(candidates),
    normalize_colname(names(dt))
  )
  idx <- idx[!is.na(idx)]

  if (length(idx) == 0L) {
    if (required) {
      stop(
        context, ": coluna esperada não encontrada. Procuramos: ",
        paste(candidates, collapse = ", "),
        ". Disponíveis: ",
        paste(names(dt), collapse = ", ")
      )
    }
    return(NA_character_)
  }

  names(dt)[idx[1]]
}

assert_one_to_one <- function(dt, key, value, label) {
  x <- dt[
    !is.na(get(key)) & !is.na(get(value)),
    .(
      n_values = uniqueN(get(value))
    ),
    by = key
  ]

  if (any(x$n_values > 1L)) {
    stop(
      label,
      ": relação não é 1:1; há ",
      sum(x$n_values > 1L),
      " chave(s) ambígua(s)."
    )
  }

  invisible(TRUE)
}

detect_bulk_individual_col <- function(dt) {
  candidates <- c(
    "individual_id", "individualID", "individualid",
    "participant_id", "subject_id"
  )
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0L) NA_character_ else hit[1]
}

make_check <- function(check, observed, expected, status, note = NA_character_) {
  data.table(
    check = check,
    observed = as.character(observed),
    expected = as.character(expected),
    status = status,
    note = note
  )
}


# ==============================================================================
# 4. FUNÇÃO PRINCIPAL
# ==============================================================================

run_prepare <- function() {

  cat("============================================================\n")
  cat("FASE 02 — 01_prepare_metadata_Mathys2019.R\n")
  cat("SELEÇÃO DOS 24 DOADORES NO-PATHOLOGY\n")
  cat("============================================================\n\n")


  # ============================================================================
  # 4.1 INPUTS
  # ==============================================================================

  required_files <- c(
    cell_metadata = cell_file,
    id_mapping = idmap_file,
    biospecimen = biospecimen_file,
    supplementary_table_1 = supp1_file,
    supplementary_table_3 = supp3_file
  )

  missing <- required_files[!file.exists(required_files)]

  if (length(missing) > 0L) {
    cat("Arquivos ausentes:\n")
    print(missing)
    stop("Inputs obrigatórios ausentes.")
  }

  input_audit <- data.table(
    role = names(required_files),
    file_name = basename(required_files),
    size_MB = round(file.info(required_files)$size / 1024^2, 3)
  )

  fwrite(
    input_audit,
    file.path(result_dir, "01_input_file_audit.tsv"),
    sep = "\t"
  )

  cat("Todos os inputs obrigatórios foram encontrados.\n")
  print(input_audit)


  # ============================================================================
  # 4.2 METADATA DOS NÚCLEOS
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("1. METADATA DOS NÚCLEOS\n")
  cat("------------------------------------------------------------\n")

  cell_meta <- fread(cell_file, na.strings = c("NA", ""))

  required_cell_cols <- c("TAG", "projid", "broad.cell.type")
  missing_cell_cols <- setdiff(required_cell_cols, names(cell_meta))

  if (length(missing_cell_cols) > 0L) {
    stop(
      "Metadata celular sem colunas: ",
      paste(missing_cell_cols, collapse = ", ")
    )
  }

  cell_meta[, matrix_col_index := seq_len(.N)]
  cell_meta[, projid_internal := clean_chr(projid)]
  cell_meta[, cell_type := clean_chr(broad.cell.type)]

  all_types <- sort(unique(na.omit(cell_meta$cell_type)))

  if (!setequal(all_types, EXPECTED_CELL_TYPES)) {
    stop(
      "Tipos celulares diferentes dos 8 esperados. Observados: ",
      paste(all_types, collapse = ", ")
    )
  }

  n_cell_donors <- uniqueN(na.omit(cell_meta$projid_internal))

  cat("Núcleos:", nrow(cell_meta), "\n")
  cat("Doadores:", n_cell_donors, "\n")
  cat("Tipos:", paste(all_types, collapse = ", "), "\n")

  if (n_cell_donors != EXPECTED_TOTAL_DONORS) {
    stop("Metadata celular não contém os 48 projid esperados.")
  }


  # ============================================================================
  # 4.3 Subject -> projid
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("2. MAPEAMENTO Subject -> projid\n")
  cat("------------------------------------------------------------\n")

  idmap <- fread(idmap_file, na.strings = c("NA", ""))

  subject_col <- find_col(idmap, "Subject", "id_mapping")
  projid_col <- find_col(idmap, "projid", "id_mapping")

  key <- unique(
    idmap[
      ,
      .(
        Subject = clean_chr(get(subject_col)),
        projid = clean_chr(get(projid_col))
      )
    ][
      !is.na(Subject) & !is.na(projid)
    ]
  )

  assert_one_to_one(key, "Subject", "projid", "Subject -> projid")
  assert_one_to_one(key, "projid", "Subject", "projid -> Subject")

  cat("Pares Subject/projid:", nrow(key), "\n")

  if (nrow(key) != EXPECTED_TOTAL_DONORS) {
    stop(
      "Esperávamos 48 pares Subject/projid; encontrados: ",
      nrow(key)
    )
  }


  # ============================================================================
  # 4.4 projid -> individualID
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("3. MAPEAMENTO projid -> individualID\n")
  cat("------------------------------------------------------------\n")

  bio <- fread(biospecimen_file, na.strings = c("NA", ""))

  bio_projid_col <- find_col(bio, "projid", "biospecimen")
  bio_individual_col <- find_col(bio, "individualID", "biospecimen")

  donor_map <- unique(
    bio[
      ,
      .(
        projid = clean_chr(get(bio_projid_col)),
        individualID = clean_chr(get(bio_individual_col))
      )
    ][
      !is.na(projid) & !is.na(individualID)
    ]
  )

  assert_one_to_one(
    donor_map, "projid", "individualID", "projid -> individualID"
  )
  assert_one_to_one(
    donor_map, "individualID", "projid", "individualID -> projid"
  )

  cat("Pares projid/individualID:", nrow(donor_map), "\n")

  if (nrow(donor_map) != EXPECTED_TOTAL_DONORS) {
    stop(
      "Esperávamos 48 pares projid/individualID; encontrados: ",
      nrow(donor_map)
    )
  }


  # ============================================================================
  # 4.5 SUPPLEMENTARY TABLE 1
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("4. SUPPLEMENTARY TABLE 1\n")
  cat("------------------------------------------------------------\n")

  sheets1 <- openxlsx::getSheetNames(supp1_file)
  cat("Sheets:", paste(sheets1, collapse = ", "), "\n")

  if (!"S1" %in% sheets1) {
    stop("Sheet 'S1' não encontrada no Supplementary Table 1.")
  }

  supp1 <- as.data.table(
    openxlsx::read.xlsx(
      supp1_file,
      sheet = "S1",
      colNames = TRUE,
      check.names = FALSE,
      skipEmptyRows = TRUE
    )
  )

  s1_subject_col <- find_col(
    supp1, "Subject", "Supplementary Table 1"
  )

  s1_pathdiag_col <- find_col(
    supp1,
    c(
      "pathologic.diagnosis.of.AD",
      "pathologic diagnosis of AD",
      "pathologic_diagnosis_of_AD"
    ),
    "Supplementary Table 1",
    required = FALSE
  )

  supp1[, Subject_internal := clean_chr(get(s1_subject_col))]

  if (anyNA(supp1$Subject_internal)) {
    stop("Supplementary Table 1 possui Subject ausente.")
  }

  if (anyDuplicated(supp1$Subject_internal)) {
    stop("Supplementary Table 1 possui Subject duplicado.")
  }

  cat("Indivíduos na S1:", nrow(supp1), "\n")

  if (nrow(supp1) != EXPECTED_TOTAL_DONORS) {
    stop("Supplementary Table 1 não contém 48 indivíduos.")
  }

  if (!is.na(s1_pathdiag_col)) {
    diag_counts <- data.table(
      value = clean_chr(supp1[[s1_pathdiag_col]])
    )[
      ,
      .N,
      by = value
    ][order(-N)]

    fwrite(
      diag_counts,
      file.path(result_dir, "01_pathologic_diagnosis_AD_counts.tsv"),
      sep = "\t",
      na = "NA"
    )

    cat("pathologic.diagnosis.of.AD:\n")
    print(diag_counts)
  }


  # ============================================================================
  # 4.6 SUPPLEMENTARY TABLE 3
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("5. SUPPLEMENTARY TABLE 3\n")
  cat("------------------------------------------------------------\n")

  sheets3 <- openxlsx::getSheetNames(supp3_file)
  cat("Sheets:", paste(sheets3, collapse = ", "), "\n")

  supp3 <- as.data.table(
    openxlsx::read.xlsx(
      supp3_file,
      sheet = 1,
      colNames = TRUE,
      check.names = FALSE,
      skipEmptyRows = TRUE
    )
  )

  s3_subject_col <- find_col(
    supp3, "Subject", "Supplementary Table 3"
  )

  s3_pathgroup_col <- find_col(
    supp3,
    c("pathology.group", "pathology group", "pathology_group"),
    "Supplementary Table 3"
  )

  supp3[, Subject_internal := clean_chr(get(s3_subject_col))]
  supp3[, pathology_group_internal := normalize_label(
    get(s3_pathgroup_col)
  )]

  if (anyNA(supp3$Subject_internal)) {
    stop("Supplementary Table 3 possui Subject ausente.")
  }

  if (anyDuplicated(supp3$Subject_internal)) {
    stop("Supplementary Table 3 possui Subject duplicado.")
  }

  pathgroup_counts_raw <- supp3[
    ,
    .(n_donors = .N),
    by = pathology_group_internal
  ][order(pathology_group_internal)]

  fwrite(
    pathgroup_counts_raw,
    file.path(
      result_dir,
      "01_pathology_group_counts_from_supplement.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  cat("pathology.group:\n")
  print(pathgroup_counts_raw)


  # ============================================================================
  # 4.7 COMBINAR SUPLEMENTOS
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("6. INTEGRAÇÃO DOS SUPLEMENTOS\n")
  cat("------------------------------------------------------------\n")

  idx3 <- match(
    supp1$Subject_internal,
    supp3$Subject_internal
  )

  if (anyNA(idx3)) {
    stop("Subjects da S1 sem correspondência na Supplementary Table 3.")
  }

  # Precisamos apenas do pathology.group da Table 3 para a seleção.
  donor_sup <- data.table(
    Subject = supp1$Subject_internal,
    pathology_group = supp3$pathology_group_internal[idx3]
  )

  if (!is.na(s1_pathdiag_col)) {
    donor_sup[
      ,
      pathologic_diagnosis_AD := clean_chr(
        supp1[[s1_pathdiag_col]]
      )
    ]
  } else {
    donor_sup[, pathologic_diagnosis_AD := NA_character_]
  }

  if (nrow(donor_sup) != EXPECTED_TOTAL_DONORS) {
    stop("Integração dos suplementos não resultou em 48 indivíduos.")
  }


  # ============================================================================
  # 4.8 Subject -> projid -> individualID
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("7. MAPEAMENTO DOS 48 INDIVÍDUOS\n")
  cat("------------------------------------------------------------\n")

  idx_key <- match(donor_sup$Subject, key$Subject)

  if (anyNA(idx_key)) {
    stop("Subjects dos suplementos sem projid no id_mapping.")
  }

  donor_sup[, projid_internal := key$projid[idx_key]]

  idx_bio <- match(donor_sup$projid_internal, donor_map$projid)

  if (anyNA(idx_bio)) {
    stop("projid dos suplementos sem individualID no biospecimen.")
  }

  donor_sup[, individualID_internal := donor_map$individualID[idx_bio]]

  if (uniqueN(donor_sup$projid_internal) != EXPECTED_TOTAL_DONORS) {
    stop("Mapeamento final não possui 48 projid únicos.")
  }

  if (uniqueN(donor_sup$individualID_internal) != EXPECTED_TOTAL_DONORS) {
    stop("Mapeamento final não possui 48 individualID únicos.")
  }

  cat("Subjects:", uniqueN(donor_sup$Subject), "\n")
  cat("projid:", uniqueN(donor_sup$projid_internal), "\n")
  cat("individualID:", uniqueN(donor_sup$individualID_internal), "\n")


  # ============================================================================
  # 4.9 SELEÇÃO NO-PATHOLOGY
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("8. SELEÇÃO DOS 24 NO-PATHOLOGY\n")
  cat("------------------------------------------------------------\n")

  pathgroup_counts <- donor_sup[
    ,
    .(
      n_donors = uniqueN(individualID_internal)
    ),
    by = pathology_group
  ][order(pathology_group)]

  fwrite(
    pathgroup_counts,
    file.path(result_dir, "01_pathology_group_counts_final.tsv"),
    sep = "\t",
    na = "NA"
  )

  print(pathgroup_counts)

  pf <- donor_sup[
    pathology_group == PATHOLOGY_FREE_LABEL
  ]

  n_pf <- uniqueN(pf$individualID_internal)

  cat("Doadores no-pathology:", n_pf, "\n")

  if (n_pf != EXPECTED_PATHOLOGY_FREE_DONORS) {
    stop(
      "Esperávamos exatamente 24 no-pathology; encontrados: ",
      n_pf,
      ". O pipeline foi interrompido."
    )
  }

  # IDs internos determinísticos. projid é usado só para ordenar em memória.
  setorder(pf, projid_internal)
  pf[, donor_id := sprintf("RefDonor_%02d", seq_len(.N))]


  # ============================================================================
  # 4.10 CROSS-CHECK COM pathologic.diagnosis.of.AD
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("9. CROSS-CHECK COM PATHOLOGIC DIAGNOSIS OF AD\n")
  cat("------------------------------------------------------------\n")

  crosscheck <- pf[
    ,
    .(
      n_donors = .N
    ),
    by = .(
      pathology_group,
      pathologic_diagnosis_AD
    )
  ]

  fwrite(
    crosscheck,
    file.path(
      result_dir,
      "01_pathology_free_diagnosis_crosscheck.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(crosscheck)

  if (
    any(
      toupper(clean_chr(pf$pathologic_diagnosis_AD)) == "YES",
      na.rm = TRUE
    )
  ) {
    warning(
      "Existe no-pathology com pathologic.diagnosis.of.AD == YES. ",
      "Revisar antes de prosseguir."
    )
  }


  # ============================================================================
  # 4.11 SELECIONAR OS NÚCLEOS DOS 24 DOADORES
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("10. SELEÇÃO DOS NÚCLEOS\n")
  cat("------------------------------------------------------------\n")

  idx_pf <- match(
    cell_meta$projid_internal,
    pf$projid_internal
  )

  cell_meta[, donor_id := pf$donor_id[idx_pf]]
  selected <- cell_meta[!is.na(donor_id)]

  if (uniqueN(selected$donor_id) != EXPECTED_PATHOLOGY_FREE_DONORS) {
    stop("Os núcleos selecionados não representam os 24 doadores.")
  }

  cat("Núcleos totais:", nrow(cell_meta), "\n")
  cat("Núcleos no-pathology:", nrow(selected), "\n")
  cat("Doadores selecionados:", uniqueN(selected$donor_id), "\n")


  # ============================================================================
  # 4.12 COBERTURA donor x cell type
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("11. COBERTURA DONOR × CELL TYPE\n")
  cat("------------------------------------------------------------\n")

  selected_types <- sort(unique(na.omit(selected$cell_type)))
  missing_types <- setdiff(EXPECTED_CELL_TYPES, selected_types)

  if (length(missing_types) > 0L) {
    stop(
      "Tipos celulares ausentes após seleção: ",
      paste(missing_types, collapse = ", ")
    )
  }

  observed <- selected[
    ,
    .(n_nuclei = .N),
    by = .(donor_id, cell_type)
  ]

  donor_cell_qc <- merge(
    CJ(
      donor_id = sort(pf$donor_id),
      cell_type = EXPECTED_CELL_TYPES,
      unique = TRUE
    ),
    observed,
    by = c("donor_id", "cell_type"),
    all.x = TRUE,
    sort = TRUE
  )

  donor_cell_qc[
    is.na(n_nuclei),
    n_nuclei := 0L
  ]

  # as.numeric() evita conflito integer/double do data.table entre grupos.
  coverage <- donor_cell_qc[
    ,
    .(
      donors_total = as.integer(uniqueN(donor_id)),
      donors_present = as.integer(sum(n_nuclei > 0)),
      donors_absent = as.integer(sum(n_nuclei == 0)),
      total_nuclei = as.numeric(sum(n_nuclei)),
      min_nuclei_including_zero = as.numeric(min(n_nuclei)),
      median_nuclei_including_zero = as.numeric(median(n_nuclei)),
      max_nuclei = as.numeric(max(n_nuclei)),
      min_nuclei_when_present = as.numeric(
        if (any(n_nuclei > 0)) min(n_nuclei[n_nuclei > 0]) else 0
      ),
      median_nuclei_when_present = as.numeric(
        if (any(n_nuclei > 0)) median(n_nuclei[n_nuclei > 0]) else 0
      )
    ),
    by = cell_type
  ][order(-total_nuclei)]

  fwrite(
    coverage,
    file.path(
      result_dir,
      "01_pathology_free_celltype_coverage.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  fwrite(
    donor_cell_qc,
    file.path(
      result_dir,
      "01_pathology_free_donor_celltype_qc.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(coverage)

  rare_qc <- coverage[
    cell_type %in% c("End", "Per")
  ]

  cat("\nCobertura específica de End/Per:\n")
  print(rare_qc)


  # ============================================================================
  # 4.13 OVERLAP COM BULK ROSMAP
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("12. OVERLAP: 24 NO-PATHOLOGY × BULK ROSMAP\n")
  cat("------------------------------------------------------------\n")

  pf[, overlap_with_bulk_ROSMAP := FALSE]

  overlap_assessed <- FALSE
  overlap_n <- NA_integer_
  bulk_n <- NA_integer_
  overlap_note <- NA_character_

  if (!file.exists(bulk_meta_file)) {

    overlap_note <- "Metadata bulk ROSMAP não encontrado."
    cat(overlap_note, "\n")

  } else {

    bulk_meta <- fread(
      bulk_meta_file,
      na.strings = c("NA", "")
    )

    bulk_id_col <- detect_bulk_individual_col(bulk_meta)

    if (is.na(bulk_id_col)) {

      overlap_note <- "Coluna individual_id não encontrada no bulk ROSMAP."
      cat(overlap_note, "\n")

    } else {

      bulk_ids <- unique(
        na.omit(
          clean_chr(
            bulk_meta[[bulk_id_col]]
          )
        )
      )

      overlap_ids <- intersect(
        pf$individualID_internal,
        bulk_ids
      )

      pf[
        ,
        overlap_with_bulk_ROSMAP :=
          individualID_internal %in% overlap_ids
      ]

      bulk_n <- length(bulk_ids)
      overlap_n <- length(overlap_ids)
      overlap_assessed <- TRUE

      overlap_note <- paste0(
        "IDs originais não foram exportados; somente contagem e flag ",
        "por donor_id anônimo."
      )

      cat("Doadores no-pathology:", n_pf, "\n")
      cat("Doadores bulk ROSMAP:", bulk_n, "\n")
      cat("Overlap:", overlap_n, "\n")
    }
  }

  overlap_summary <- data.table(
    overlap_assessed = overlap_assessed,
    pathology_free_donors = n_pf,
    bulk_ROSMAP_donors = bulk_n,
    overlap_n = overlap_n,
    note = overlap_note
  )

  fwrite(
    overlap_summary,
    file.path(
      result_dir,
      "01_pathology_free_overlap_ROSMAP_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  fwrite(
    pf[
      ,
      .(
        donor_id,
        overlap_with_bulk_ROSMAP
      )
    ],
    file.path(
      result_dir,
      "01_pathology_free_overlap_by_anonymous_donor.tsv"
    ),
    sep = "\t",
    na = "NA"
  )


  # ============================================================================
  # 4.14 OUTPUTS PARA AS PRÓXIMAS ETAPAS
  # ==============================================================================

  cat("\n------------------------------------------------------------\n")
  cat("13. SALVAR METADATA PREPARADO\n")
  cat("------------------------------------------------------------\n")

  # Metadata de doadores anônimo.
  donor_out <- pf[
    ,
    .(
      donor_id,
      pathology_group,
      pathologic_diagnosis_AD,
      overlap_with_bulk_ROSMAP
    )
  ]

  nuclei_per_donor <- selected[
    ,
    .(n_nuclei = .N),
    by = donor_id
  ]

  donor_out <- merge(
    donor_out,
    nuclei_per_donor,
    by = "donor_id",
    all.x = TRUE,
    sort = FALSE
  )

  setorder(donor_out, donor_id)

  fwrite(
    donor_out,
    file.path(
      prepared_dir,
      "mathys2019_pathology_free_donor_metadata.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  # Adicionar flag de overlap ao metadata dos núcleos.
  selected[
    ,
    overlap_with_bulk_ROSMAP :=
      pf$overlap_with_bulk_ROSMAP[
        match(donor_id, pf$donor_id)
      ]
  ]

  cell_output_cols <- c(
    "matrix_col_index",
    "TAG",
    "donor_id",
    "cell_type",
    "Subcluster",
    "pre.cluster",
    "overlap_with_bulk_ROSMAP"
  )

  cell_output_cols <- cell_output_cols[
    cell_output_cols %in% names(selected)
  ]

  cell_out <- selected[, ..cell_output_cols]

  fwrite(
    cell_out,
    file.path(
      prepared_dir,
      "mathys2019_pathology_free_cell_metadata.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  # Índices das colunas da matriz original a manter.
  selected_columns <- selected$matrix_col_index

  saveRDS(
    selected_columns,
    file.path(
      prepared_dir,
      "mathys2019_pathology_free_matrix_columns.rds"
    ),
    compress = FALSE
  )

  cat("Metadata de doadores:", nrow(donor_out), "linhas\n")
  cat("Metadata de núcleos:", nrow(cell_out), "linhas\n")
  cat("Índices da matriz:", length(selected_columns), "\n")


  # ============================================================================
  # 4.15 RESUMO E CHECKS
  # ==============================================================================

  selection_summary <- data.table(
    metric = c(
      "total_Mathys_donors",
      "pathology_free_donors",
      "total_nuclei",
      "pathology_free_nuclei",
      "cell_types_full_reference",
      "cell_types_pathology_free",
      "bulk_ROSMAP_overlap_assessed",
      "bulk_ROSMAP_overlap_n"
    ),
    value = c(
      EXPECTED_TOTAL_DONORS,
      n_pf,
      nrow(cell_meta),
      nrow(selected),
      length(all_types),
      length(selected_types),
      overlap_assessed,
      overlap_n
    )
  )

  fwrite(
    selection_summary,
    file.path(
      result_dir,
      "01_reference_selection_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  checks <- rbindlist(
    list(
      make_check(
        "mapped_total_donors",
        uniqueN(donor_sup$individualID_internal),
        EXPECTED_TOTAL_DONORS,
        if (uniqueN(donor_sup$individualID_internal) == 48L) "PASS" else "FAIL"
      ),
      make_check(
        "pathology_free_donors",
        n_pf,
        EXPECTED_PATHOLOGY_FREE_DONORS,
        if (n_pf == 24L) "PASS" else "FAIL"
      ),
      make_check(
        "pathology_free_donors_in_cell_metadata",
        uniqueN(selected$donor_id),
        EXPECTED_PATHOLOGY_FREE_DONORS,
        if (uniqueN(selected$donor_id) == 24L) "PASS" else "FAIL"
      ),
      make_check(
        "eight_cell_types_present_after_selection",
        setequal(selected_types, EXPECTED_CELL_TYPES),
        TRUE,
        if (setequal(selected_types, EXPECTED_CELL_TYPES)) "PASS" else "FAIL"
      ),
      make_check(
        "matrix_column_indices_unique",
        anyDuplicated(selected_columns) == 0L,
        TRUE,
        if (anyDuplicated(selected_columns) == 0L) "PASS" else "FAIL"
      ),
      make_check(
        "bulk_ROSMAP_overlap_assessed",
        overlap_assessed,
        TRUE,
        if (overlap_assessed) "PASS" else "WARN",
        overlap_note
      )
    ),
    fill = TRUE
  )

  fwrite(
    checks,
    file.path(
      result_dir,
      "01_prepare_metadata_status.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  cat("\nResumo:\n")
  print(selection_summary)

  cat("\nChecks finais:\n")
  print(checks)


  # ============================================================================
  # 4.16 SESSION INFO E CONCLUSÃO
  # ==============================================================================

  sink(session_file)
  print(sessionInfo())
  sink()

  cat("\n------------------------------------------------------------\n")
  cat("14. CONCLUSÃO\n")
  cat("------------------------------------------------------------\n")

  cat("Doadores totais Mathys:", EXPECTED_TOTAL_DONORS, "\n")
  cat("Doadores no-pathology:", n_pf, "\n")
  cat("Núcleos no-pathology:", nrow(selected), "\n")
  cat("Tipos mantidos:", paste(selected_types, collapse = ", "), "\n")
  cat(
    "Overlap com bulk ROSMAP:",
    if (overlap_assessed) overlap_n else "não avaliado",
    "\n"
  )

  cat("\nNenhuma matriz de counts foi modificada.\n")
  cat("Nenhum gene foi harmonizado nesta etapa.\n")
  cat("Nenhuma normalização foi aplicada.\n")

  cat("\nPróxima etapa:\n")
  cat("  02_harmonize_genes_Mathys2019.R\n")

  invisible(TRUE)
}


# ==============================================================================
# 5. EXECUÇÃO COM LOG
# ==============================================================================

log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)

execution_error <- NULL

tryCatch(
  {
    run_prepare()
  },
  error = function(e) {
    execution_error <<- e
    cat("\n============================================================\n")
    cat("ERRO NA PREPARAÇÃO DO METADATA\n")
    cat("============================================================\n")
    cat(conditionMessage(e), "\n")
  }
)

sink()
close(log_con)

if (!is.null(execution_error)) {
  stop(conditionMessage(execution_error), call. = FALSE)
}

cat("\n============================================================\n")
cat("01_prepare_metadata_Mathys2019.R CONCLUÍDO\n")
cat("============================================================\n\n")

cat("Log:\n  ", log_file, "\n\n", sep = "")
cat(
  "Resumo:\n  ",
  file.path(result_dir, "01_reference_selection_summary.tsv"),
  "\n\n",
  sep = ""
)
cat(
  "Cobertura celular:\n  ",
  file.path(result_dir, "01_pathology_free_celltype_coverage.tsv"),
  "\n\n",
  sep = ""
)
cat(
  "Overlap com ROSMAP:\n  ",
  file.path(result_dir, "01_pathology_free_overlap_ROSMAP_summary.tsv"),
  "\n\n",
  sep = ""
)
cat(
  "Metadata preparado:\n  ",
  file.path(
    prepared_dir,
    "mathys2019_pathology_free_cell_metadata.tsv"
  ),
  "\n",
  sep = ""
)
