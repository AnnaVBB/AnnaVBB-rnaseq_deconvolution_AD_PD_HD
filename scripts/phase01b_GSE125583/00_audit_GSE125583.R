#!/usr/bin/env Rscript

# ==============================================================================
# 00_audit_GSE125583.R
#
# FASE 01b — AUDITORIA DO DATASET CANDIDATO DE ALZHEIMER (GSE125583)
#
# OBJETIVOS
# ------------------------------------------------------------------------------
# 1. Auditar a matriz de raw counts padronizada pelo NCBI.
# 2. Verificar GeneID, integridade dos counts e estrutura das amostras.
# 3. Auditar informações essenciais do SOFT sem ainda construir o metadata final.
# 4. Validar correspondência entre amostras da matriz e amostras do SOFT.
# 5. Comparar, em nível BRUTO, o universo de GeneIDs de GSE125583 com HD e PD.
#
# IMPORTANTE
# ------------------------------------------------------------------------------
# - Este script NÃO filtra genes.
# - Este script NÃO normaliza counts.
# - Este script NÃO executa DESeq2.
# - Este script NÃO constrói o metadata final.
# - Este script NÃO decide ainda se GSE125583 substituirá GSE53697.
#
# A comparação de GeneIDs feita aqui é uma auditoria estrutural da entrada.
# A comparação formal dos universos após preparação/filtro será feita depois.
#
# ENTRADAS
# ------------------------------------------------------------------------------
# data/raw/GSE125583/GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz
# data/raw/GSE125583/GSE125583_family.soft.gz
# data/raw/GSE64810/GSE64810_raw_counts_GRCh38.p13_NCBI.tsv.gz
# data/raw/GSE68719/GSE68719_raw_counts_GRCh38.p13_NCBI.tsv.gz
#
# SAÍDAS
# ------------------------------------------------------------------------------
# results/candidate_GSE125583/audit/00_audit_GSE125583.txt
# results/candidate_GSE125583/audit/00_gene_universe_raw.tsv
# ==============================================================================

required_packages <- c("data.table")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Pacotes necessários não instalados: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# 1. CAMINHOS
# ==============================================================================

ad_counts_file <- file.path(
  "data", "raw", "GSE125583",
  "GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)

ad_soft_file <- file.path(
  "data", "raw", "GSE125583",
  "GSE125583_family.soft.gz"
)

hd_counts_file <- file.path(
  "data", "raw", "GSE64810",
  "GSE64810_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)

pd_counts_file <- file.path(
  "data", "raw", "GSE68719",
  "GSE68719_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)

output_dir <- file.path(
  "results", "candidate_GSE125583", "audit"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

audit_file <- file.path(
  output_dir,
  "00_audit_GSE125583.txt"
)

gene_universe_file <- file.path(
  output_dir,
  "00_gene_universe_raw.tsv"
)

# ==============================================================================
# 2. FUNÇÕES AUXILIARES
# ==============================================================================

validate_file <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo não encontrado: ", path)
  }

  info <- file.info(path)

  if (is.na(info$size) || info$size == 0) {
    stop("Arquivo vazio ou tamanho indisponível: ", path)
  }

  invisible(TRUE)
}

validate_gene_id <- function(gene_id, dataset_name) {
  gene_id <- as.character(gene_id)

  missing_gene_id <- (
    is.na(gene_id) |
      trimws(gene_id) == ""
  )

  if (any(missing_gene_id)) {
    stop(dataset_name, ": existem GeneIDs ausentes.")
  }

  if (anyDuplicated(gene_id)) {
    stop(dataset_name, ": existem GeneIDs duplicados.")
  }

  invisible(TRUE)
}

read_gene_ids <- function(path, dataset_name) {
  dt <- fread(
    path,
    select = "GeneID",
    showProgress = FALSE
  )

  if (!"GeneID" %in% names(dt)) {
    stop(dataset_name, ": coluna GeneID não encontrada.")
  }

  gene_id <- as.character(dt$GeneID)

  validate_gene_id(
    gene_id,
    dataset_name
  )

  unique(gene_id)
}

normalize_missing_character <- function(x) {
  x <- trimws(as.character(x))

  x[
    x %in% c("", "NA", "N/A", "na", "n/a")
  ] <- NA_character_

  x
}

extract_characteristic <- function(block, key) {
  prefix <- paste0(
    "!Sample_characteristics_ch1 = ",
    key,
    ":"
  )

  hits <- block[
    startsWith(block, prefix)
  ]

  if (length(hits) == 0) {
    return(NA_character_)
  }

  trimws(
    substring(
      hits[1],
      nchar(prefix) + 1L
    )
  )
}

extract_single <- function(block, prefix) {
  hits <- block[
    startsWith(block, prefix)
  ]

  if (length(hits) == 0) {
    return(NA_character_)
  }

  trimws(
    substring(
      hits[1],
      nchar(prefix) + 1L
    )
  )
}

parse_soft_metadata <- function(path) {
  con <- gzfile(
    path,
    open = "rt"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  lines <- readLines(
    con,
    warn = FALSE
  )

  sample_start <- grep(
    "^\\^SAMPLE = ",
    lines
  )

  if (length(sample_start) == 0) {
    stop("Nenhuma seção ^SAMPLE encontrada no SOFT.")
  }

  sample_end <- c(
    sample_start[-1] - 1L,
    length(lines)
  )

  out <- vector(
    "list",
    length(sample_start)
  )

  for (i in seq_along(sample_start)) {
    block <- lines[
      sample_start[i]:sample_end[i]
    ]

    sample_id <- sub(
      "^\\^SAMPLE = ",
      "",
      block[1]
    )

    out[[i]] <- data.table(
      sample_id = trimws(sample_id),
      title = extract_single(
        block,
        "!Sample_title = "
      ),
      tissue = extract_characteristic(
        block,
        "tissue"
      ),
      case_id = extract_characteristic(
        block,
        "case id"
      ),
      diagnosis = extract_characteristic(
        block,
        "diagnosis"
      ),
      age = extract_characteristic(
        block,
        "age"
      ),
      sex = extract_characteristic(
        block,
        "Sex"
      ),
      braak_score = extract_characteristic(
        block,
        "braak.score"
      )
    )
  }

  metadata <- rbindlist(
    out,
    fill = TRUE
  )

  for (
    col_name in c(
      "title",
      "tissue",
      "case_id",
      "diagnosis",
      "age",
      "sex",
      "braak_score"
    )
  ) {
    set(
      metadata,
      j = col_name,
      value = normalize_missing_character(
        metadata[[col_name]]
      )
    )
  }

  metadata[
    ,
    age_numeric := suppressWarnings(
      as.numeric(age)
    )
  ]

  metadata
}

format_int <- function(x) {
  format(
    x,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

format_num <- function(x, digits = 4) {
  format(
    round(x, digits),
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

# ==============================================================================
# 3. VALIDAR ARQUIVOS DE ENTRADA
# ==============================================================================

input_files <- c(
  ad_counts_file,
  ad_soft_file,
  hd_counts_file,
  pd_counts_file
)

invisible(
  lapply(
    input_files,
    validate_file
  )
)

# ==============================================================================
# 4. LER MATRIZ DE COUNTS DO GSE125583
# ==============================================================================

counts <- fread(
  ad_counts_file,
  showProgress = TRUE
)

if (ncol(counts) < 2) {
  stop("GSE125583: matriz de counts possui menos de duas colunas.")
}

if (names(counts)[1] != "GeneID") {
  stop(
    "GSE125583: a primeira coluna deveria ser GeneID, mas é: ",
    names(counts)[1]
  )
}

gene_id <- as.character(
  counts$GeneID
)

validate_gene_id(
  gene_id,
  "GSE125583"
)

sample_ids_counts <- setdiff(
  names(counts),
  "GeneID"
)

if (anyDuplicated(sample_ids_counts)) {
  stop("GSE125583: existem IDs de amostras duplicados na matriz.")
}

# ==============================================================================
# 5. AUDITAR COUNTS
# ==============================================================================

sample_dt <- counts[
  ,
  setdiff(names(counts), "GeneID"),
  with = FALSE
]

non_numeric_columns <- names(sample_dt)[
  !vapply(
    sample_dt,
    is.numeric,
    logical(1)
  )
]

if (length(non_numeric_columns) > 0) {
  stop(
    "GSE125583: colunas de amostras não numéricas: ",
    paste(non_numeric_columns, collapse = ", ")
  )
}

count_matrix <- as.matrix(sample_dt)

rm(sample_dt)

n_missing_counts <- sum(
  is.na(count_matrix)
)

n_nonfinite_counts <- sum(
  !is.finite(count_matrix),
  na.rm = TRUE
)

n_negative_counts <- sum(
  count_matrix < 0,
  na.rm = TRUE
)

fractional_mask <- (
  !is.na(count_matrix) &
    is.finite(count_matrix) &
    abs(
      count_matrix -
        round(count_matrix)
    ) > 1e-8
)

n_fractional_counts <- sum(
  fractional_mask
)

if (n_missing_counts > 0) {
  stop(
    "GSE125583: foram encontrados ",
    n_missing_counts,
    " counts ausentes."
  )
}

if (n_nonfinite_counts > 0) {
  stop(
    "GSE125583: foram encontrados ",
    n_nonfinite_counts,
    " counts não finitos."
  )
}

if (n_negative_counts > 0) {
  stop(
    "GSE125583: foram encontrados ",
    n_negative_counts,
    " counts negativos."
  )
}

if (n_fractional_counts > 0) {
  stop(
    "GSE125583: foram encontrados ",
    n_fractional_counts,
    " counts fracionários. ",
    "Esperavam-se raw counts inteiros da pipeline NCBI."
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

all_zero_genes <- sum(
  row_totals == 0
)

all_zero_samples <- sum(
  library_sizes == 0
)

n_zero_counts <- sum(
  count_matrix == 0
)

total_count_cells <- length(
  count_matrix
)

zero_fraction <- (
  n_zero_counts /
    total_count_cells
)

# ==============================================================================
# 6. LER E AUDITAR SOFT
# ==============================================================================

metadata_audit <- parse_soft_metadata(
  ad_soft_file
)

if (anyDuplicated(metadata_audit$sample_id)) {
  stop("GSE125583 SOFT: existem sample IDs duplicados.")
}

sample_ids_soft <- metadata_audit$sample_id

samples_only_counts <- setdiff(
  sample_ids_counts,
  sample_ids_soft
)

samples_only_soft <- setdiff(
  sample_ids_soft,
  sample_ids_counts
)

if (
  length(samples_only_counts) > 0 ||
    length(samples_only_soft) > 0
) {
  stop(
    paste0(
      "GSE125583: counts e SOFT não contêm exatamente as mesmas amostras.\n",
      "Somente counts: ",
      paste(samples_only_counts, collapse = ", "),
      "\nSomente SOFT: ",
      paste(samples_only_soft, collapse = ", ")
    )
  )
}

same_sample_order <- identical(
  sample_ids_counts,
  sample_ids_soft
)

unexpected_diagnosis <- setdiff(
  unique(
    na.omit(
      metadata_audit$diagnosis
    )
  ),
  c(
    "Alzheimer's disease",
    "control"
  )
)

if (length(unexpected_diagnosis) > 0) {
  stop(
    "GSE125583: diagnósticos inesperados no SOFT: ",
    paste(unexpected_diagnosis, collapse = ", ")
  )
}

# ==============================================================================
# 7. COMPARAR UNIVERSO BRUTO DE GeneID COM HD E PD
# ==============================================================================

ids_ad <- gene_id

ids_hd <- read_gene_ids(
  hd_counts_file,
  "GSE64810"
)

ids_pd <- read_gene_ids(
  pd_counts_file,
  "GSE68719"
)

intersection_ad_hd <- intersect(
  ids_ad,
  ids_hd
)

intersection_ad_pd <- intersect(
  ids_ad,
  ids_pd
)

intersection_hd_pd <- intersect(
  ids_hd,
  ids_pd
)

intersection_all <- Reduce(
  intersect,
  list(
    ids_ad,
    ids_hd,
    ids_pd
  )
)

union_all <- Reduce(
  union,
  list(
    ids_ad,
    ids_hd,
    ids_pd
  )
)

gene_universe_audit <- data.table(
  metric = c(
    "n_genes_GSE125583",
    "n_genes_GSE64810",
    "n_genes_GSE68719",
    "GSE125583_intersect_GSE64810",
    "GSE125583_intersect_GSE68719",
    "GSE64810_intersect_GSE68719",
    "triple_intersection",
    "union_all",
    "GSE125583_equals_GSE64810",
    "GSE125583_equals_GSE68719",
    "GSE64810_equals_GSE68719"
  ),
  value = c(
    length(ids_ad),
    length(ids_hd),
    length(ids_pd),
    length(intersection_ad_hd),
    length(intersection_ad_pd),
    length(intersection_hd_pd),
    length(intersection_all),
    length(union_all),
    as.character(
      setequal(ids_ad, ids_hd)
    ),
    as.character(
      setequal(ids_ad, ids_pd)
    ),
    as.character(
      setequal(ids_hd, ids_pd)
    )
  )
)

fwrite(
  gene_universe_audit,
  gene_universe_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ==============================================================================
# 8. GERAR RELATÓRIO DE AUDITORIA
# ==============================================================================

diagnosis_counts <- metadata_audit[
  ,
  .N,
  by = diagnosis
][
  order(diagnosis)
]

sex_counts <- metadata_audit[
  ,
  .N,
  by = sex
][
  order(sex)
]

braak_levels <- c(
  "I", "II", "III", "IV", "V", "VI"
)

braak_counts <- metadata_audit[
  ,
  .N,
  by = braak_score
]

braak_counts[
  ,
  braak_order := match(
    braak_score,
    braak_levels
  )
]

setorder(
  braak_counts,
  braak_order
)

braak_counts[
  ,
  braak_order := NULL
]

tissue_counts <- metadata_audit[
  ,
  .N,
  by = tissue
][
  order(tissue)
]

diagnosis_sex_counts <- metadata_audit[
  ,
  .N,
  by = .(
    diagnosis,
    sex
  )
][
  order(
    diagnosis,
    sex
  )
]

diagnosis_braak_counts <- metadata_audit[
  ,
  .N,
  by = .(
    diagnosis,
    braak_score
  )
]

diagnosis_braak_counts[
  ,
  braak_order := match(
    braak_score,
    braak_levels
  )
]

setorder(
  diagnosis_braak_counts,
  diagnosis,
  braak_order
)

diagnosis_braak_counts[
  ,
  braak_order := NULL
]

case_id_missing <- sum(
  is.na(
    metadata_audit$case_id
  )
)

case_id_nonmissing <- metadata_audit[
  !is.na(case_id),
  case_id
]

case_id_duplicated <- sum(
  duplicated(
    case_id_nonmissing
  )
)

age_missing <- sum(
  is.na(
    metadata_audit$age_numeric
  )
)

age_valid <- metadata_audit[
  !is.na(age_numeric),
  age_numeric
]

sink(audit_file)

cat("============================================================\n")
cat("00 — AUDITORIA GSE125583\n")
cat("============================================================\n\n")

cat("DATASET CANDIDATO DE AD\n")
cat("GEO: GSE125583\n")
cat("Entrada principal: NCBI raw counts GRCh38.p13\n\n")

cat("------------------------------------------------------------\n")
cat("1. ESTRUTURA DA MATRIZ GSE125583\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes: ",
  format_int(nrow(counts)),
  "\n",
  sep = ""
)

cat(
  "Amostras: ",
  format_int(length(sample_ids_counts)),
  "\n",
  sep = ""
)

cat(
  "Colunas totais: ",
  format_int(ncol(counts)),
  "\n",
  sep = ""
)

cat(
  "Primeira coluna: ",
  names(counts)[1],
  "\n",
  sep = ""
)

cat("GeneIDs ausentes: 0\n")
cat("GeneIDs duplicados: 0\n\n")

cat("------------------------------------------------------------\n")
cat("2. INTEGRIDADE DOS RAW COUNTS\n")
cat("------------------------------------------------------------\n")

cat(
  "Counts ausentes: ",
  format_int(n_missing_counts),
  "\n",
  sep = ""
)

cat(
  "Counts não finitos: ",
  format_int(n_nonfinite_counts),
  "\n",
  sep = ""
)

cat(
  "Counts negativos: ",
  format_int(n_negative_counts),
  "\n",
  sep = ""
)

cat(
  "Counts fracionários: ",
  format_int(n_fractional_counts),
  "\n",
  sep = ""
)

cat(
  "Genes totalmente zerados: ",
  format_int(all_zero_genes),
  "\n",
  sep = ""
)

cat(
  "Amostras totalmente zeradas: ",
  format_int(all_zero_samples),
  "\n",
  sep = ""
)

cat(
  "Proporção global de zeros: ",
  format_num(
    100 * zero_fraction,
    digits = 2
  ),
  "%\n",
  sep = ""
)

cat(
  "Menor library size: ",
  format_int(min(library_sizes)),
  "\n",
  sep = ""
)

cat(
  "Mediana library size: ",
  format_int(median(library_sizes)),
  "\n",
  sep = ""
)

cat(
  "Maior library size: ",
  format_int(max(library_sizes)),
  "\n",
  sep = ""
)

cat(
  "Genes detectados/amostra — mínimo: ",
  format_int(min(genes_detected_per_sample)),
  "\n",
  sep = ""
)

cat(
  "Genes detectados/amostra — mediana: ",
  format_int(median(genes_detected_per_sample)),
  "\n",
  sep = ""
)

cat(
  "Genes detectados/amostra — máximo: ",
  format_int(max(genes_detected_per_sample)),
  "\n\n",
  sep = ""
)

cat("------------------------------------------------------------\n")
cat("3. SOFT / AMOSTRAS / METADADOS ESSENCIAIS\n")
cat("------------------------------------------------------------\n")

cat(
  "Amostras no SOFT: ",
  format_int(nrow(metadata_audit)),
  "\n",
  sep = ""
)

cat("Mesmos sample IDs em counts e SOFT: TRUE\n")

cat(
  "Mesma ordem em counts e SOFT: ",
  same_sample_order,
  "\n",
  sep = ""
)

cat("\nDiagnóstico:\n")
print(diagnosis_counts)

cat("\nSexo:\n")
print(sex_counts)

cat("\nBraak score:\n")
print(braak_counts)

cat("\nTecido:\n")
print(tissue_counts)

cat("\nDiagnóstico × sexo:\n")
print(diagnosis_sex_counts)

cat("\nDiagnóstico × Braak score:\n")
print(diagnosis_braak_counts)

cat(
  "\nCase IDs ausentes: ",
  format_int(case_id_missing),
  "\n",
  sep = ""
)

cat(
  "Case IDs duplicados entre valores não ausentes: ",
  format_int(case_id_duplicated),
  "\n",
  sep = ""
)

cat(
  "\nIdades ausentes/não numéricas: ",
  format_int(age_missing),
  "\n",
  sep = ""
)

if (length(age_valid) > 0) {
  cat(
    "Idade mínima: ",
    min(age_valid),
    "\n",
    sep = ""
  )

  cat(
    "Idade mediana: ",
    median(age_valid),
    "\n",
    sep = ""
  )

  cat(
    "Idade máxima: ",
    max(age_valid),
    "\n",
    sep = ""
  )
}

cat("\nIdade por diagnóstico:\n")
print(
  metadata_audit[
    ,
    .(
      n = .N,
      n_age_missing = sum(
        is.na(age_numeric)
      ),
      age_min = if (
        all(is.na(age_numeric))
      ) {
        NA_real_
      } else {
        min(
          age_numeric,
          na.rm = TRUE
        )
      },
      age_median = if (
        all(is.na(age_numeric))
      ) {
        NA_real_
      } else {
        median(
          age_numeric,
          na.rm = TRUE
        )
      },
      age_max = if (
        all(is.na(age_numeric))
      ) {
        NA_real_
      } else {
        max(
          age_numeric,
          na.rm = TRUE
        )
      }
    ),
    by = diagnosis
  ]
)

cat("\n")

cat("------------------------------------------------------------\n")
cat("4. UNIVERSO BRUTO DE GeneID — AD NOVO vs HD vs PD\n")
cat("------------------------------------------------------------\n")

cat(
  "GSE125583: ",
  format_int(length(ids_ad)),
  " genes\n",
  sep = ""
)

cat(
  "GSE64810: ",
  format_int(length(ids_hd)),
  " genes\n",
  sep = ""
)

cat(
  "GSE68719: ",
  format_int(length(ids_pd)),
  " genes\n\n",
  sep = ""
)

cat(
  "GSE125583 ∩ GSE64810: ",
  format_int(length(intersection_ad_hd)),
  "\n",
  sep = ""
)

cat(
  "GSE125583 ∩ GSE68719: ",
  format_int(length(intersection_ad_pd)),
  "\n",
  sep = ""
)

cat(
  "GSE64810 ∩ GSE68719: ",
  format_int(length(intersection_hd_pd)),
  "\n",
  sep = ""
)

cat(
  "GSE125583 ∩ GSE64810 ∩ GSE68719: ",
  format_int(length(intersection_all)),
  "\n",
  sep = ""
)

cat(
  "União dos três universos: ",
  format_int(length(union_all)),
  "\n\n",
  sep = ""
)

cat(
  "GSE125583 == GSE64810: ",
  setequal(ids_ad, ids_hd),
  "\n",
  sep = ""
)

cat(
  "GSE125583 == GSE68719: ",
  setequal(ids_ad, ids_pd),
  "\n",
  sep = ""
)

cat(
  "GSE64810 == GSE68719: ",
  setequal(ids_hd, ids_pd),
  "\n",
  sep = ""
)

cat("\n")

cat("------------------------------------------------------------\n")
cat("5. CONCLUSÃO ESTRUTURAL DA AUDITORIA\n")
cat("------------------------------------------------------------\n")

cat(
  "- A matriz GSE125583 contém raw counts inteiros, não negativos e sem valores ausentes.\n"
)

cat(
  "- Counts e SOFT contêm o mesmo conjunto de amostras.\n"
)

cat(
  "- O universo bruto de GeneIDs de GSE125583 pode ser comparado diretamente aos universos NCBI de HD e PD.\n"
)

if (
  setequal(ids_ad, ids_hd) &&
    setequal(ids_ad, ids_pd)
) {
  cat(
    "- Os três datasets possuem exatamente o mesmo conjunto bruto de GeneIDs.\n"
  )
} else {
  cat(
    "- Os três datasets NÃO possuem universos brutos de GeneIDs idênticos; revisar 00_gene_universe_raw.tsv.\n"
  )
}

cat(
  "- Nenhuma decisão sobre exclusão de amostras, covariáveis, filtragem ou substituição de GSE53697 foi tomada nesta etapa.\n"
)

cat("\n")
cat("============================================================\n")
cat("00_audit_GSE125583.R CONCLUÍDO\n")
cat("============================================================\n")

sink()

# ==============================================================================
# 9. RESUMO NO TERMINAL
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("00 — AUDITORIA GSE125583 CONCLUÍDA\n")
cat("============================================================\n")

cat(
  "Genes: ",
  format_int(nrow(counts)),
  "\n",
  sep = ""
)

cat(
  "Amostras: ",
  format_int(length(sample_ids_counts)),
  "\n",
  sep = ""
)

cat(
  "Genes totalmente zerados: ",
  format_int(all_zero_genes),
  "\n",
  sep = ""
)

cat(
  "Counts fracionários: ",
  format_int(n_fractional_counts),
  "\n",
  sep = ""
)

cat(
  "Universo bruto idêntico a HD: ",
  setequal(ids_ad, ids_hd),
  "\n",
  sep = ""
)

cat(
  "Universo bruto idêntico a PD: ",
  setequal(ids_ad, ids_pd),
  "\n",
  sep = ""
)

cat(
  "Relatório: ",
  audit_file,
  "\n",
  sep = ""
)

cat(
  "Auditoria de GeneIDs: ",
  gene_universe_file,
  "\n",
  sep = ""
)

cat("============================================================\n")
