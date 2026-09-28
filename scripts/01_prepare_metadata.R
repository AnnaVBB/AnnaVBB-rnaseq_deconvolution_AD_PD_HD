#!/usr/bin/env Rscript

# ==============================================================================
# 01_prepare_metadata.R
#
# Objetivo:
#   Preparar e padronizar os metadados das amostras dos três datasets:
#
#   AD - GSE53697
#   HD - GSE64810
#   PD - GSE68719
#
# Este script:
#   1. identifica as amostras presentes nas matrizes de expressão;
#   2. lê os arquivos GEO family.soft.gz;
#   3. associa cada coluna da matriz ao respectivo registro GSM;
#   4. padroniza os grupos experimentais;
#   5. extrai covariáveis clínicas/técnicas disponíveis;
#   6. preserva variáveis clínicas específicas de HD;
#   7. identifica reanálises de controles HD no dataset de PD;
#   8. valida correspondência, grupos e ordem das amostras;
#   9. gera metadados padronizados.
#
# IMPORTANTE:
#   - Nenhuma amostra é removida.
#   - Nenhuma matriz de contagens é modificada.
#   - Nenhuma covariável é incluída em modelo estatístico nesta etapa.
#   - Valores ausentes no GEO permanecem como NA.
#
# Saídas:
#   results/metadata/metadata_AD.tsv
#   results/metadata/metadata_HD.tsv
#   results/metadata/metadata_PD.tsv
#   results/metadata/metadata_all.tsv
# ==============================================================================


# ==============================================================================
# 1. PACOTES
# ==============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(stringr)
})


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

ad_expr <- "data/raw/GSE53697/GSE53697_RNAseq_AD.txt.gz"
hd_expr <- "data/raw/GSE64810/GSE64810_raw_counts_GRCh38.p13_NCBI.tsv.gz"
pd_expr <- "data/raw/GSE68719/GSE68719_raw_counts_GRCh38.p13_NCBI.tsv.gz"

ad_soft <- "data/raw/GSE53697/GSE53697_family.soft.gz"
hd_soft <- "data/raw/GSE64810/GSE64810_family.soft.gz"
pd_soft <- "data/raw/GSE68719/GSE68719_family.soft.gz"

output_dir <- "results/metadata"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 3. VERIFICAÇÃO DOS ARQUIVOS DE ENTRADA
# ==============================================================================

input_files <- c(
  ad_expr,
  hd_expr,
  pd_expr,
  ad_soft,
  hd_soft,
  pd_soft
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Arquivos de entrada não encontrados:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}


# ==============================================================================
# 4. FUNÇÕES AUXILIARES
# ==============================================================================

# ------------------------------------------------------------------------------
# Colapsar um campo de metadados do GEO em uma única string
# ------------------------------------------------------------------------------

collapse_meta <- function(x) {

  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }

  x <- as.character(x)
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = " | ")
}


# ------------------------------------------------------------------------------
# Extrair uma característica do campo characteristics_ch1
#
# Exemplo:
#   "age of death: 73 | diagnosis: Neurologically normal | pmi: 2"
#
# keys permite diferentes nomes possíveis para a mesma variável.
# ------------------------------------------------------------------------------

extract_characteristic <- function(text, keys) {

  if (
    length(text) == 0 ||
    is.na(text) ||
    trimws(text) == ""
  ) {
    return(NA_character_)
  }

  parts <- trimws(unlist(strsplit(text, "\\s*\\|\\s*")))

  for (part in parts) {

    if (!grepl(":", part, fixed = TRUE)) {
      next
    }

    # Campo antes do primeiro ":"
    field <- trimws(sub(":.*$", "", part))

    # Remove eventual texto entre parênteses no nome do campo.
    # Exemplo:
    # "post-mortem interval (pmi)" -> "post-mortem interval"
    field_without_parentheses <- trimws(
      sub("\\s*\\([^)]*\\)\\s*$", "", field)
    )

    # Valor depois do primeiro ":"
    value <- trimws(sub("^[^:]+:\\s*", "", part))

    field_candidates <- tolower(
      c(field, field_without_parentheses)
    )

    if (any(field_candidates %in% tolower(keys))) {

      if (
        value == "" ||
        tolower(value) %in% c("na", "n/a", "nan", "unknown")
      ) {
        return(NA_character_)
      }

      return(value)
    }
  }

  NA_character_
}


# ------------------------------------------------------------------------------
# Extrair primeiro número de uma string
#
# Exemplos:
#   "73 years" -> 73
#   "37.25"    -> 37.25
#   "NA"       -> NA
# ------------------------------------------------------------------------------

extract_number <- function(x) {

  if (length(x) == 0 || is.na(x)) {
    return(NA_real_)
  }

  x <- trimws(as.character(x))

  if (
    x == "" ||
    tolower(x) %in% c("na", "n/a", "nan", "unknown")
  ) {
    return(NA_real_)
  }

  number <- stringr::str_extract(
    x,
    "-?[0-9]+(?:\\.[0-9]+)?"
  )

  suppressWarnings(as.numeric(number))
}


# ------------------------------------------------------------------------------
# Padronizar sexo
# ------------------------------------------------------------------------------

standardize_sex <- function(x) {

  if (length(x) == 0 || is.na(x)) {
    return(NA_character_)
  }

  value <- tolower(trimws(as.character(x)))

  if (value %in% c("male", "m", "man")) {
    return("Male")
  }

  if (value %in% c("female", "f", "woman")) {
    return("Female")
  }

  if (
    value == "" ||
    value %in% c("na", "n/a", "unknown")
  ) {
    return(NA_character_)
  }

  as.character(x)
}


# ------------------------------------------------------------------------------
# Extrair metadados básicos de um objeto GSM
# ------------------------------------------------------------------------------

extract_gsm_metadata <- function(gsm_object, gsm_id) {

  meta <- Meta(gsm_object)

  title <- collapse_meta(meta$title)
  source_name <- collapse_meta(meta$source_name_ch1)
  characteristics <- collapse_meta(meta$characteristics_ch1)

  data.table(
    geo_accession = gsm_id,
    title = title,
    source_name = source_name,
    characteristics = characteristics
  )
}


# ------------------------------------------------------------------------------
# Ler todos os GSMs de um arquivo SOFT
# ------------------------------------------------------------------------------

read_soft_metadata <- function(path) {

  cat("Lendo SOFT:", path, "\n")

  geo <- GEOquery::getGEO(
    filename = path,
    GSEMatrix = FALSE
  )

  gsms <- GEOquery::GSMList(geo)

  cat("  GSMs encontrados:", length(gsms), "\n")

  result <- rbindlist(
    lapply(
      names(gsms),
      function(gsm_id) {
        extract_gsm_metadata(
          gsms[[gsm_id]],
          gsm_id
        )
      }
    ),
    fill = TRUE
  )

  if (anyDuplicated(result$geo_accession)) {
    stop("GSM duplicado encontrado no SOFT: ", path)
  }

  result
}


# ==============================================================================
# 5. IDENTIFICAR AS AMOSTRAS NAS MATRIZES
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("IDENTIFICAÇÃO DAS AMOSTRAS NAS MATRIZES\n")
cat("============================================\n")


# ------------------------------------------------------------------------------
# AD
# ------------------------------------------------------------------------------

ad_header <- names(
  fread(
    ad_expr,
    nrows = 0,
    check.names = FALSE
  )
)

ad_raw_samples <- grep(
  "_raw$",
  ad_header,
  value = TRUE
)

ad_rpkm_samples <- grep(
  "_rpkm$",
  ad_header,
  value = TRUE
)

cat("\nAD / GSE53697\n")
cat("  Colunas raw:", length(ad_raw_samples), "\n")
cat("  Colunas RPKM:", length(ad_rpkm_samples), "\n")

if (length(ad_raw_samples) != 17) {
  stop(
    "AD deveria possuir 17 colunas _raw, mas foram encontradas ",
    length(ad_raw_samples),
    "."
  )
}

if (length(ad_rpkm_samples) != 17) {
  stop(
    "AD deveria possuir 17 colunas _rpkm, mas foram encontradas ",
    length(ad_rpkm_samples),
    "."
  )
}


# ------------------------------------------------------------------------------
# HD
# ------------------------------------------------------------------------------

hd_header <- names(
  fread(
    hd_expr,
    nrows = 0,
    check.names = FALSE
  )
)

hd_samples <- hd_header[-1]

cat("\nHD / GSE64810\n")
cat("  Amostras:", length(hd_samples), "\n")

if (length(hd_samples) != 69) {
  stop(
    "HD deveria possuir 69 amostras, mas foram encontradas ",
    length(hd_samples),
    "."
  )
}

if (!all(grepl("^GSM[0-9]+$", hd_samples))) {
  stop("Existem colunas inesperadas na matriz HD.")
}


# ------------------------------------------------------------------------------
# PD
# ------------------------------------------------------------------------------

pd_header <- names(
  fread(
    pd_expr,
    nrows = 0,
    check.names = FALSE
  )
)

pd_samples <- pd_header[-1]

cat("\nPD / GSE68719\n")
cat("  Amostras:", length(pd_samples), "\n")

if (length(pd_samples) != 73) {
  stop(
    "PD deveria possuir 73 amostras, mas foram encontradas ",
    length(pd_samples),
    "."
  )
}

if (!all(grepl("^GSM[0-9]+$", pd_samples))) {
  stop("Existem colunas inesperadas na matriz PD.")
}


# ==============================================================================
# 6. LER OS ARQUIVOS SOFT
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("LEITURA DOS METADADOS GEO\n")
cat("============================================\n")

ad_geo <- read_soft_metadata(ad_soft)
hd_geo <- read_soft_metadata(hd_soft)
pd_geo <- read_soft_metadata(pd_soft)


# ==============================================================================
# 7. EXTRAIR COVARIÁVEIS COMUNS
# ==============================================================================

add_common_covariates <- function(dt) {

  dt[, raw_diagnosis := vapply(
    characteristics,
    extract_characteristic,
    character(1),
    keys = c(
      "diagnosis",
      "disease status",
      "disease",
      "disease state",
      "neurological status",
      "condition"
    )
  )]

  dt[, age := vapply(
    characteristics,
    function(x) {
      extract_number(
        extract_characteristic(
          x,
          c(
            "age",
            "age at death",
            "age of death"
          )
        )
      )
    },
    numeric(1)
  )]

  dt[, sex := vapply(
    characteristics,
    function(x) {
      standardize_sex(
        extract_characteristic(
          x,
          c(
            "sex",
            "gender"
          )
        )
      )
    },
    character(1)
  )]

  dt[, rin := vapply(
    characteristics,
    function(x) {
      extract_number(
        extract_characteristic(
          x,
          c(
            "rin",
            "rna integrity number"
          )
        )
      )
    },
    numeric(1)
  )]

  dt[, pmi := vapply(
    characteristics,
    function(x) {
      extract_number(
        extract_characteristic(
          x,
          c(
            "pmi",
            "post mortem interval",
            "post-mortem interval"
          )
        )
      )
    },
    numeric(1)
  )]

  dt[, region := vapply(
    characteristics,
    extract_characteristic,
    character(1),
    keys = c(
      "tissue",
      "tissue subtype",
      "brain region",
      "region"
    )
  )]

  dt
}


ad_geo <- add_common_covariates(ad_geo)
hd_geo <- add_common_covariates(hd_geo)
pd_geo <- add_common_covariates(pd_geo)


# ==============================================================================
# 8. PREPARAR METADATA AD
# ==============================================================================

# No arquivo de expressão:
#   C1_raw ... C8_raw  = controles
#   A9_raw ... A17_raw = AD
#
# No SOFT:
#   RNAseq_Ctrl_1
#   RNAseq_AD_9
#   etc.

ad_geo[, study_sample_id := fifelse(
  str_detect(title, regex("^RNAseq_Ctrl_", ignore_case = TRUE)),
  paste0(
    "C",
    str_extract(title, "[0-9]+$")
  ),
  fifelse(
    str_detect(title, regex("^RNAseq_AD_", ignore_case = TRUE)),
    paste0(
      "A",
      str_extract(title, "[0-9]+$")
    ),
    NA_character_
  )
)]

ad_geo[, sample_id := fifelse(
  !is.na(study_sample_id),
  paste0(study_sample_id, "_raw"),
  NA_character_
)]

ad_geo[, group := fifelse(
  str_detect(study_sample_id, "^C"),
  "Control",
  fifelse(
    str_detect(study_sample_id, "^A"),
    "AD",
    NA_character_
  )
)]

ad_geo[, dataset := "GSE53697"]
ad_geo[, original_gsm := NA_character_]


# Associar exatamente às colunas raw da matriz AD

ad_metadata <- merge(
  data.table(
    sample_id = ad_raw_samples,
    matrix_order = seq_along(ad_raw_samples)
  ),
  ad_geo,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

setorder(ad_metadata, matrix_order)


# ==============================================================================
# 9. PREPARAR METADATA HD
# ==============================================================================

# Na matriz HD, as colunas já são GSM.
# O title contém identificadores como:
#   C_0002
#   H_0726

hd_geo[, study_sample_id := str_extract(
  title,
  "^[CH]_[0-9]+"
)]

hd_geo[, group := fifelse(
  str_detect(study_sample_id, "^C_"),
  "Control",
  fifelse(
    str_detect(study_sample_id, "^H_"),
    "HD",
    NA_character_
  )
)]

hd_geo[, sample_id := geo_accession]
hd_geo[, dataset := "GSE64810"]
hd_geo[, original_gsm := NA_character_]


# ------------------------------------------------------------------------------
# Variáveis clínicas específicas de HD
#
# A auditoria direta do SOFT confirmou:
#   age of onset       = 20 registros
#   cag                = 20 registros
#   duration           = 20 registros
#   h-v cortical score = 20 registros
#   h-v striatal score = 20 registros
#   vonsattel grade    = 20 registros
#
# São preservadas no metadata, mas NÃO são utilizadas automaticamente
# no modelo estatístico.
# ------------------------------------------------------------------------------

hd_geo[, age_of_onset := vapply(
  characteristics,
  function(x) {
    extract_number(
      extract_characteristic(
        x,
        "age of onset"
      )
    )
  },
  numeric(1)
)]

hd_geo[, cag := vapply(
  characteristics,
  function(x) {
    extract_number(
      extract_characteristic(
        x,
        "cag"
      )
    )
  },
  numeric(1)
)]

hd_geo[, duration := vapply(
  characteristics,
  function(x) {
    extract_number(
      extract_characteristic(
        x,
        "duration"
      )
    )
  },
  numeric(1)
)]

hd_geo[, hv_cortical_score := vapply(
  characteristics,
  function(x) {
    extract_number(
      extract_characteristic(
        x,
        "h-v cortical score"
      )
    )
  },
  numeric(1)
)]

hd_geo[, hv_striatal_score := vapply(
  characteristics,
  function(x) {
    extract_number(
      extract_characteristic(
        x,
        "h-v striatal score"
      )
    )
  },
  numeric(1)
)]

hd_geo[, vonsattel_grade := vapply(
  characteristics,
  function(x) {
    extract_number(
      extract_characteristic(
        x,
        "vonsattel grade"
      )
    )
  },
  numeric(1)
)]


hd_metadata <- merge(
  data.table(
    sample_id = hd_samples,
    matrix_order = seq_along(hd_samples)
  ),
  hd_geo,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

setorder(hd_metadata, matrix_order)


# ==============================================================================
# 10. PREPARAR METADATA PD
# ==============================================================================

# Exemplos de títulos:
#
# Controle:
#   C_0002 [reanalysis of GSM1580869]
#
# PD:
#   P_0003

pd_geo[, study_sample_id := str_extract(
  title,
  "^[CP]_[0-9]+"
)]

pd_geo[, group := fifelse(
  str_detect(study_sample_id, "^C_"),
  "Control",
  fifelse(
    str_detect(study_sample_id, "^P_"),
    "PD",
    NA_character_
  )
)]

pd_geo[, original_gsm := str_extract(
  title,
  regex(
    "(?<=reanalysis of )GSM[0-9]+",
    ignore_case = TRUE
  )
)]

pd_geo[, sample_id := geo_accession]
pd_geo[, dataset := "GSE68719"]


pd_metadata <- merge(
  data.table(
    sample_id = pd_samples,
    matrix_order = seq_along(pd_samples)
  ),
  pd_geo,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

setorder(pd_metadata, matrix_order)


# ==============================================================================
# 11. VALIDAÇÕES
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("VALIDAÇÃO DOS GRUPOS\n")
cat("============================================\n")


# ------------------------------------------------------------------------------
# Função geral de validação
# ------------------------------------------------------------------------------

validate_metadata <- function(
  metadata,
  expected_samples,
  dataset_name
) {

  if (nrow(metadata) != length(expected_samples)) {
    stop(
      dataset_name,
      ": número de linhas do metadata diferente do número de amostras."
    )
  }

  if (any(is.na(metadata$geo_accession))) {
    stop(
      dataset_name,
      ": existem amostras da matriz sem correspondência GSM."
    )
  }

  if (any(is.na(metadata$group))) {
    stop(
      dataset_name,
      ": existem amostras sem grupo definido."
    )
  }

  if (anyDuplicated(metadata$sample_id)) {
    stop(
      dataset_name,
      ": sample_id duplicado."
    )
  }

  if (!identical(metadata$sample_id, expected_samples)) {
    stop(
      dataset_name,
      ": ordem do metadata não corresponde à ordem da matriz."
    )
  }

  invisible(TRUE)
}


validate_metadata(
  ad_metadata,
  ad_raw_samples,
  "AD"
)

validate_metadata(
  hd_metadata,
  hd_samples,
  "HD"
)

validate_metadata(
  pd_metadata,
  pd_samples,
  "PD"
)


# ------------------------------------------------------------------------------
# Validar números esperados dos grupos
# ------------------------------------------------------------------------------

ad_groups <- table(ad_metadata$group)
hd_groups <- table(hd_metadata$group)
pd_groups <- table(pd_metadata$group)

cat("\nAD:\n")
print(ad_groups)

cat("\nHD:\n")
print(hd_groups)

cat("\nPD:\n")
print(pd_groups)


if (
  unname(ad_groups["Control"]) != 8 ||
  unname(ad_groups["AD"]) != 9
) {
  stop("Distribuição inesperada dos grupos em AD.")
}

if (
  unname(hd_groups["Control"]) != 49 ||
  unname(hd_groups["HD"]) != 20
) {
  stop("Distribuição inesperada dos grupos em HD.")
}

if (
  unname(pd_groups["Control"]) != 44 ||
  unname(pd_groups["PD"]) != 29
) {
  stop("Distribuição inesperada dos grupos em PD.")
}


# ==============================================================================
# 12. AUDITORIA DAS REANÁLISES PD / HD
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("REANÁLISES PD / HD\n")
cat("============================================\n")

pd_reanalysis <- pd_metadata[
  !is.na(original_gsm)
]

n_reanalysis <- nrow(pd_reanalysis)

n_reanalysis_in_hd <- sum(
  pd_reanalysis$original_gsm %in%
    hd_metadata$geo_accession
)

cat(
  "Amostras PD com GSM original explicitamente informado:",
  n_reanalysis,
  "\n"
)

cat(
  "GSMs originais também presentes no metadata HD:",
  n_reanalysis_in_hd,
  "/",
  n_reanalysis,
  "\n"
)

if (n_reanalysis != 44) {
  stop(
    "Esperadas 44 reanálises no dataset PD; encontradas ",
    n_reanalysis,
    "."
  )
}

if (n_reanalysis_in_hd != n_reanalysis) {
  stop(
    "Nem todos os GSMs originais do PD foram encontrados no metadata HD."
  )
}


# ==============================================================================
# 13. AUDITORIA DOS CAMPOS CLÍNICOS ESPECÍFICOS DE HD
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("CAMPOS CLÍNICOS ESPECÍFICOS DE HD\n")
cat("============================================\n")

hd_specific_fields <- c(
  "age_of_onset",
  "cag",
  "duration",
  "hv_cortical_score",
  "hv_striatal_score",
  "vonsattel_grade"
)

for (field in hd_specific_fields) {

  n_available <- sum(
    !is.na(hd_metadata[[field]])
  )

  cat(
    sprintf(
      "%-20s %d / %d\n",
      paste0(field, ":"),
      n_available,
      nrow(hd_metadata)
    )
  )
}


# Como a auditoria direta do SOFT mostrou 20 registros para cada
# variável clínica específica de HD, validamos esse número.

for (field in hd_specific_fields) {

  n_available <- sum(
    !is.na(hd_metadata[[field]])
  )

  if (n_available != 20) {
    stop(
      "Campo HD '",
      field,
      "' deveria possuir 20 valores disponíveis, mas foram encontrados ",
      n_available,
      "."
    )
  }
}


# ==============================================================================
# 14. AUDITORIA DO PMI DE HD
# ==============================================================================

# O SOFT contém o campo PMI para as 69 amostras, porém uma delas
# apresenta explicitamente "pmi: NA".
#
# Portanto, esperamos exatamente 68 valores numéricos disponíveis.

n_hd_pmi <- sum(!is.na(hd_metadata$pmi))

if (n_hd_pmi != 68) {
  stop(
    "Esperados 68 valores numéricos de PMI em HD; encontrados ",
    n_hd_pmi,
    "."
  )
}


# ==============================================================================
# 15. SELECIONAR COLUNAS FINAIS
# ==============================================================================

# Campos comuns aos três datasets

common_columns <- c(
  "dataset",
  "sample_id",
  "geo_accession",
  "study_sample_id",
  "original_gsm",
  "title",
  "group",
  "raw_diagnosis",
  "age",
  "sex",
  "rin",
  "pmi",
  "region"
)


# AD

ad_final <- copy(ad_metadata)

for (field in hd_specific_fields) {
  ad_final[, (field) := NA_real_]
}

ad_final <- ad_final[
  ,
  c(
    common_columns,
    hd_specific_fields
  ),
  with = FALSE
]


# HD

hd_final <- hd_metadata[
  ,
  c(
    common_columns,
    hd_specific_fields
  ),
  with = FALSE
]


# PD

pd_final <- copy(pd_metadata)

for (field in hd_specific_fields) {
  pd_final[, (field) := NA_real_]
}

pd_final <- pd_final[
  ,
  c(
    common_columns,
    hd_specific_fields
  ),
  with = FALSE
]


# ==============================================================================
# 16. METADATA COMBINADO
# ==============================================================================

metadata_all <- rbindlist(
  list(
    ad_final,
    hd_final,
    pd_final
  ),
  use.names = TRUE,
  fill = TRUE
)

if (nrow(metadata_all) != 159) {
  stop(
    "metadata_all deveria possuir 159 amostras, mas possui ",
    nrow(metadata_all),
    "."
  )
}

if (anyDuplicated(paste(
  metadata_all$dataset,
  metadata_all$sample_id,
  sep = "::"
))) {
  stop(
    "Existem combinações dataset/sample_id duplicadas."
  )
}


# ==============================================================================
# 17. RESUMO DOS METADADOS
# ==============================================================================

summarize_metadata <- function(dt, label) {

  cat("\n")
  cat("--------------------------------------------\n")
  cat(label, "\n")
  cat("--------------------------------------------\n")

  cat(
    "Amostras:",
    nrow(dt),
    "\n"
  )

  cat("\nGrupos:\n\n")
  print(table(dt$group))

  cat(
    "\nGSM:",
    sum(!is.na(dt$geo_accession)),
    "/",
    nrow(dt),
    "\n"
  )

  cat(
    "Diagnóstico informado:",
    sum(!is.na(dt$raw_diagnosis)),
    "/",
    nrow(dt),
    "\n"
  )

  cat(
    "Idade:",
    sum(!is.na(dt$age)),
    "/",
    nrow(dt),
    "\n"
  )

  cat(
    "Sexo:",
    sum(!is.na(dt$sex)),
    "/",
    nrow(dt),
    "\n"
  )

  cat(
    "RIN:",
    sum(!is.na(dt$rin)),
    "/",
    nrow(dt),
    "\n"
  )

  cat(
    "PMI:",
    sum(!is.na(dt$pmi)),
    "/",
    nrow(dt),
    "\n"
  )

  cat(
    "Região:",
    sum(!is.na(dt$region)),
    "/",
    nrow(dt),
    "\n"
  )
}


cat("\n")
cat("============================================\n")
cat("RESUMO DOS METADADOS\n")
cat("============================================\n")

summarize_metadata(
  ad_final,
  "AD / GSE53697"
)

summarize_metadata(
  hd_final,
  "HD / GSE64810"
)

summarize_metadata(
  pd_final,
  "PD / GSE68719"
)


# ==============================================================================
# 18. GRAVAR ARQUIVOS
# ==============================================================================

fwrite(
  ad_final,
  file.path(
    output_dir,
    "metadata_AD.tsv"
  ),
  sep = "\t",
  na = "NA"
)

fwrite(
  hd_final,
  file.path(
    output_dir,
    "metadata_HD.tsv"
  ),
  sep = "\t",
  na = "NA"
)

fwrite(
  pd_final,
  file.path(
    output_dir,
    "metadata_PD.tsv"
  ),
  sep = "\t",
  na = "NA"
)

fwrite(
  metadata_all,
  file.path(
    output_dir,
    "metadata_all.tsv"
  ),
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 19. FINALIZAÇÃO
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("01_prepare_metadata.R CONCLUÍDO\n")
cat("============================================\n")

cat(
  "AD:",
  nrow(ad_final),
  "amostras\n"
)

cat(
  "HD:",
  nrow(hd_final),
  "amostras\n"
)

cat(
  "PD:",
  nrow(pd_final),
  "amostras\n"
)

cat(
  "Total:",
  nrow(metadata_all),
  "amostras\n"
)

cat("\nArquivos gerados:\n")

cat(
  " -",
  file.path(
    output_dir,
    "metadata_AD.tsv"
  ),
  "\n"
)

cat(
  " -",
  file.path(
    output_dir,
    "metadata_HD.tsv"
  ),
  "\n"
)

cat(
  " -",
  file.path(
    output_dir,
    "metadata_PD.tsv"
  ),
  "\n"
)

cat(
  " -",
  file.path(
    output_dir,
    "metadata_all.tsv"
  ),
  "\n"
)

cat("\nNenhuma amostra foi removida.\n")
cat("Nenhuma matriz de contagens foi modificada.\n")
cat("Nenhuma covariável foi incluída em modelo estatístico nesta etapa.\n")