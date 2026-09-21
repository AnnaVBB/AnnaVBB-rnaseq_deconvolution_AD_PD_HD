#!/usr/bin/env Rscript

# ============================================================
# 01_parse_metadata.R
# Construção e padronização dos metadados das amostras para os três datasets do núcleo:
#   GSE53697 -> Alzheimer (AD)
#   GSE64810 -> Huntington (HD)
#   GSE68719 -> Parkinson (PD)
#
# Os sample_id produzidos por este script correspondem aos nomes das amostras presentes nos arquivos processados disponibilizados pelos autores.
# Uso:
# Rscript scripts/01_parse_metadata.R \
#   <dataset> \
#   <expression_file> \
#   <soft_file> \
#   <output_csv>
# ============================================================

# ------------------------------------------------------------
# 1. Pacotes
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
})


# ------------------------------------------------------------
# 2. Argumentos
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    paste0(
      "Uso:\n",
      "Rscript scripts/01_parse_metadata.R ",
      "<dataset> <expression_file> <soft_file> <output_csv>"
    )
  )
}

dataset <- args[1]
expression_file <- args[2]
soft_file <- args[3]
output_csv <- args[4]


# ------------------------------------------------------------
# 3. Verificações iniciais
# ------------------------------------------------------------

valid_datasets <- c(
  "GSE53697",
  "GSE64810",
  "GSE68719"
)

if (!dataset %in% valid_datasets) {
  stop(
    "Dataset não reconhecido: ",
    dataset
  )
}

if (!file.exists(expression_file)) {
  stop(
    "Arquivo de expressão não encontrado: ",
    expression_file
  )
}

if (!file.exists(soft_file)) {
  stop(
    "Arquivo SOFT não encontrado: ",
    soft_file
  )
}


# ------------------------------------------------------------
# 4. Criar diretório de saída
# ------------------------------------------------------------

output_dir <- dirname(output_csv)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------
# 5. Ler matriz de expressão
# ------------------------------------------------------------
cat("Lendo matriz de expressão...\n")
expr <- fread(
  expression_file,
  check.names = FALSE
)

# ------------------------------------------------------------
# 6. Identificar amostras presentes na matriz
# ------------------------------------------------------------
if (dataset == "GSE53697") {
  control_samples <- grep(
    "^C[0-9]+_raw$",
    colnames(expr),
    value = TRUE
  )
  disease_samples <- grep(
    "^A[0-9]+_raw$",
    colnames(expr),
    value = TRUE
  )
  expected_control <- 8
  expected_disease <- 9

} else if (dataset == "GSE64810") {
  control_samples <- grep(
    "^C_",
    colnames(expr),
    value = TRUE
  )
  disease_samples <- grep(
    "^H_",
    colnames(expr),
    value = TRUE
  )
  expected_control <- 49
  expected_disease <- 20

} else if (dataset == "GSE68719") {
  control_samples <- grep(
    "^C_",
    colnames(expr),
    value = TRUE
  )
  disease_samples <- grep(
    "^P_",
    colnames(expr),
    value = TRUE
  )
  expected_control <- 44
  expected_disease <- 29
}
sample_ids <- c(
  control_samples,
  disease_samples
)
groups <- case_when(

  dataset == "GSE53697" & str_detect(sample_ids, "^C[0-9]+_raw$") ~ "Control",
  dataset == "GSE53697" & str_detect(sample_ids, "^A[0-9]+_raw$") ~ "AD",

  dataset == "GSE64810" & str_detect(sample_ids, "^C_") ~ "Control",
  dataset == "GSE64810" & str_detect(sample_ids, "^H_") ~ "HD",

  dataset == "GSE68719" & str_detect(sample_ids, "^C_") ~ "Control",
  dataset == "GSE68719" & str_detect(sample_ids, "^P_") ~ "PD",

  TRUE ~ NA_character_
)

if (any(is.na(groups))) {
  stop(
    paste0(
      "Não foi possível definir o grupo para: ",
      paste(sample_ids[is.na(groups)], collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 7. Verificar número de amostras
# ------------------------------------------------------------
if (length(control_samples) != expected_control) {
  stop(
    "Número inesperado de controles: ",
    length(control_samples),
    ". Esperado: ",
    expected_control
  )
}

if (length(disease_samples) != expected_disease) {
  stop(
    "Número inesperado de casos: ",
    length(disease_samples),
    ". Esperado: ",
    expected_disease
  )
}

# ------------------------------------------------------------
# 8. Definir grupo a partir da matriz
# ------------------------------------------------------------
if (dataset == "GSE53697") {
  group <- ifelse(
    grepl("^C", sample_ids),
    "Control",
    "AD"
  )

} else if (dataset == "GSE64810") {
  group <- ifelse(
    grepl("^C_", sample_ids),
    "Control",
    "HD"
  )

} else if (dataset == "GSE68719") {
  group <- ifelse(
    grepl("^C_", sample_ids),
    "Control",
    "PD"
  )
}

# ------------------------------------------------------------
# 9. Ler arquivo SOFT
# ------------------------------------------------------------
cat("Lendo arquivo SOFT...\n")

gse <- GEOquery::getGEO(
  filename = soft_file,
  getGPL = FALSE
)

gsm_list <- GEOquery::GSMList(gse)

cat(
  "Amostras encontradas no SOFT:",
  length(gsm_list),
  "\n"
)


# ------------------------------------------------------------
# 10. Funções auxiliares
# ------------------------------------------------------------

collapse_meta <- function(x) {

  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }

  paste(
    as.character(x),
    collapse = " | "
  )
}


# Extrai uma característica procurando múltiplos nomes
# possíveis para o mesmo campo.
extract_characteristic <- function(text, keys) {

  if (
    is.na(text) ||
    text == ""
  ) {
    return(NA_character_)
  }

  parts <- unlist(
    strsplit(
      text,
      "\\s*\\|\\s*"
    )
  )

  parts <- trimws(parts)

  for (key in keys) {

    # Escapa caracteres especiais da chave para uso em regex
    key_regex <- stringr::str_replace_all(
      key,
      "([.\\\\+*?\\[\\](){}^$|])",
      "\\\\\\1"
    )

    pattern <- paste0(
      "^\\s*",
      key_regex,
      "(?:\\s*\\([^)]*\\))?",
      "\\s*:"
    )

    hit <- grep(
      pattern,
      parts,
      ignore.case = TRUE,
      value = TRUE
    )

    if (length(hit) > 0) {
      value <- sub(
        "^[^:]+:\\s*",
        "",
        hit[1]
      )

      return(
        trimws(value)
      )
    }
  }

  NA_character_
}

extract_number <- function(x) {
  if (
    is.na(x) ||
    x == ""
  ) {
    return(NA_real_)
  }

  value <- str_extract(
    x,
    "-?[0-9]+(?:\\.[0-9]+)?"
  )

  suppressWarnings(
    as.numeric(value)
  )
}

# ------------------------------------------------------------
# 11. Extrair metadata bruto de cada GSM
# ------------------------------------------------------------

extract_gsm_metadata <- function(gsm) {

  meta <- Meta(gsm)

  geo_accession <- collapse_meta(
    meta$geo_accession
  )

  title <- collapse_meta(
    meta$title
  )

  source <- collapse_meta(
    meta$source_name_ch1
  )

  characteristics <- collapse_meta(
    meta$characteristics_ch1
  )


  # ----------------------------------------------------------
  # Diagnóstico
  # ----------------------------------------------------------

  diagnosis <- extract_characteristic(
    characteristics,
    c(
      "disease status",
      "disease",
      "diagnosis",
      "disease state",
      "neurological status",
      "condition"
    )
  )


  # ----------------------------------------------------------
  # Idade
  # ----------------------------------------------------------
  age_raw <- extract_characteristic(
  characteristics,
  c(
    "age",
    "age at death",
    "age of death",
    "age years"
  )
)


  # ----------------------------------------------------------
  # Sexo
  # ----------------------------------------------------------

  sex_raw <- extract_characteristic(
    characteristics,
    c(
      "sex",
      "gender"
    )
  )


  # ----------------------------------------------------------
  # RIN
  #
  # Incluímos variações observadas em diferentes GEOs.
  # ----------------------------------------------------------

  rin_raw <- extract_characteristic(
    characteristics,
    c(
      "rin",
      "rna integrity number",
      "rna integrity"
    )
  )


  # ----------------------------------------------------------
  # PMI
  # ----------------------------------------------------------
  pmi_raw <- extract_characteristic(
    characteristics,
    c(
      "pmi",
      "post mortem interval",
      "post-mortem interval",
      "postmortem interval"
    )
  )

  # ----------------------------------------------------------
  # Região
  # ----------------------------------------------------------
    region_raw <- extract_characteristic(
    characteristics,
    c(
      "tissue subtype",
      "brain region",
      "region",
      "tissue"
    )
  )

  age_num <- extract_number(age_raw)
  rin_num <- extract_number(rin_raw)
  pmi_num <- extract_number(pmi_raw)


  result <- base::data.frame(
    geo_accession = geo_accession,
    title = title,
    source_name = source,
    characteristics = characteristics,
    raw_diagnosis = diagnosis,
    age = age_num,
    sex_raw = sex_raw,
    rin = rin_num,
    pmi = pmi_num,
    region = region_raw,
    stringsAsFactors = FALSE
  )

  return(result)
}


soft_metadata <- bind_rows(
  lapply(
    gsm_list,
    extract_gsm_metadata
  )
)

# ------------------------------------------------------------
# 12. Padronizar sexo
# ------------------------------------------------------------
soft_metadata <- soft_metadata %>%
  mutate(

    sex = case_when(

      str_to_lower(
        trimws(sex_raw)
      ) %in% c(
        "m",
        "male"
      ) ~ "M",

      str_to_lower(
        trimws(sex_raw)
      ) %in% c(
        "f",
        "female"
      ) ~ "F",

      TRUE ~ NA_character_
    )
  )


# ------------------------------------------------------------
# 13. Mapear IDs do arquivo dos autores para GEO
# ------------------------------------------------------------
# ------------------------------------------------------------
# AD — GSE53697
#
# GEO: RNAseq_Ctrl_1 x RNAseq_AD_9
# Matriz: C1_raw x A9_raw
# ------------------------------------------------------------

if (dataset == "GSE53697") {
  soft_metadata <- soft_metadata %>%
    mutate(
      sample_id = case_when(
        str_detect(
          title,
          regex(
            "^RNAseq_Ctrl_[0-9]+$",
            ignore_case = TRUE
          )
        ) ~ paste0(
          "C",
          str_extract(
            title,
            "[0-9]+$"
          ),
          "_raw"
        ),

        str_detect(
          title,
          regex(
            "^RNAseq_AD_[0-9]+$",
            ignore_case = TRUE
          )
        ) ~ paste0(
          "A",
          str_extract(
            title,
            "[0-9]+$"
          ),
          "_raw"
        ),

        TRUE ~ NA_character_
      )
    )
}


# ------------------------------------------------------------
# HD — GSE64810
# ------------------------------------------------------------

if (dataset == "GSE64810") {

  soft_metadata <- soft_metadata %>%
    mutate(

      combined_text = paste(
        title,
        source_name,
        characteristics,
        sep = " | "
      ),

      sample_id = str_extract(
        combined_text,
        regex(
          "(C|H)_[0-9]+",
          ignore_case = TRUE
        )
      ),

      sample_id = toupper(
        sample_id
      )
    )
}


# ------------------------------------------------------------
# PD — GSE68719
# ------------------------------------------------------------
if (dataset == "GSE68719") {
  soft_metadata <- soft_metadata %>%
    mutate(

      combined_text = paste(
        title,
        source_name,
        characteristics,
        sep = " | "
      ),

      sample_id = str_extract(
        combined_text,
        regex(
          "(C|P)_[0-9]+",
          ignore_case = TRUE
        )
      ),
      sample_id = toupper(
        sample_id
      )
    )
}

# ------------------------------------------------------------
# 14. Manter apenas amostras presentes na matriz
# ------------------------------------------------------------
soft_matched <- soft_metadata %>%
  filter(
    !is.na(sample_id),
    sample_id %in% sample_ids
  )


# ------------------------------------------------------------
# 15. Verificar duplicações
# ------------------------------------------------------------
duplicated_mapping <- soft_matched %>%
  count(sample_id) %>%
  filter(n > 1)

if (nrow(duplicated_mapping) > 0) {
  stop(
    paste0(
      "Mapeamento duplicado para: ",
      paste(
        duplicated_mapping$sample_id,
        collapse = ", "
      )
    )
  )
}

# ------------------------------------------------------------
# 16. Metadata-base
# ------------------------------------------------------------

metadata <- data.frame(
  sample_id = sample_ids,
  group = groups,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 17. Adicionar informações GEO
# ------------------------------------------------------------
metadata <- metadata %>%
  left_join(
    soft_matched %>%
      select(
        sample_id,
        geo_accession,
        title,
        raw_diagnosis,
        age,
        sex,
        rin,
        pmi,
        region
      ),
    by = "sample_id"
  )

# ------------------------------------------------------------
# 18. Ordenar colunas
# ------------------------------------------------------------
metadata <- metadata %>%
  select(
    sample_id,
    geo_accession,
    title,
    group,
    raw_diagnosis,
    age,
    sex,
    rin,
    pmi,
    region
  )

# ------------------------------------------------------------
# 19. Garantir ordem idêntica à matriz
# ------------------------------------------------------------
metadata <- metadata %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = sample_ids
    )
  ) %>%

  arrange(sample_id) %>%
  mutate(
    sample_id = as.character(
      sample_id
    )
  )

# ------------------------------------------------------------
# 20. Validações finais
# ------------------------------------------------------------
if (nrow(metadata) != length(sample_ids)) {
  stop(
    "Número de linhas do metadata diferente ",
    "do número de amostras da matriz."
  )
}

if (anyDuplicated(metadata$sample_id)) {
  stop(
    "Há sample_id duplicado."
  )
}

if (!identical(
  metadata$sample_id,
  sample_ids
)) {
  stop(
    "A ordem do metadata não corresponde ",
    "à ordem da matriz."
  )
}

if (any(is.na(metadata$group))) {
  stop(
    "Há amostras sem grupo."
  )
}


# ------------------------------------------------------------
# 21. Verificar correspondência GEO
# ------------------------------------------------------------
n_geo <- sum(
  !is.na(metadata$geo_accession)
)

if (n_geo != nrow(metadata)) {

  warning(
    nrow(metadata) - n_geo,
    " amostras não foram ligadas a GSM."
  )
}


# ------------------------------------------------------------
# 22. Resumo
# ------------------------------------------------------------
cat("\n")
cat("============================================\n")
cat("METADATA\n")
cat("============================================\n")

cat("Dataset:", dataset, "\n")
cat("Amostras:", nrow(metadata), "\n")

cat("\nGrupos:\n")

print(
  table(metadata$group)
)

cat(
  "\nAmostras ligadas a GSM:",
  n_geo,
  "/",
  nrow(metadata),
  "\n"
)


cat("\nCovariáveis disponíveis:\n")

cat(
  "Diagnosis:",
  sum(!is.na(metadata$raw_diagnosis)),
  "/",
  nrow(metadata),
  "\n"
)

cat(
  "Age:",
  sum(!is.na(metadata$age)),
  "/",
  nrow(metadata),
  "\n"
)

cat(
  "Sex:",
  sum(!is.na(metadata$sex)),
  "/",
  nrow(metadata),
  "\n"
)

cat(
  "RIN:",
  sum(!is.na(metadata$rin)),
  "/",
  nrow(metadata),
  "\n"
)

cat(
  "PMI:",
  sum(!is.na(metadata$pmi)),
  "/",
  nrow(metadata),
  "\n"
)

cat(
  "Region:",
  sum(!is.na(metadata$region)),
  "/",
  nrow(metadata),
  "\n"
)


# ------------------------------------------------------------
# 23. Salvar
# ------------------------------------------------------------
write_csv(
  metadata,
  output_csv,
  na = "NA"
)


cat("\nMetadata salvo em:\n")
cat(output_csv, "\n")

cat("============================================\n")