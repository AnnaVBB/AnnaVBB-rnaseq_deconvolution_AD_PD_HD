#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(stringr)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Uso: Rscript 01_parse_metadata.R <soft.gz> <counts.tsv.gz> <out_samples.csv>")
}

soft_file   <- args[1]
counts_file <- args[2]
output_csv  <- args[3]

counts_header <- colnames(fread(cmd = paste("zcat", counts_file), nrows = 1, data.table = FALSE))
valid_samples <- counts_header[-1]

gds <- getGEO(filename = soft_file)
gsm_list <- GSMList(gds)

extract_char_value <- function(chars, key_pattern) {
  matched <- grep(key_pattern, chars, value = TRUE, ignore.case = TRUE)
  if (length(matched) == 0) return(NA_character_)
  val <- str_split_fixed(matched[1], ":", 2)[, 2]
  return(str_trim(val))
}

harmonize_diagnosis <- function(raw_diag) {
  if (is.na(raw_diag)) return(NA_character_)
  diag_clean <- tolower(raw_diag)
  
  if (str_detect(diag_clean, "control|normal|non-demented|neurologically normal")) {
    return("Control")
  } else if (str_detect(diag_clean, "alzheimer|ad")) {
    return("AD")
  } else if (str_detect(diag_clean, "huntington|hd")) {
    return("HD")
  } else if (str_detect(diag_clean, "parkinson|pd")) {
    return("PD")
  } else {
    return("Other")
  }
}

meta_rows <- list()

for (gsm_id in names(gsm_list)) {
  if (!gsm_id %in% valid_samples) next
  
  gsm <- gsm_list[[gsm_id]]
  meta <- Meta(gsm)
  chars <- meta$characteristics_ch1
  
  raw_disease <- extract_char_value(chars, "disease|diagnosis|status|condition")
  age_val     <- extract_char_value(chars, "age")
  sex_val     <- extract_char_value(chars, "sex|gender")
  rin_val     <- extract_char_value(chars, "rin|rna integrity number")
  pmi_val     <- extract_char_value(chars, "pmi|post-mortem interval")
  region_val  <- extract_char_value(chars, "brain region|region|tissue")
  
  meta_rows[[gsm_id]] <- data.frame(
    sample_id     = gsm_id,
    title         = meta$title,
    group         = harmonize_diagnosis(raw_disease),
    raw_diagnosis = raw_disease,
    age           = as.numeric(str_extract(age_val, "\\d+(\\.\\d+)?")),
    sex           = ifelse(str_detect(tolower(sex_val), "^m"), "M", ifelse(str_detect(tolower(sex_val), "^f"), "F", NA_character_)),
    rin           = as.numeric(str_extract(rin_val, "\\d+(\\.\\d+)?")),
    pmi           = as.numeric(str_extract(pmi_val, "\\d+(\\.\\d+)?")),
    region        = region_val,
    stringsAsFactors = FALSE
  )
}

sample_metadata <- do.call(rbind, meta_rows)
write_csv(sample_metadata, output_csv)
cat("[SUCCESS] Metadados extraídos e salvos em: ", output_csv, "\n")
