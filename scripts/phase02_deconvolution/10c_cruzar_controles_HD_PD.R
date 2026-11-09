#!/usr/bin/env Rscript
# ==============================================================================
# 10c | TRANSFERENCIA AUDITADA DE METADADOS DE CONTROLES COMPARTILHADOS HD/PD
# Projeto TCC: GSE64810 (HD) e GSE68719 (PD), BA9.
#
# FINALIDADE
# - Cruzar EXCLUSIVAMENTE pelo sample_id EXATO (por exemplo C_0002),
#   verificando que ambos os metadados da etapa 05 classificam a amostra
#   como Control.
# - Aproveitar campos efetivamente observados no GEO de HD para os 44
#   controles compartilhados com PD; nao inventar sexo, idade, RIN, PMI.
# - Produzir extensoes CANDIDATAS e relatorios de proveniencia/NA.
# - NAO modificar data/processed nem autorizar modelos automaticamente.
#
# PRE-REQUISITOS
# - Ter executado a versao v2 do 10b, que inclui "age_of_death".
# - SOFT original de GSE64810 local em data/raw/GSE64810/.
#
# USO (na raiz do repositorio):
# Rscript scripts/phase02_deconvolution/10c_cruzar_controles_HD_PD.R
# ==============================================================================
options(stringsAsFactors = FALSE)
if (!requireNamespace("data.table", quietly = TRUE)) stop("Pacote data.table ausente")
suppressPackageStartupMessages(library(data.table))
ROOT <- normalizePath(Sys.getenv("DECONV_ROOT", "."), mustWork = TRUE)
P <- function(...) file.path(ROOT, ...)
DIR10 <- P("results", "deconvolution", "Integrated_Comparison", "10_metadata", "10b_GEO")
OUT <- P("results", "deconvolution", "Integrated_Comparison", "10_metadata", "10c_HD_PD")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
WR <- function(x, name) fwrite(x, file.path(OUT, name), sep = "\t", quote = FALSE, na = "NA")

load_checked <- function(path, required, key = "sample_id") {
  if (!file.exists(path)) stop("Arquivo ausente: ", path)
  x <- fread(path, na.strings = c("NA", "", "NULL", "NaN"))
  if (!all(required %in% names(x))) stop(basename(path), ": campos ausentes: ",
      paste(setdiff(required, names(x)), collapse = ", "))
  x[, (key) := trimws(as.character(get(key)))]
  if (anyNA(x[[key]]) || any(!nzchar(x[[key]])) || anyDuplicated(x, by = key))
    stop("Chaves invalidas/duplicadas no arquivo: ", basename(path))
  x
}

read_meta <- function(ds) {
  f <- P("data", "processed", "deconvolution", ds, "metadata_deconvolution.tsv")
  x <- load_checked(f, c("sample_id", "group"))
  x[, group := trimws(as.character(group))]
  expected <- if (ds == "GSE64810") c("Control", "HD") else c("Control", "PD")
  if (!setequal(unique(x$group), expected)) stop(ds, ": grupos inesperados")
  expected_n <- if (ds == "GSE64810") 69L else 73L
  if (nrow(x) != expected_n) stop(ds, ": numero de amostras inesperado: ", nrow(x))
  x[, .(sample_id, group)]
}

add_candidate <- function(meta, ds) {
  f <- file.path(DIR10, paste0("10b_", ds, "_candidate_extension_REVIEW_FIRST.tsv"))
  c <- load_checked(f, "sample_id")
  extra <- setdiff(c$sample_id, meta$sample_id)
  if (length(extra)) stop("Candidata ", ds, ": IDs nao presentes no metadata 05: ",
      paste(head(extra, 5), collapse = ", "))
  idx <- match(meta$sample_id, c$sample_id)
  out <- copy(meta)
  for (nm in setdiff(names(c), "sample_id")) out[, (nm) := c[[nm]][idx]]
  out
}

# Primeiro valor efetivamente observado conforme a prioridade dos campos.
first_raw <- function(x, columns) {
  out <- rep(NA_character_, nrow(x))
  for (col in intersect(columns, names(x))) {
    z <- trimws(as.character(x[[col]]))
    ok <- !is.na(z) & nzchar(z) & !toupper(z) %in% c("NA", "NULL", "N/A")
    use <- is.na(out) & ok
    out[use] <- z[use]
  }
  out
}

# Conversoes CONSERVADORAS; faixa etaria, "90+", unidades e texto -> NA.
num_exact <- function(x) {
  z <- trimws(as.character(x))
  ok <- !is.na(z) & grepl("^[0-9]+([.][0-9]+)?$", z)
  out <- rep(NA_real_, length(z))
  out[ok] <- as.numeric(z[ok])
  out
}
sex_exact <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(z))
  out[z %in% c("m", "male", "masculino")] <- "M"
  out[z %in% c("f", "female", "feminino")] <- "F"
  out
}

HD <- add_candidate(read_meta("GSE64810"), "GSE64810")
PD <- add_candidate(read_meta("GSE68719"), "GSE68719")

# A nomenclatura do arquivo original do GEO de HD e "age of death"
# (normalizada pelo 10b v2 para age_of_death), nao apenas age_at_death.
HD[, age_raw := first_raw(HD, c("age_of_death", "age_at_death", "age"))]
HD[, sex_raw := first_raw(HD, c("sex", "gender", "sex_at_birth"))]
HD[, pmi_raw := first_raw(HD, c("pmi", "postmortem_interval", "post_mortem_interval"))]
HD[, rin_raw := first_raw(HD, c("rin", "rna_integrity_number"))]
PD[, age_raw := first_raw(PD, c("age_at_death", "age_of_death", "age"))]
PD[, sex_raw := first_raw(PD, c("sex", "gender", "sex_at_birth"))]
PD[, pmi_raw := first_raw(PD, c("pmi", "postmortem_interval", "post_mortem_interval"))]
PD[, rin_raw := first_raw(PD, c("rin", "rna_integrity_number"))]

# Cruzamento um-para-um e rotulo Control em ambas as coortes.
PD_control <- PD[group == "Control"]
HD_control <- HD[group == "Control"]
M <- merge(PD_control[, .(sample_id, group_PD = group, age_PD = age_raw,
                         sex_PD = sex_raw, pmi_PD = pmi_raw, rin_PD = rin_raw)],
           HD_control[, .(sample_id, group_HD = group, age_HD = age_raw,
                         sex_HD = sex_raw, pmi_HD = pmi_raw, rin_HD = rin_raw)],
           by = "sample_id", all.x = TRUE, sort = FALSE)
M[, status_match := fifelse(!is.na(group_HD), "EXACT_CONTROL_IN_BOTH", "NOT_FOUND_IN_HD_CONTROL")]
M[, age_conflict := !is.na(age_PD) & !is.na(age_HD) & num_exact(age_PD) != num_exact(age_HD)]
M[is.na(age_conflict), age_conflict := FALSE]
M[, sex_conflict := !is.na(sex_PD) & !is.na(sex_HD) &
    !is.na(sex_exact(sex_PD)) & !is.na(sex_exact(sex_HD)) & sex_exact(sex_PD) != sex_exact(sex_HD)]
M[is.na(sex_conflict), sex_conflict := FALSE]
WR(M, "10c_PD_HD_shared_controls_crosswalk.tsv")

n_shared <- sum(M$status_match == "EXACT_CONTROL_IN_BOTH")
cat("Controles PD encontrados como controles HD (sample_id exato): ",
    n_shared, "/", nrow(PD_control), "\n", sep = "")
if (n_shared != 44L || nrow(PD_control) != 44L)
  stop("Nao foram confirmados os 44 controles compartilhados. Revise 10c_*crosswalk.tsv")
if (any(M$age_conflict) || any(M$sex_conflict))
  stop("Conflitos de idade/sexo entre registros! Revise crosswalk antes de transferir.")

# Produz extensoes separadas; NA continua NA.
HD_ext <- data.table(sample_id = HD$sample_id,
                     age = num_exact(HD$age_raw),
                     sex = sex_exact(HD$sex_raw),
                     pmi = num_exact(HD$pmi_raw),
                     rin = num_exact(HD$rin_raw))

lookup_hd <- match(PD$sample_id, HD$sample_id)
is_control <- PD$group == "Control"
get_hd <- function(field) HD[[field]][lookup_hd]
choose_source <- function(pd, hd, use_hd) {
  result <- pd
  result[use_hd] <- hd[use_hd]
  result
}
PD_ext <- data.table(sample_id = PD$sample_id,
  age = num_exact(choose_source(PD$age_raw, get_hd("age_raw"), is_control)),
  sex = sex_exact(choose_source(PD$sex_raw, get_hd("sex_raw"), is_control)),
  pmi = num_exact(choose_source(PD$pmi_raw, get_hd("pmi_raw"), is_control)),
  rin = num_exact(choose_source(PD$rin_raw, get_hd("rin_raw"), is_control)))

WR(HD_ext, "10c_GSE64810_covariates_REVIEW_FIRST.tsv")
WR(PD_ext, "10c_GSE68719_covariates_REVIEW_FIRST.tsv")

# Proveniencia por amostra/campo e indicadores de texto nao convertido.
prov_HD <- data.table(dataset = "GSE64810", sample_id = HD$sample_id, group = HD$group,
  age_raw = HD$age_raw, sex_raw = HD$sex_raw, pmi_raw = HD$pmi_raw, rin_raw = HD$rin_raw,
  source_age = "GSE64810 GEO SOFT", source_sex = "GSE64810 GEO SOFT")
prov_PD <- data.table(dataset = "GSE68719", sample_id = PD$sample_id, group = PD$group,
  age_raw = choose_source(PD$age_raw, get_hd("age_raw"), is_control),
  sex_raw = choose_source(PD$sex_raw, get_hd("sex_raw"), is_control),
  pmi_raw = choose_source(PD$pmi_raw, get_hd("pmi_raw"), is_control),
  rin_raw = choose_source(PD$rin_raw, get_hd("rin_raw"), is_control),
  source_age = fifelse(is_control, "GSE64810 GEO SOFT: matching control",
                       "GSE68719 GEO SOFT: PD case"),
  source_sex = fifelse(is_control, "GSE64810 GEO SOFT: matching control",
                       "GSE68719 GEO SOFT: PD case"))
WR(rbindlist(list(prov_HD, prov_PD)), "10c_provenance_raw_values.tsv")

coverage_one <- function(ext, meta, ds) {
  xx <- merge(meta, ext, by = "sample_id", all.x = TRUE, sort = FALSE)
  rbindlist(lapply(c("age", "sex", "pmi", "rin"), function(f) {
    xx[, .(dataset = ds, field = f, n = .N,
            available = sum(!is.na(get(f))),
            missing = sum(is.na(get(f)))), by = group]
  }))
}
cover <- rbindlist(list(coverage_one(HD_ext, HD[, .(sample_id, group)], "GSE64810"),
                        coverage_one(PD_ext, PD[, .(sample_id, group)], "GSE68719")))
WR(cover, "10c_covariate_coverage_by_group.tsv")
print(cover[order(dataset, group, field)])

notes <- data.table(note = c(
  "A associacao dos controles foi por sample_id EXATO com grupo Control em ambas as coortes.",
  "Consulte o registro GEO para confirmar que os IDs designam os mesmos individuos.",
  "O script nao estima nem inventa metadados ausentes.",
  "NA em sex indica ausencia de anotacao; nao inferir sexo por gene ou diagnostico.",
  "Nao copiar extensoes para metadata_extension sem rever proveniencia e plausibilidade.",
  "Nao autorizar modelos ajustados sem cobertura suficiente de ambos os grupos.",
  "As 44 amostras controle de PD e HD sao compartilhadas: modelos por coorte podem ser separados, mas efeitos transdiagnosticos nao sao independentes."
))
WR(notes, "10c_REVIEW_NOTES.tsv")
cat("10c finalizado: candidatos e auditorias em ", OUT, "\n", sep = "")
