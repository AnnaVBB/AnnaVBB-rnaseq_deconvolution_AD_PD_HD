#!/usr/bin/env Rscript
# ============================================================================
# 10b | Recuperar metadados originais das series GEO (arquivos SOFT locais)
# TCC AD / PD / HD | SOMENTE AUDITORIA, SEM ALTERAR DADOS DO PROJETO
# ============================================================================
# Objetivo:
# - Examinar os arquivos GSE..._family.soft[.gz] ja baixados em data/raw/;
# - Recuperar os campos GEO !Sample_*, com procedencia para cada campo;
# - Fazer crosswalk exato entre GSM e sample_id dos metadados da etapa 05;
# - Entregar EXTENSOES CANDIDATAS, nunca aplicar automaticamente.
#
# SEGURANCA METODOLOGICA:
# - Nao consulta a internet, nao inventa idade/sexo/regiao e nao imputa NAs;
# - Nao transforma 'braak_score' ou gravidade em covariavel de ajuste;
# - Nao confunde regiao constante por estudo com ajuste intraestudo;
# - Nao promove tabela candidata para metadata_extension automaticamente;
# - ROSMAP exige obtencao licita de metadados originais e ID crosswalk;
# - Uma correspondencia por titulo so vale quando titulo e ID sao EXATAMENTE
#   identicos e unicos, registrando a modalidade de match.
#
# Uso: Rscript scripts/phase02_deconvolution/10b_recuperar_covariaveis_GEO_SOFT.R
# Saida: results/deconvolution/Integrated_Comparison/10_metadata/10b_GEO/
# ============================================================================
options(stringsAsFactors = FALSE)
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Requer o pacote R data.table (renv).");
}
suppressPackageStartupMessages(library(data.table))
ROOT <- normalizePath(Sys.getenv("DECONV_ROOT", "."), mustWork = TRUE)
P <- function(...) file.path(ROOT, ...)
OUT <- P("results", "deconvolution", "Integrated_Comparison", "10_metadata", "10b_GEO")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DATASETS <- c("GSE125583", "GSE64810", "GSE68719")

WR <- function(x, name) {
  fwrite(x, file.path(OUT, name), sep = "\t", na = "NA", quote = TRUE)
  invisible(x)
}

clean_value <- function(x) trimws(as.character(x))
clean_field <- function(x) {
  y <- tolower(trimws(x))
  y <- gsub("[^a-z0-9]+", "_", y)
  y <- gsub("^_+|_+$", "", y)
  y
}

# Procurar somente nomes de familia GEO em diretorios esperados.
find_geo_soft <- function(ds) {
  bases <- c(P("data", "raw", ds), P("data", "raw"))
  res <- character()
  for (dir in bases) {
    if (!dir.exists(dir)) next
    paths <- list.files(dir, recursive = TRUE, full.names = TRUE,
                        pattern = "\\.soft(\\.gz)?$", ignore.case = TRUE)
    res <- c(res, paths[grepl(ds, basename(paths), fixed = TRUE) &
                           grepl("family", basename(paths), ignore.case = TRUE)])
  }
  unique(normalizePath(res, winslash = "/", mustWork = TRUE))
}

# Leitura progressiva: ignora tabelas de plataformas e expressoes do SOFT.
read_soft_fields <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  current_gsm <- NA_character_
  rows <- list()
  pos <- 0L
  repeat {
    block <- readLines(con, n = 20000L, warn = FALSE)
    if (!length(block)) break
    for (line in block) {
      if (startsWith(line, "^SAMPLE = ")) {
        current_gsm <- clean_value(sub("^\\^SAMPLE = ", "", line))
        next
      }
      if (is.na(current_gsm) || !startsWith(line, "!Sample_")) next
      pieces <- regexpr(" = ", line, fixed = TRUE)
      if (pieces[1L] < 1L) next
      header <- substr(line, 9L, pieces[1L] - 1L)
      val <- clean_value(substr(line, pieces[1L] + 3L, nchar(line)))
      if (!nzchar(val)) next
      src <- header
      if (grepl("^characteristics_ch[0-9]+$", header)) {
        sep <- regexpr(":", val, fixed = TRUE)
        if (sep[1L] < 1L) next
        key <- clean_field(substr(val, 1L, sep[1L] - 1L))
        val <- clean_value(substr(val, sep[1L] + 1L, nchar(val)))
        src <- paste0("characteristics:", key)
      } else if (header == "title") {
        key <- "geo_title"
      } else if (grepl("^source_name", header)) {
        key <- "geo_source_name"
      } else if (header == "description") {
        key <- "geo_description"
      } else {
        # Registra o campo original, mas nao interpreta arbitrariamente.
        key <- paste0("geo_", clean_field(header))
      }
      if (!nzchar(key) || !nzchar(val)) next
      pos <- pos + 1L
      rows[[pos]] <- list(gsm = current_gsm, field = key, value = val,
                          source = src, soft_file = path)
    }
  }
  if (!length(rows)) return(data.table(gsm = character(), field = character(),
                                      value = character(), source = character(),
                                      soft_file = character()))
  rbindlist(rows, fill = TRUE)
}

# As fontes GEO podem ter caracteristicas repetidas; nao remove a proveniencia.
# Valores diferentes para o MESMO campo e GSM sao sinalizados como conflito,
# em vez de escolher silenciosamente um deles.
canonicalize <- function(long) {
  z <- unique(long[, .(gsm, field, value)])
  conf <- z[, .(n_values = uniqueN(value),
                values = paste(unique(value), collapse = " | ")), by = .(gsm, field)]
  duplicates <- conf[n_values > 1L]
  safe <- conf[n_values == 1L, .(gsm, field, value = values)]
  wide <- dcast(safe, gsm ~ field, value.var = "value")
  list(wide = wide, conflicts = duplicates)
}

# Correspondencia de titulo aceita SOMENTE titulo exatamente igual a sample_id
# e unico nos GEO, nao tenta adivinhar a relacao por doenca/ordem.
match_samples <- function(base, wide) {
  stopifnot("sample_id" %in% names(base), "gsm" %in% names(wide))
  base <- copy(base[, .(sample_id, group)])
  base[, sample_id := clean_value(sample_id)]
  gsm <- wide[, .(gsm)]
  gsm[, geo_sample_id := gsm]
  direct <- merge(base, gsm[, .(sample_id = gsm, geo_sample_id)],
                  by = "sample_id", all.x = TRUE, sort = FALSE)
  direct[, match_method := fifelse(!is.na(geo_sample_id), "GSM_EXACT", NA_character_)]
  if ("geo_title" %in% names(wide)) {
    titles <- wide[!is.na(geo_title) & nzchar(geo_title), .(geo_title, gsm)]
    titles <- titles[, .(n_gsm = .N, gsm = first(gsm)), by = geo_title][n_gsm == 1L]
    tmatch <- merge(direct[is.na(geo_sample_id), .(sample_id)], titles,
                    by.x = "sample_id", by.y = "geo_title", all.x = TRUE, sort = FALSE)
    if (nrow(tmatch)) {
      direct[match(tmatch$sample_id, direct$sample_id),
             `:=`(geo_sample_id = tmatch$gsm,
                  match_method = fifelse(!is.na(tmatch$gsm), "TITLE_EXACT_UNIQUE", NA_character_))]
    }
  }
  direct[is.na(match_method), match_method := "NOT_MATCHED"]
  if (anyDuplicated(direct[!is.na(geo_sample_id)]$geo_sample_id)) {
    stop("Mais de um sample_id foi associado ao mesmo GSM: verificar crosswalk manualmente")
  }
  direct
}

manifest <- list()
for (ds in DATASETS) {
  cat("\n10b | ", ds, "\n", sep = "")
  basefile <- P("data", "processed", "deconvolution", ds, "metadata_deconvolution.tsv")
  if (!file.exists(basefile)) {
    manifest[[ds]] <- data.table(dataset = ds, status = "BASE_METADATA_MISSING",
                                 soft_file = NA_character_, n_matched = 0L, n_total = NA_integer_)
    next
  }
  base <- fread(basefile, na.strings = c("NA", ""))
  if (!all(c("sample_id", "group") %in% names(base)) || anyDuplicated(base$sample_id)) {
    stop(ds, ": arquivo de etapa 05 nao tem chaves validas")
  }
  files <- find_geo_soft(ds)
  if (!length(files)) {
    manifest[[ds]] <- data.table(dataset = ds, status = "SOFT_NOT_FOUND",
                                 soft_file = NA_character_, n_matched = 0L, n_total = nrow(base))
    WR(base[, .(sample_id, group)], paste0("10b_", ds, "_need_crosswalk.tsv"))
    cat("  SOFT ausente em data/raw/. Baixe o arquivo da serie GEO para data/raw/", ds, "/\n", sep = "")
    next
  }
  if (length(files) > 1L) {
    # Preferir o primeiro por caminho ordenado; nunca misturar versoes sem registro.
    files <- sort(files)
    cat("  Mais de um SOFT encontrado. Selecionado: ", files[1L], "\n", sep = "")
  }
  file <- files[1L]
  long <- read_soft_fields(file)
  if (!nrow(long)) {
    manifest[[ds]] <- data.table(dataset = ds, status = "SOFT_WITHOUT_FIELDS",
                                 soft_file = file, n_matched = 0L, n_total = nrow(base))
    next
  }
  WR(long, paste0("10b_", ds, "_raw_fields_with_provenance.tsv"))
  parsed <- canonicalize(long)
  WR(parsed$conflicts, paste0("10b_", ds, "_field_conflicts.tsv"))
  wide <- parsed$wide
  map <- match_samples(base, wide)
  WR(map, paste0("10b_", ds, "_sample_crosswalk.tsv"))
  matched <- merge(map, wide, by.x = "geo_sample_id", by.y = "gsm", all.x = TRUE, sort = FALSE)
  candidate_keys <- c("age", "age_at_death", "sex", "gender", "pmi", "postmortem_interval",
                       "rin", "rna_integrity_number", "batch", "sequencing_batch",
                       "brain_region", "region", "tissue_subtype", "tissue", "braak_score",
                       "apoe_status", "geo_source_name")
  available <- intersect(candidate_keys, names(matched))
  cand <- matched[, c("sample_id", available), with = FALSE]
  if ("geo_source_name" %in% names(cand)) {
    setnames(cand, "geo_source_name", "geo_source_name_provisional")
  }
  # Nao promover GEO title/diagnosis para covariavel; preservar diagnostico
  # apenas no arquivo bruto para revisao de consistencia com group.
  WR(cand, paste0("10b_", ds, "_candidate_extension_REVIEW_FIRST.tsv"))
  summary_fields <- if (length(available)) {
    rbindlist(lapply(available, function(nm) {
      data.table(dataset = ds, column = nm,
                 n_available = sum(!is.na(matched[[nm]]) & nzchar(matched[[nm]])),
                 n_missing = sum(is.na(matched[[nm]]) | !nzchar(matched[[nm]])),
                 n_distinct = uniqueN(matched[[nm]][!is.na(matched[[nm]])]))
    }))
  } else data.table(dataset = character(), column = character(), n_available = integer(),
                    n_missing = integer(), n_distinct = integer())
  WR(summary_fields, paste0("10b_", ds, "_candidate_coverage.tsv"))
  nmatch <- sum(map$match_method != "NOT_MATCHED")
  manifest[[ds]] <- data.table(dataset = ds,
      status = if (nmatch == nrow(base) && nrow(parsed$conflicts) == 0L) "CROSSWALK_COMPLETE_REVIEW_REQUIRED" else "REVIEW_MATCHES_OR_CONFLICTS",
      soft_file = file, n_matched = nmatch, n_total = nrow(base),
      n_fields = ncol(wide)-1L, n_field_conflicts = nrow(parsed$conflicts))
  cat("  Amostras correspondentes: ", nmatch, "/", nrow(base),
      " | campos disponiveis: ", paste(available, collapse = ", "), "\n", sep = "")
}
WR(rbindlist(manifest, fill = TRUE), "10b_source_manifest.tsv")
WR(data.table(
  note = c(
    "Nenhuma extensao foi instalada automaticamente em data/processed/deconvolution/metadata_extension/.",
    "Revise valores e unidades dos campos; por exemplo age pode ser texto ou faixa etaria.",
    "Tissue ou source_name nao e automaticamente equivalente a brain_region.",
    "BRAAK, estagio e gravidade de doenca NAO devem entrar como ajuste sem justificativa causal especifica.",
    "ROSMAP nao e GEO: obtenha dados clinicos/biospecime sob os termos de acesso autorizados.",
    "Somente depois de conferencia de IDs, diagnostico, valores e origem copie colunas aprovadas para metadata_extension/.",
    "Rode novamente o script 10 apos incorporar extensoes verificadas; 11 continua bloqueado ate aprovacao."
  )), "10b_REVIEW_NOTES.tsv")
cat("\n10b CONCLUIDO (extracao exploratoria, nao aprova ajuste). Saidas: ", OUT, "\n", sep = "")
