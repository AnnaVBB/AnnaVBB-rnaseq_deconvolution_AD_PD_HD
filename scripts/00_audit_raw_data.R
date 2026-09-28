#!/usr/bin/env Rscript

# ==============================================================================
# 00_audit_raw_data.R
#
# Auditoria inicial e reprodutivel dos dados de bulk RNA-seq.
#
# Datasets:
#   GSE53697 -> Alzheimer (AD)
#   GSE64810 -> Huntington (HD)
#   GSE68719 -> Parkinson (PD)
#
# Objetivos:
#   1. Verificar integridade e estrutura dos arquivos de entrada.
#   2. Auditar raw counts sem modificar os dados.
#   3. Auditar GeneIDs e a anotacao GRCh38.p13.
#   4. Quantificar a sobreposicao genica entre AD, HD e PD.
#   5. Registrar genes all-zero e valores fracionarios.
#   6. Produzir uma tabela por GeneID para uso posterior na harmonizacao.
#
# IMPORTANTE:
#   Este script NAO:
#     - arredonda counts;
#     - remove genes;
#     - filtra baixa expressao;
#     - harmoniza identificadores;
#     - normaliza counts;
#     - executa DESeq2.
#
# Outputs:
#   results/audit/00_audit_raw_data.txt
#   results/audit/00_gene_universe_audit.tsv
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Pacotes
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})


# ------------------------------------------------------------------------------
# 2. Caminhos
# ------------------------------------------------------------------------------

file_ad <- "data/raw/GSE53697/GSE53697_RNAseq_AD.txt.gz"

file_hd <- paste0(
  "data/raw/GSE64810/",
  "GSE64810_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)

file_pd <- paste0(
  "data/raw/GSE68719/",
  "GSE68719_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)

file_annotation <- "data/raw/Human.GRCh38.p13.annot.tsv.gz"

output_dir <- "results/audit"

output_report <- file.path(
  output_dir,
  "00_audit_raw_data.txt"
)

output_gene_universe <- file.path(
  output_dir,
  "00_gene_universe_audit.tsv"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 3. Funcoes auxiliares
# ------------------------------------------------------------------------------

section <- function(title) {

  cat("\n")
  cat(paste(rep("=", 78), collapse = ""), "\n")
  cat(title, "\n")
  cat(paste(rep("=", 78), collapse = ""), "\n")
}


check_file <- function(path, label) {

  if (!file.exists(path)) {
    stop(
      "Arquivo ausente [", label, "]: ",
      path,
      call. = FALSE
    )
  }

  info <- file.info(path)

  if (is.na(info$size) || info$size == 0) {
    stop(
      "Arquivo vazio [", label, "]: ",
      path,
      call. = FALSE
    )
  }

  cat(
    sprintf(
      "%-12s OK | %.2f MB | %s\n",
      label,
      info$size / 1024^2,
      path
    )
  )
}


normalize_gene_id <- function(x) {

  x <- trimws(as.character(x))

  x[
    is.na(x) |
    x == "" |
    toupper(x) == "NA"
  ] <- NA_character_

  x
}


clean_annotation_field <- function(x) {

  x <- trimws(as.character(x))

  x[
    is.na(x) |
    x == "" |
    toupper(x) == "NA" |
    x == "-"
  ] <- NA_character_

  x
}


safe_range <- function(x) {

  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(c(NA_real_, NA_real_))
  }

  range(x)
}


audit_counts <- function(count_matrix, dataset) {

  values <- as.vector(count_matrix)

  rg <- safe_range(values)

  data.table(
    dataset = dataset,
    n_genes = nrow(count_matrix),
    n_samples = ncol(count_matrix),
    n_values = length(values),
    n_na = sum(is.na(values)),
    n_negative = sum(values < 0, na.rm = TRUE),
    n_non_integer = sum(
      abs(values - round(values)) > 1e-8,
      na.rm = TRUE
    ),
    n_zero = sum(values == 0, na.rm = TRUE),
    pct_zero = 100 *
      sum(values == 0, na.rm = TRUE) /
      sum(!is.na(values)),
    min_count = rg[1],
    max_count = rg[2],
    total_counts = sum(values, na.rm = TRUE)
  )
}


audit_fractional <- function(count_matrix, dataset) {

  values <- as.vector(count_matrix)

  fractional <- values[
    !is.na(values) &
    abs(values - round(values)) > 1e-8
  ]

  if (length(fractional) == 0L) {

    return(
      data.table(
        dataset = dataset,
        n_fractional = 0L,
        min_fractional = NA_real_,
        median_fractional = NA_real_,
        max_fractional = NA_real_,
        max_distance_to_integer = 0
      )
    )
  }

  data.table(
    dataset = dataset,
    n_fractional = length(fractional),
    min_fractional = min(fractional),
    median_fractional = median(fractional),
    max_fractional = max(fractional),
    max_distance_to_integer =
      max(abs(fractional - round(fractional)))
  )
}


# ------------------------------------------------------------------------------
# 4. Verificar arquivos
# ------------------------------------------------------------------------------

section("1. VERIFICACAO DOS ARQUIVOS DE ENTRADA")

check_file(file_ad, "AD")
check_file(file_hd, "HD")
check_file(file_pd, "PD")
check_file(file_annotation, "ANNOTATION")


# ------------------------------------------------------------------------------
# 5. Carregar arquivos
# ------------------------------------------------------------------------------

section("2. CARREGANDO RAW DATA")

cat("Lendo AD...\n")
ad <- fread(file_ad)

cat("Lendo HD...\n")
hd <- fread(file_hd)

cat("Lendo PD...\n")
pd <- fread(file_pd)

cat("Lendo anotacao GRCh38.p13...\n")
annotation <- fread(file_annotation)

cat("Arquivos carregados com sucesso.\n")


# ------------------------------------------------------------------------------
# 6. Validar estrutura
# ------------------------------------------------------------------------------

section("3. VALIDACAO DA ESTRUTURA")

if (!all(c("GeneID", "GeneSymbol") %in% names(ad))) {
  stop(
    "AD nao possui GeneID e/ou GeneSymbol.",
    call. = FALSE
  )
}

if (!"GeneID" %in% names(hd)) {
  stop(
    "HD nao possui GeneID.",
    call. = FALSE
  )
}

if (!"GeneID" %in% names(pd)) {
  stop(
    "PD nao possui GeneID.",
    call. = FALSE
  )
}


ad_raw_cols <- grep(
  "_raw$",
  names(ad),
  value = TRUE
)

ad_rpkm_cols <- grep(
  "_rpkm$",
  names(ad),
  value = TRUE
)

hd_sample_cols <- setdiff(
  names(hd),
  "GeneID"
)

pd_sample_cols <- setdiff(
  names(pd),
  "GeneID"
)


if (length(ad_raw_cols) != 17L) {
  stop(
    "AD: numero inesperado de colunas raw: ",
    length(ad_raw_cols),
    ". Esperado: 17.",
    call. = FALSE
  )
}

if (length(ad_rpkm_cols) != 17L) {
  stop(
    "AD: numero inesperado de colunas RPKM: ",
    length(ad_rpkm_cols),
    ". Esperado: 17.",
    call. = FALSE
  )
}

if (length(hd_sample_cols) != 69L) {
  stop(
    "HD: numero inesperado de amostras: ",
    length(hd_sample_cols),
    ". Esperado: 69.",
    call. = FALSE
  )
}

if (length(pd_sample_cols) != 73L) {
  stop(
    "PD: numero inesperado de amostras: ",
    length(pd_sample_cols),
    ". Esperado: 73.",
    call. = FALSE
  )
}


file_dimensions <- data.table(
  dataset = c("AD", "HD", "PD"),
  n_genes = c(
    nrow(ad),
    nrow(hd),
    nrow(pd)
  ),
  n_columns_total = c(
    ncol(ad),
    ncol(hd),
    ncol(pd)
  ),
  n_count_samples = c(
    length(ad_raw_cols),
    length(hd_sample_cols),
    length(pd_sample_cols)
  )
)

print(file_dimensions)

cat("\nAD RPKM columns:", length(ad_rpkm_cols), "\n")


# ------------------------------------------------------------------------------
# 7. Construir matrizes temporarias de raw counts
# ------------------------------------------------------------------------------

ad_counts <- as.matrix(
  ad[, ..ad_raw_cols]
)

hd_counts <- as.matrix(
  hd[, ..hd_sample_cols]
)

pd_counts <- as.matrix(
  pd[, ..pd_sample_cols]
)

storage.mode(ad_counts) <- "numeric"
storage.mode(hd_counts) <- "numeric"
storage.mode(pd_counts) <- "numeric"


# ------------------------------------------------------------------------------
# 8. Auditoria dos counts
# ------------------------------------------------------------------------------

section("4. AUDITORIA DOS RAW COUNTS")

count_summary <- rbindlist(
  list(
    audit_counts(ad_counts, "AD"),
    audit_counts(hd_counts, "HD"),
    audit_counts(pd_counts, "PD")
  )
)

print(count_summary)


if (any(count_summary$n_na > 0)) {
  stop(
    "Foram detectados valores NA nos raw counts.",
    call. = FALSE
  )
}

if (any(count_summary$n_negative > 0)) {
  stop(
    "Foram detectados valores negativos nos raw counts.",
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# 9. Valores fracionarios
# ------------------------------------------------------------------------------

section("5. VALORES FRACIONARIOS")

fractional_summary <- rbindlist(
  list(
    audit_fractional(ad_counts, "AD"),
    audit_fractional(hd_counts, "HD"),
    audit_fractional(pd_counts, "PD")
  )
)

print(fractional_summary)

if (any(count_summary$n_non_integer > 0)) {

  cat(
    "\nNOTE: valores fracionarios foram detectados. ",
    "Nenhum arredondamento foi realizado nesta etapa.\n",
    sep = ""
  )
}


# ------------------------------------------------------------------------------
# 10. GeneIDs
# ------------------------------------------------------------------------------

section("6. AUDITORIA DOS GeneIDs")

ad_gene_id <- normalize_gene_id(ad$GeneID)
hd_gene_id <- normalize_gene_id(hd$GeneID)
pd_gene_id <- normalize_gene_id(pd$GeneID)

gene_id_summary <- data.table(
  dataset = c("AD", "HD", "PD"),

  n_rows = c(
    length(ad_gene_id),
    length(hd_gene_id),
    length(pd_gene_id)
  ),

  n_missing = c(
    sum(is.na(ad_gene_id)),
    sum(is.na(hd_gene_id)),
    sum(is.na(pd_gene_id))
  ),

  n_unique = c(
    uniqueN(ad_gene_id[!is.na(ad_gene_id)]),
    uniqueN(hd_gene_id[!is.na(hd_gene_id)]),
    uniqueN(pd_gene_id[!is.na(pd_gene_id)])
  ),

  n_duplicated = c(
    sum(duplicated(ad_gene_id) & !is.na(ad_gene_id)),
    sum(duplicated(hd_gene_id) & !is.na(hd_gene_id)),
    sum(duplicated(pd_gene_id) & !is.na(pd_gene_id))
  )
)

print(gene_id_summary)


if (
  any(gene_id_summary$n_missing > 0) ||
  any(gene_id_summary$n_duplicated > 0)
) {

  stop(
    "GeneIDs ausentes ou duplicados foram encontrados.",
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# 11. Simbolos do AD
# ------------------------------------------------------------------------------

section("7. SIMBOLOS DO AD")

ad_symbol <- clean_annotation_field(
  ad$GeneSymbol
)

ad_symbol_summary <- data.table(
  n_rows = length(ad_symbol),
  n_present = sum(!is.na(ad_symbol)),
  n_missing = sum(is.na(ad_symbol)),
  n_unique = uniqueN(ad_symbol[!is.na(ad_symbol)]),
  n_duplicated = sum(
    duplicated(ad_symbol) &
    !is.na(ad_symbol)
  )
)

print(ad_symbol_summary)


# ------------------------------------------------------------------------------
# 12. Anotacao GRCh38.p13
# ------------------------------------------------------------------------------

section("8. ANOTACAO GRCh38.p13")

required_annotation <- c(
  "GeneID",
  "Symbol",
  "GeneType",
  "EnsemblGeneID"
)

missing_annotation <- setdiff(
  required_annotation,
  names(annotation)
)

if (length(missing_annotation) > 0L) {

  stop(
    "Colunas ausentes na anotacao: ",
    paste(missing_annotation, collapse = ", "),
    call. = FALSE
  )
}


annotation_gene_id <- normalize_gene_id(
  annotation$GeneID
)

annotation_symbol <- clean_annotation_field(
  annotation$Symbol
)

annotation_ensembl <- clean_annotation_field(
  annotation$EnsemblGeneID
)

annotation_gene_type <- clean_annotation_field(
  annotation$GeneType
)


annotation_summary <- data.table(
  n_rows = nrow(annotation),

  n_unique_gene_id =
    uniqueN(
      annotation_gene_id[
        !is.na(annotation_gene_id)
      ]
    ),

  n_duplicated_gene_id =
    sum(
      duplicated(annotation_gene_id) &
      !is.na(annotation_gene_id)
    ),

  n_symbol_present =
    sum(!is.na(annotation_symbol)),

  n_ensembl_present =
    sum(!is.na(annotation_ensembl)),

  n_ensembl_missing =
    sum(is.na(annotation_ensembl)),

  n_gene_type_present =
    sum(!is.na(annotation_gene_type)),

  n_gene_type_missing =
    sum(is.na(annotation_gene_type))
)

print(annotation_summary)


# ------------------------------------------------------------------------------
# 13. Universo genico
# ------------------------------------------------------------------------------

section("9. UNIVERSO GENICO")

genes_ad <- unique(
  ad_gene_id[!is.na(ad_gene_id)]
)

genes_hd <- unique(
  hd_gene_id[!is.na(hd_gene_id)]
)

genes_pd <- unique(
  pd_gene_id[!is.na(pd_gene_id)]
)

annotation_ids <- unique(
  annotation_gene_id[
    !is.na(annotation_gene_id)
  ]
)

all_gene_ids <- sort(
  unique(
    c(
      genes_ad,
      genes_hd,
      genes_pd
    )
  )
)


gene_universe <- data.table(
  gene_id = all_gene_ids
)

gene_universe[
  ,
  in_AD := gene_id %in% genes_ad
]

gene_universe[
  ,
  in_HD := gene_id %in% genes_hd
]

gene_universe[
  ,
  in_PD := gene_id %in% genes_pd
]

gene_universe[
  ,
  in_annotation := gene_id %in% annotation_ids
]


gene_universe[
  ,
  presence_class := fifelse(
    in_AD & in_HD & in_PD,
    "AD_HD_PD",

    fifelse(
      in_AD & in_HD & !in_PD,
      "AD_HD",

      fifelse(
        in_AD & !in_HD & in_PD,
        "AD_PD",

        fifelse(
          !in_AD & in_HD & in_PD,
          "HD_PD",

          fifelse(
            in_AD & !in_HD & !in_PD,
            "AD_only",

            fifelse(
              !in_AD & in_HD & !in_PD,
              "HD_only",

              fifelse(
                !in_AD & !in_HD & in_PD,
                "PD_only",
                "other"
              )
            )
          )
        )
      )
    )
  )
]


overlap_summary <- gene_universe[
  ,
  .N,
  by = presence_class
][
  order(presence_class)
]

setnames(
  overlap_summary,
  "N",
  "n_genes"
)

print(overlap_summary)


cat(
  "\nGenes comuns AD/HD/PD:",
  sum(
    gene_universe$in_AD &
    gene_universe$in_HD &
    gene_universe$in_PD
  ),
  "\n"
)

cat(
  "GeneIDs AD ausentes da anotacao:",
  sum(
    gene_universe$in_AD &
    !gene_universe$in_annotation
  ),
  "\n"
)


# ------------------------------------------------------------------------------
# 14. Genes all-zero
# ------------------------------------------------------------------------------

section("10. GENES ALL-ZERO")

all_zero_summary <- data.table(
  dataset = c("AD", "HD", "PD"),

  n_genes = c(
    nrow(ad_counts),
    nrow(hd_counts),
    nrow(pd_counts)
  ),

  n_all_zero = c(
    sum(rowSums(ad_counts) == 0),
    sum(rowSums(hd_counts) == 0),
    sum(rowSums(pd_counts) == 0)
  )
)

all_zero_summary[
  ,
  pct_all_zero :=
    100 * n_all_zero / n_genes
]

print(all_zero_summary)


# ------------------------------------------------------------------------------
# 15. Salvar tabela do universo genico
# ------------------------------------------------------------------------------

fwrite(
  gene_universe,
  output_gene_universe,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------------------------
# 16. Criar relatorio humano
# ------------------------------------------------------------------------------

capture.output(
  {

    cat("RAW RNA-SEQ DATA AUDIT\n")
    cat("======================\n\n")

    cat(
      "Audit performed before count preparation, ",
      "gene harmonization, filtering, normalization ",
      "or differential expression analysis.\n\n",
      sep = ""
    )


    cat("1. FILE DIMENSIONS\n")
    cat("------------------\n")
    print(file_dimensions)
    cat("\nAD RPKM columns:", length(ad_rpkm_cols), "\n\n")


    cat("2. RAW COUNT SUMMARY\n")
    cat("--------------------\n")
    print(count_summary)
    cat("\n")


    cat("3. FRACTIONAL VALUES\n")
    cat("--------------------\n")
    print(fractional_summary)

    if (any(count_summary$n_non_integer > 0)) {

      cat(
        "\nNOTE: fractional values were detected. ",
        "No rounding was performed during audit.\n"
      )
    }

    cat("\n")


    cat("4. GENE ID SUMMARY\n")
    cat("------------------\n")
    print(gene_id_summary)
    cat("\n")


    cat("5. AD SYMBOL SUMMARY\n")
    cat("--------------------\n")
    print(ad_symbol_summary)
    cat("\n")


    cat("6. GRCh38.p13 ANNOTATION\n")
    cat("------------------------\n")
    print(annotation_summary)
    cat("\n")


    cat("7. GENE UNIVERSE\n")
    cat("----------------\n")
    print(overlap_summary)

    cat(
      "\nCommon AD/HD/PD GeneIDs:",
      sum(
        gene_universe$in_AD &
        gene_universe$in_HD &
        gene_universe$in_PD
      ),
      "\n"
    )

    cat(
      "AD GeneIDs absent from GRCh38.p13 annotation:",
      sum(
        gene_universe$in_AD &
        !gene_universe$in_annotation
      ),
      "\n\n"
    )


    cat("8. ALL-ZERO GENES\n")
    cat("-----------------\n")
    print(all_zero_summary)
    cat("\n")


    cat("9. AUDIT STATUS\n")
    cat("---------------\n")

    cat(
      "PASS: no missing, negative, or duplicated GeneIDs ",
      "were detected in the input matrices.\n"
    )

    cat(
      "PASS: no missing or negative raw-count values ",
      "were detected.\n"
    )

    if (any(count_summary$n_non_integer > 0)) {

      cat(
        "NOTE: fractional raw values were detected and ",
        "must be handled explicitly during count preparation ",
        "before DESeq2.\n"
      )
    }

    cat(
      "\nNo genes were removed, harmonized, filtered ",
      "or normalized during this audit.\n"
    )

  },
  file = output_report
)


# ------------------------------------------------------------------------------
# 17. Final
# ------------------------------------------------------------------------------

section("11. AUDITORIA CONCLUIDA")

cat("Outputs permanentes:\n")
cat("  ", output_report, "\n")
cat("  ", output_gene_universe, "\n\n")

cat(
  "Nenhum gene foi removido, harmonizado, ",
  "filtrado ou normalizado.\n",
  sep = ""
)



