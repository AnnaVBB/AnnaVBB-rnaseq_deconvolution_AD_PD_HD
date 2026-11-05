#!/usr/bin/env Rscript

# ==============================================================================
# 03_harmonize_genes_GSE125583.R
#
# FASE 01b — HARMONIZAÇÃO E ANOTAÇÃO DO UNIVERSO GÊNICO
#
# Objetivos:
#   1) comparar formalmente os GeneIDs brutos entre:
#        - GSE125583 (AD candidato)
#        - GSE64810  (HD)
#        - GSE68719  (PD)
#   2) confirmar interseções, união, igualdade dos conjuntos e igualdade da ordem;
#   3) integrar a anotação Human.GRCh38.p13.annot.tsv.gz;
#   4) criar um gene master para uso nas etapas posteriores.
#
# IMPORTANTE:
#   - este script NÃO altera counts;
#   - NÃO filtra genes por expressão;
#   - NÃO define ainda o universo testável para DESeq2;
#   - NÃO executa QC nem DESeq2.
#
# A etapa 03 descreve o universo gênico BRUTO.
# O universo retido para expressão diferencial será definido na etapa 04.
#
# Entradas:
#   results/candidate_GSE125583/counts/counts_GSE125583.tsv
#   data/raw/GSE64810/GSE64810_raw_counts_GRCh38.p13_NCBI.tsv.gz
#   data/raw/GSE68719/GSE68719_raw_counts_GRCh38.p13_NCBI.tsv.gz
#   data/raw/Human.GRCh38.p13.annot.tsv.gz
#
# Saídas:
#   results/candidate_GSE125583/harmonization/
#     gene_master_GSE125583_HD_PD.tsv
#     gene_presence_audit_GSE125583_HD_PD.tsv
#     gene_annotation_missingness.tsv
#     03_harmonize_genes_GSE125583.txt
# ==============================================================================


# ==============================================================================
# 1. PACOTES
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
# 2. CAMINHOS
# ==============================================================================

ad_counts_file <- file.path(
  "results", "candidate_GSE125583", "counts", "counts_GSE125583.tsv"
)

hd_counts_file <- file.path(
  "data", "raw", "GSE64810",
  "GSE64810_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)

pd_counts_file <- file.path(
  "data", "raw", "GSE68719",
  "GSE68719_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)

annotation_file <- file.path(
  "data", "raw", "Human.GRCh38.p13.annot.tsv.gz"
)

output_dir <- file.path(
  "results", "candidate_GSE125583", "harmonization"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gene_master_file <- file.path(
  output_dir, "gene_master_GSE125583_HD_PD.tsv"
)

presence_audit_file <- file.path(
  output_dir, "gene_presence_audit_GSE125583_HD_PD.tsv"
)

annotation_missingness_file <- file.path(
  output_dir, "gene_annotation_missingness.tsv"
)

audit_file <- file.path(
  output_dir, "03_harmonize_genes_GSE125583.txt"
)

EXPECTED_N_GENES <- 39376L


# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================

# Interrompe se um arquivo não existir ou estiver vazio.
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


# Padroniza valores textuais ausentes.
normalize_missing_character <- function(x) {
  x <- trimws(as.character(x))

  x[x %in% c("", "NA", "N/A", "na", "n/a")] <- NA_character_

  x
}


# Verifica GeneID ausente ou duplicado.
validate_gene_ids <- function(gene_ids, dataset_name) {
  gene_ids <- as.character(gene_ids)

  missing_gene_id <- is.na(gene_ids) | trimws(gene_ids) == ""

  if (any(missing_gene_id)) {
    stop(dataset_name, ": existem GeneIDs ausentes.")
  }

  if (anyDuplicated(gene_ids)) {
    stop(dataset_name, ": existem GeneIDs duplicados.")
  }

  invisible(TRUE)
}


# Lê apenas GeneID, evitando carregar todas as colunas de counts.
read_gene_ids <- function(path, dataset_name) {
  dt <- fread(
    path,
    select = "GeneID",
    showProgress = FALSE
  )

  if (!"GeneID" %in% names(dt)) {
    stop(dataset_name, ": coluna GeneID não encontrada.")
  }

  gene_ids <- as.character(dt$GeneID)

  validate_gene_ids(gene_ids, dataset_name)

  gene_ids
}


# Formatação para relatório em padrão pt-BR.
format_int <- function(x) {
  format(
    x,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


# Conta campos de anotação realmente preenchidos.
count_non_missing <- function(x) {
  sum(
    !is.na(x) &
      trimws(as.character(x)) != ""
  )
}


# ==============================================================================
# 4. VALIDAR ENTRADAS
# ==============================================================================

input_files <- c(
  ad_counts_file,
  hd_counts_file,
  pd_counts_file,
  annotation_file
)

invisible(lapply(input_files, validate_file))


# ==============================================================================
# 5. LER GeneIDs DOS TRÊS DATASETS
# ==============================================================================

ids_ad <- read_gene_ids(ad_counts_file, "GSE125583")
ids_hd <- read_gene_ids(hd_counts_file, "GSE64810")
ids_pd <- read_gene_ids(pd_counts_file, "GSE68719")


# As três matrizes já foram auditadas anteriormente como contendo 39.376 genes.
# Aqui esse número funciona como trava para detectar mudanças inesperadas.
for (item in list(
  GSE125583 = ids_ad,
  GSE64810 = ids_hd,
  GSE68719 = ids_pd
)) {
  if (length(item) != EXPECTED_N_GENES) {
    stop(
      "Número inesperado de GeneIDs. Esperado: ",
      EXPECTED_N_GENES,
      "; observado: ",
      length(item)
    )
  }
}


# ==============================================================================
# 6. COMPARAR UNIVERSOS BRUTOS
# ==============================================================================

intersection_ad_hd <- intersect(ids_ad, ids_hd)
intersection_ad_pd <- intersect(ids_ad, ids_pd)
intersection_hd_pd <- intersect(ids_hd, ids_pd)

intersection_all <- Reduce(
  intersect,
  list(ids_ad, ids_hd, ids_pd)
)

union_all <- Reduce(
  union,
  list(ids_ad, ids_hd, ids_pd)
)

# setequal() ignora ordem.
setequal_ad_hd <- setequal(ids_ad, ids_hd)
setequal_ad_pd <- setequal(ids_ad, ids_pd)
setequal_hd_pd <- setequal(ids_hd, ids_pd)

# identical() exige mesmo conteúdo E mesma ordem.
identical_order_ad_hd <- identical(ids_ad, ids_hd)
identical_order_ad_pd <- identical(ids_ad, ids_pd)
identical_order_hd_pd <- identical(ids_hd, ids_pd)


# ==============================================================================
# 7. TABELA DE PRESENÇA GÊNICA
# ==============================================================================

gene_presence <- data.table(
  GeneID = union_all
)

gene_presence[
  ,
  present_GSE125583 := GeneID %in% ids_ad
]

gene_presence[
  ,
  present_GSE64810 := GeneID %in% ids_hd
]

gene_presence[
  ,
  present_GSE68719 := GeneID %in% ids_pd
]

gene_presence[
  ,
  present_all_three :=
    present_GSE125583 &
    present_GSE64810 &
    present_GSE68719
]


# Classifica o padrão de presença.
# Neste caso esperamos AD+HD+PD para todos os genes, mas mantemos a lógica
# completa para auditoria e possível reutilização futura.
gene_presence[
  ,
  presence_pattern := fcase(
    present_GSE125583 & present_GSE64810 & present_GSE68719, "AD+HD+PD",
    present_GSE125583 & present_GSE64810 & !present_GSE68719, "AD+HD",
    present_GSE125583 & !present_GSE64810 & present_GSE68719, "AD+PD",
    !present_GSE125583 & present_GSE64810 & present_GSE68719, "HD+PD",
    present_GSE125583 & !present_GSE64810 & !present_GSE68719, "AD_only",
    !present_GSE125583 & present_GSE64810 & !present_GSE68719, "HD_only",
    !present_GSE125583 & !present_GSE64810 & present_GSE68719, "PD_only",
    default = "none"
  )
]

presence_summary <- gene_presence[
  ,
  .N,
  by = presence_pattern
][
  order(-N)
]


# ==============================================================================
# 8. LER E VALIDAR A ANOTAÇÃO GRCh38.p13
# ==============================================================================

annotation <- fread(
  annotation_file,
  showProgress = FALSE
)

if (!"GeneID" %in% names(annotation)) {
  stop("Arquivo de anotação sem coluna GeneID.")
}

annotation[
  ,
  GeneID := as.character(GeneID)
]

validate_gene_ids(
  annotation$GeneID,
  "Human.GRCh38.p13.annot"
)

expected_annotation_columns <- c(
  "GeneID",
  "Symbol",
  "Description",
  "Synonyms",
  "GeneType",
  "EnsemblGeneID"
)

missing_annotation_columns <- setdiff(
  expected_annotation_columns,
  names(annotation)
)

if (length(missing_annotation_columns) > 0) {
  stop(
    "Arquivo de anotação não contém as colunas esperadas: ",
    paste(missing_annotation_columns, collapse = ", ")
  )
}

annotation_text_columns <- c(
  "Symbol",
  "Description",
  "Synonyms",
  "GeneType",
  "EnsemblGeneID"
)

for (col_name in annotation_text_columns) {
  set(
    annotation,
    j = col_name,
    value = normalize_missing_character(annotation[[col_name]])
  )
}


# ==============================================================================
# 9. COBERTURA DA ANOTAÇÃO
# ==============================================================================

annotation_ids <- annotation$GeneID

genes_without_annotation_row <- setdiff(
  union_all,
  annotation_ids
)

annotation_not_in_datasets <- setdiff(
  annotation_ids,
  union_all
)


# ==============================================================================
# 10. CONSTRUIR O GENE MASTER
# ==============================================================================
#
# O left join parte de gene_presence para garantir que nenhum GeneID dos
# datasets seja perdido por ausência de anotação.

gene_master <- merge(
  gene_presence,
  annotation[
    ,
    ..expected_annotation_columns
  ],
  by = "GeneID",
  all.x = TRUE,
  sort = FALSE
)


# Restaurar uma ordem canônica.
# Como GSE125583 faz parte do conjunto principal candidato, usamos sua ordem.
canonical_order <- c(
  ids_ad,
  setdiff(gene_master$GeneID, ids_ad)
)

gene_master[
  ,
  order_index := match(GeneID, canonical_order)
]

setorder(
  gene_master,
  order_index
)

gene_master[
  ,
  order_index := NULL
]


# ==============================================================================
# 11. VALIDAR O GENE MASTER
# ==============================================================================

if (nrow(gene_master) != length(union_all)) {
  stop(
    "Gene master não possui o mesmo número de genes da união."
  )
}

if (anyDuplicated(gene_master$GeneID)) {
  stop(
    "Gene master contém GeneIDs duplicados."
  )
}

if (!setequal(gene_master$GeneID, union_all)) {
  stop(
    "Gene master não preservou exatamente o universo gênico."
  )
}


# ==============================================================================
# 12. AUDITAR MISSINGNESS DE ANOTAÇÃO
# ==============================================================================

annotation_missingness <- data.table(
  field = c(
    "annotation_row",
    "Symbol",
    "Description",
    "Synonyms",
    "GeneType",
    "EnsemblGeneID"
  ),

  n_present = c(
    sum(gene_master$GeneID %in% annotation_ids),
    count_non_missing(gene_master$Symbol),
    count_non_missing(gene_master$Description),
    count_non_missing(gene_master$Synonyms),
    count_non_missing(gene_master$GeneType),
    count_non_missing(gene_master$EnsemblGeneID)
  )
)

annotation_missingness[
  ,
  n_missing := nrow(gene_master) - n_present
]

annotation_missingness[
  ,
  percent_present := round(
    100 * n_present / nrow(gene_master),
    2
  )
]


# ==============================================================================
# 13. SALVAR SAÍDAS
# ==============================================================================

fwrite(
  gene_master,
  gene_master_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  gene_presence,
  presence_audit_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  annotation_missingness,
  annotation_missingness_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 14. RELATÓRIO TEXTUAL
# ==============================================================================

sink(audit_file)

cat("============================================================\n")
cat("03 — HARMONIZAÇÃO GÊNICA GSE125583 + HD + PD\n")
cat("============================================================\n\n")


cat("------------------------------------------------------------\n")
cat("1. UNIVERSOS GÊNICOS BRUTOS\n")
cat("------------------------------------------------------------\n")

cat("GSE125583: ", format_int(length(ids_ad)), " genes\n", sep = "")
cat("GSE64810: ", format_int(length(ids_hd)), " genes\n", sep = "")
cat("GSE68719: ", format_int(length(ids_pd)), " genes\n\n", sep = "")


cat("------------------------------------------------------------\n")
cat("2. INTERSEÇÕES BRUTAS\n")
cat("------------------------------------------------------------\n")

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
  "Interseção tripla: ",
  format_int(length(intersection_all)),
  "\n",
  sep = ""
)

cat(
  "União dos três datasets: ",
  format_int(length(union_all)),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("3. IGUALDADE DOS CONJUNTOS\n")
cat("------------------------------------------------------------\n")

cat(
  "GSE125583 == GSE64810 (setequal): ",
  setequal_ad_hd,
  "\n",
  sep = ""
)

cat(
  "GSE125583 == GSE68719 (setequal): ",
  setequal_ad_pd,
  "\n",
  sep = ""
)

cat(
  "GSE64810 == GSE68719 (setequal): ",
  setequal_hd_pd,
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("4. IGUALDADE DA ORDEM DOS GeneIDs\n")
cat("------------------------------------------------------------\n")

cat(
  "GSE125583 == GSE64810 (identical): ",
  identical_order_ad_hd,
  "\n",
  sep = ""
)

cat(
  "GSE125583 == GSE68719 (identical): ",
  identical_order_ad_pd,
  "\n",
  sep = ""
)

cat(
  "GSE64810 == GSE68719 (identical): ",
  identical_order_hd_pd,
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("5. PADRÕES DE PRESENÇA\n")
cat("------------------------------------------------------------\n")

print(presence_summary)
cat("\n")


cat("------------------------------------------------------------\n")
cat("6. ANOTAÇÃO GRCh38.p13\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes na anotação: ",
  format_int(nrow(annotation)),
  "\n",
  sep = ""
)

cat(
  "Genes dos datasets sem linha na anotação: ",
  format_int(length(genes_without_annotation_row)),
  "\n",
  sep = ""
)

cat(
  "Genes da anotação ausentes nos três datasets: ",
  format_int(length(annotation_not_in_datasets)),
  "\n\n",
  sep = ""
)

print(annotation_missingness)
cat("\n")


cat("------------------------------------------------------------\n")
cat("7. GENE MASTER\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes no gene master: ",
  format_int(nrow(gene_master)),
  "\n",
  sep = ""
)

cat(
  "GeneIDs duplicados no gene master: ",
  format_int(anyDuplicated(gene_master$GeneID)),
  "\n",
  sep = ""
)

cat(
  "Mesmo universo da união: ",
  setequal(gene_master$GeneID, union_all),
  "\n",
  sep = ""
)

cat(
  "Genes presentes nos três datasets: ",
  format_int(sum(gene_master$present_all_three)),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("8. INTERPRETAÇÃO DA ETAPA\n")
cat("------------------------------------------------------------\n")

if (setequal_ad_hd && setequal_ad_pd && setequal_hd_pd) {
  cat(
    "- GSE125583, GSE64810 e GSE68719 possuem exatamente o mesmo conjunto bruto de GeneIDs.\n"
  )
} else {
  cat(
    "- Os três datasets NÃO possuem universos brutos de GeneIDs idênticos.\n"
  )
}

if (
  identical_order_ad_hd &&
  identical_order_ad_pd &&
  identical_order_hd_pd
) {
  cat(
    "- Os GeneIDs também aparecem na mesma ordem nos três arquivos.\n"
  )
} else {
  cat(
    "- Os conjuntos são comparáveis, mas a ordem dos GeneIDs não é idêntica em todos os arquivos.\n"
  )
}

cat("- Nenhuma matriz de counts foi alterada.\n")
cat("- Nenhum gene foi removido por expressão.\n")
cat("- O gene master registra presença e anotação do universo bruto.\n")
cat(
  "- O universo testável para DESeq2 será definido somente após o pré-filtro da etapa 04.\n"
)
cat(
  "- Ter o mesmo universo bruto não garante que os mesmos genes serão retidos após o filtro de expressão.\n"
)

cat("\n")
cat("============================================================\n")
cat("03_harmonize_genes_GSE125583.R CONCLUÍDO\n")
cat("============================================================\n")

sink()


# ==============================================================================
# 15. RESUMO NO TERMINAL
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("03 — HARMONIZAÇÃO GÊNICA CONCLUÍDA\n")
cat("============================================================\n")

cat("GSE125583: ", format_int(length(ids_ad)), " genes\n", sep = "")
cat("GSE64810: ", format_int(length(ids_hd)), " genes\n", sep = "")
cat("GSE68719: ", format_int(length(ids_pd)), " genes\n", sep = "")

cat(
  "Interseção tripla bruta: ",
  format_int(length(intersection_all)),
  "\n",
  sep = ""
)

cat(
  "União bruta: ",
  format_int(length(union_all)),
  "\n",
  sep = ""
)

cat(
  "Universos exatamente iguais: ",
  setequal_ad_hd && setequal_ad_pd && setequal_hd_pd,
  "\n",
  sep = ""
)

cat(
  "Ordem dos GeneIDs idêntica nos três: ",
  identical_order_ad_hd &&
    identical_order_ad_pd &&
    identical_order_hd_pd,
  "\n",
  sep = ""
)

cat("Gene master: ", gene_master_file, "\n", sep = "")
cat("Auditoria de presença: ", presence_audit_file, "\n", sep = "")
cat(
  "Auditoria de anotação: ",
  annotation_missingness_file,
  "\n",
  sep = ""
)
cat("Relatório: ", audit_file, "\n", sep = "")

cat("============================================================\n")
