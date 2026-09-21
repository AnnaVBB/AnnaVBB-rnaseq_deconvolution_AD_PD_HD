#!/usr/bin/env Rscript

# ============================================================
# 02a_prepare_AD.R
#
# Prepara a matriz de expressão do GSE53697 (Alzheimer) para a análise diferencial da Fase 1.
#
# Entrada principal: matriz processada disponibilizada pelos autores contendo GeneID, GeneSymbol, *_raw e *_rpkm
#
# Etapas:
#   1. Selecionar as 17 colunas RPKM
#   2. Remover genes sem expressão em todas as amostras
#   3. Aplicar log2(RPKM + 1)
#   4. Mapear GeneID (Entrez) -> Ensembl
#   5. Remover versão dos IDs Ensembl, se houver
#   6. Preservar GeneID e GeneSymbol originais
#   7. Salvar matriz preparada
#
# Uso:
# Rscript scripts/02a_prepare_AD.R \
#   data/raw/GSE53697/GSE53697_RNAseq_AD.txt.gz \
#   data/raw/GSE53697/Human.GRCh38.p13.annot.tsv.gz \
#   data/processed/GSE53697_expression_prepared.csv
# ============================================================


# ------------------------------------------------------------
# 1. Pacotes
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})


# ------------------------------------------------------------
# 2. Argumentos
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "\nUso:\n\n",
      "Rscript scripts/02a_prepare_AD.R ",
      "<expression_file> ",
      "<annotation_file> ",
      "<output_file>\n"
    )
  )
}

expression_file <- args[1]
annotation_file <- args[2]
output_file <- args[3]


# ------------------------------------------------------------
# 3. Verificar arquivos
# ------------------------------------------------------------

if (!file.exists(expression_file)) {
  stop("Arquivo de expressão não encontrado: ", expression_file)
}

if (!file.exists(annotation_file)) {
  stop("Arquivo de anotação não encontrado: ", annotation_file)
}

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Carregar matriz original
# ------------------------------------------------------------

message("Carregando matriz de expressão...")

expr <- fread(
  expression_file,
  data.table = FALSE
)

required_columns <- c(
  "GeneID",
  "GeneSymbol"
)

missing_columns <- setdiff(
  required_columns,
  colnames(expr)
)

if (length(missing_columns) > 0) {
  stop(
    "Colunas obrigatórias ausentes: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 5. Identificar colunas RPKM
# ------------------------------------------------------------

rpkm_cols <- grep(
  "_rpkm$",
  colnames(expr),
  value = TRUE
)

if (length(rpkm_cols) != 17) {
  stop(
    "Esperadas 17 colunas RPKM, mas foram encontradas ",
    length(rpkm_cols),
    "."
  )
}

message(
  "Colunas RPKM encontradas: ",
  length(rpkm_cols)
)


# ------------------------------------------------------------
# 6. Construir matriz RPKM
# ------------------------------------------------------------

rpkm <- as.matrix(
  expr[, rpkm_cols, drop = FALSE]
)

storage.mode(rpkm) <- "numeric"

if (anyNA(rpkm)) {
  stop("Foram encontrados valores NA na matriz RPKM.")
}

if (any(rpkm < 0)) {
  stop("Foram encontrados valores negativos na matriz RPKM.")
}

n_genes_initial <- nrow(rpkm)


# ------------------------------------------------------------
# 7. Remover genes all-zero
# ------------------------------------------------------------

keep <- rowSums(rpkm > 0) > 0

n_all_zero <- sum(!keep)

message(
  "Genes sem expressão em todas as amostras: ",
  n_all_zero
)

expr_filtered <- expr[keep, , drop = FALSE]
rpkm_filtered <- rpkm[keep, , drop = FALSE]

n_genes_filtered <- nrow(rpkm_filtered)


# ------------------------------------------------------------
# 8. Transformação log2(RPKM + 1)
# ------------------------------------------------------------

log_expr <- log2(
  rpkm_filtered + 1
)

colnames(log_expr) <- rpkm_cols


# ------------------------------------------------------------
# 9. Converter nomes das amostras
# C1_rpkm  -> C1_raw
# A9_rpkm  -> A9_raw
#
# Isso faz a matriz preparada usar os mesmos sample_id
# presentes em metadata/GSE53697_metadata.csv.
# ------------------------------------------------------------

sample_ids <- sub(
  "_rpkm$",
  "_raw",
  colnames(log_expr)
)

colnames(log_expr) <- sample_ids


# ------------------------------------------------------------
# 10. Carregar anotação NCBI
# ------------------------------------------------------------

message("Carregando anotação...")

annotation <- fread(
  annotation_file,
  data.table = FALSE
)

annotation_required <- c(
  "GeneID",
  "Symbol",
  "EnsemblGeneID"
)

missing_annotation <- setdiff(
  annotation_required,
  colnames(annotation)
)

if (length(missing_annotation) > 0) {
  stop(
    "Colunas ausentes na anotação: ",
    paste(missing_annotation, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 11. Padronizar GeneID para junção
# ------------------------------------------------------------

expr_filtered$GeneID <- as.character(
  expr_filtered$GeneID
)

annotation$GeneID <- as.character(
  annotation$GeneID
)


# ------------------------------------------------------------
# 12. Selecionar anotação necessária
# ------------------------------------------------------------

annotation_small <- annotation %>%
  select(
    GeneID,
    annotation_symbol = Symbol,
    EnsemblGeneID
  ) %>%
  distinct()


# ------------------------------------------------------------
# 13. Juntar expressão e anotação
# ------------------------------------------------------------

gene_info <- expr_filtered %>%
  select(
    GeneID,
    GeneSymbol
  ) %>%
  left_join(
    annotation_small,
    by = "GeneID"
  )


# ------------------------------------------------------------
# 14. Padronizar Ensembl
#
# Alguns arquivos podem conter: # ENSG00000123456.7; O identificador harmonizado será: ENSG00000123456
# ------------------------------------------------------------

gene_info$gene_id <- sub(
  "\\..*$",
  "",
  gene_info$EnsemblGeneID
)

gene_info$gene_id[
  gene_info$gene_id == ""
] <- NA_character_


# ------------------------------------------------------------
# 15. Escolher símbolo
# Prioridade: GeneSymbol fornecido pelo estudo e Symbol da anotação NCBI
# ------------------------------------------------------------

gene_info$symbol <- ifelse(
  !is.na(gene_info$GeneSymbol) &
    gene_info$GeneSymbol != "",
  gene_info$GeneSymbol,
  gene_info$annotation_symbol
)


# ------------------------------------------------------------
# 16. Montar tabela final
# ------------------------------------------------------------

expression_df <- as.data.frame(
  log_expr,
  check.names = FALSE
)

prepared <- bind_cols(
  gene_info %>%
    transmute(
      gene_id = gene_id,
      original_id = GeneID,
      symbol = symbol
    ),
  expression_df
)


# ------------------------------------------------------------
# 17. Auditoria de IDs
# ------------------------------------------------------------

n_ensembl <- sum(
  !is.na(prepared$gene_id)
)

n_missing_ensembl <- sum(
  is.na(prepared$gene_id)
)

n_duplicated_ensembl <- sum(
  duplicated(
    prepared$gene_id[
      !is.na(prepared$gene_id)
    ]
  )
)

# ------------------------------------------------------------
# 17b. Auditoria do identificador original
# ------------------------------------------------------------

n_duplicated_original <- sum(
  duplicated(prepared$original_id)
)

if (n_duplicated_original > 0) {
  warning(
    "Foram encontrados ",
    n_duplicated_original,
    " original_id duplicados."
  )
}

# ------------------------------------------------------------
# 18. Resumo
# ------------------------------------------------------------
cat(
  "\n",
  "============================================\n",
  "PREPARAÇÃO DA EXPRESSÃO — GSE53697\n",
  "============================================\n",
  "Genes no arquivo original: ",
  n_genes_initial,
  "\n",
  "Genes all-zero removidos: ",
  n_all_zero,
  "\n",
  "Genes após filtro: ",
  n_genes_filtered,
  "\n",
  "Original IDs duplicados: ",
  n_duplicated_original,
  "\n",
  "Amostras: ",
  ncol(log_expr),
  "\n",
  "Genes com Ensembl: ",
  n_ensembl,
  " / ",
  nrow(prepared),
  "\n",
  "Genes sem Ensembl: ",
  n_missing_ensembl,
  "\n",
  "Ensembl duplicados: ",
  n_duplicated_ensembl,
  "\n",
  "Transformação: log2(RPKM + 1)\n",
  "============================================\n",
  sep = ""
)


# ------------------------------------------------------------
# 19. Salvar
# ------------------------------------------------------------
write_csv(
  prepared,
  output_file,
  na = "NA"
)

message(
  "\nMatriz preparada salva em:\n",
  output_file
)