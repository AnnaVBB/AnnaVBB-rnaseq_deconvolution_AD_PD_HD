#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(readr)
  library(dplyr)
  library(tibble)
})

# ============================================================
# 02c_pathway_enrichment.R
#
# Realiza análise de enriquecimento funcional GO Biological
# Process (ORA - Over-Representation Analysis) a partir dos
# genes diferencialmente expressos identificados no script 02b.
#
# O conjunto de genes testados na análise diferencial é usado
# como universo/background do enriquecimento.
# ============================================================


# ------------------------------------------------------------
# 0. Receber argumentos do terminal
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Uso: Rscript 02c_pathway_enrichment.R ",
    "<deg_results.csv> <out_prefix>"
  )
}

deg_file   <- args[1]
out_prefix <- args[2]


# ------------------------------------------------------------
# 1. Carregar resultados da análise diferencial
# ------------------------------------------------------------

deg_res <- readr::read_csv(
  deg_file,
  show_col_types = FALSE
)

# Verificar se as colunas necessárias estão presentes.

required_cols <- c(
  "gene_symbol",
  "logFC",
  "adj.P.Val"
)

missing_cols <- setdiff(
  required_cols,
  colnames(deg_res)
)

if (length(missing_cols) > 0) {
  stop(
    "Colunas necessárias ausentes no arquivo DEG: ",
    paste(missing_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 2. Definir arquivo de saída
# ------------------------------------------------------------

go_out <- paste0(
  out_prefix,
  "_pathways_GO_BP.csv"
)

# Criar diretório de saída caso ainda não exista.

dir.create(
  dirname(go_out),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Selecionar DEGs significativos
#
# Critérios definidos no projeto:
# FDR < 0.05
# |log2FC| > 0.58
# ------------------------------------------------------------

sig_res <- deg_res %>%
  dplyr::filter(
    !is.na(adj.P.Val),
    !is.na(logFC),
    adj.P.Val < 0.05,
    abs(logFC) > 0.58
  )

cat(
  "[INFO] Total de DEGs significativos:",
  nrow(sig_res),
  "\n"
)


# ------------------------------------------------------------
# 4. Obter símbolos dos DEGs
# ------------------------------------------------------------

sig_genes <- sig_res %>%
  dplyr::filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    gene_symbol != "-"
  ) %>%
  dplyr::pull(gene_symbol) %>%
  unique()

cat(
  "[INFO] Símbolos disponíveis para enriquecimento:",
  length(sig_genes),
  "\n"
)


# ------------------------------------------------------------
# 5. Obter universo/background
#
# O arquivo produzido pelo 02b contém todos os genes que
# permaneceram após a filtragem por expressão e foram testados
# pelo limma-voom.
#
# Portanto, esses genes representam um background mais adequado
# para a análise de enriquecimento do que utilizar todos os genes
# presentes no banco org.Hs.eg.db.
# ------------------------------------------------------------

background_symbols <- deg_res %>%
  dplyr::filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    gene_symbol != "-"
  ) %>%
  dplyr::pull(gene_symbol) %>%
  unique()

cat(
  "[INFO] Símbolos disponíveis no background:",
  length(background_symbols),
  "\n"
)


# ------------------------------------------------------------
# 6. Caso existam poucos DEGs, gerar arquivo vazio
# ------------------------------------------------------------

# Uma ORA com pouquíssimos genes tende a ser pouco informativa.
# Mantemos o limite mínimo de 5 genes para prosseguir.

if (length(sig_genes) < 5) {

  warning(
    "Poucos genes significativos disponíveis para enriquecimento."
  )

  readr::write_csv(
    tibble::tibble(),
    go_out
  )

  quit(
    save = "no",
    status = 0
  )
}


# ------------------------------------------------------------
# 7. Converter DEGs de SYMBOL para ENTREZID
# ------------------------------------------------------------

gene_ids <- clusterProfiler::bitr(
  sig_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
) %>%
  dplyr::distinct(
    ENTREZID,
    .keep_all = TRUE
  )

cat(
  "[INFO] DEGs convertidos para ENTREZID:",
  nrow(gene_ids),
  "\n"
)


# ------------------------------------------------------------
# 8. Converter background de SYMBOL para ENTREZID
# ------------------------------------------------------------

background_ids <- clusterProfiler::bitr(
  background_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
) %>%
  dplyr::distinct(
    ENTREZID,
    .keep_all = TRUE
  )

cat(
  "[INFO] Genes do background convertidos para ENTREZID:",
  nrow(background_ids),
  "\n"
)


# ------------------------------------------------------------
# 9. Verificar se há genes suficientes após conversão
# ------------------------------------------------------------

if (nrow(gene_ids) < 5) {

  warning(
    "Poucos DEGs puderam ser convertidos para ENTREZID."
  )

  readr::write_csv(
    tibble::tibble(),
    go_out
  )

  quit(
    save = "no",
    status = 0
  )
}


# ------------------------------------------------------------
# 10. GO Biological Process
#
# gene:
#   DEGs significativos.
#
# universe:
#   todos os genes testados na análise diferencial e que
#   possuem anotação válida.
#
# ont = "BP":
#   utiliza a ontologia Biological Process.
#
# pAdjustMethod = "BH":
#   aplica correção de Benjamini-Hochberg.
# ------------------------------------------------------------

ego <- clusterProfiler::enrichGO(
  gene = gene_ids$ENTREZID,
  universe = background_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# ------------------------------------------------------------
# 11. Converter resultado para data.frame
# ------------------------------------------------------------

go_res <- as.data.frame(
  ego
)

cat(
  "[INFO] Termos GO BP significativos encontrados:",
  nrow(go_res),
  "\n"
)


# ------------------------------------------------------------
# 12. Salvar resultado
# ------------------------------------------------------------

# Mesmo quando nenhum termo significativo é encontrado,
# o arquivo é criado. Isso é importante para manter o pipeline
# reprodutível e permitir que o Snakemake reconheça a execução
# da etapa.

readr::write_csv(
  go_res,
  go_out
)

cat(
  "[SUCCESS] GO Biological Process salvo em:",
  go_out,
  "\n"
)