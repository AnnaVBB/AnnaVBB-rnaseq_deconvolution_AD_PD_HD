#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(stringr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Uso: Rscript 00_auditoria.R <counts.tsv.gz> <annot.tsv.gz> <soft.gz> <out_report.txt>")
}

counts_file <- args[1]
annot_file  <- args[2]
soft_file   <- args[3]
output_txt  <- args[4]

report_lines <- c()
add_to_report <- function(...) {
  report_lines <<- c(report_lines, paste(...))
}

add_to_report("==================================================================")
add_to_report("           RELATÓRIO DE AUDITORIA DE INTEGRIDADE DE DADOS         ")
add_to_report("==================================================================")
add_to_report("Data: ", as.character(Sys.time()))
add_to_report("Counts: ", counts_file)
add_to_report("Annot:  ", annot_file)
add_to_report("SOFT:   ", soft_file)
add_to_report("------------------------------------------------------------------\n")

# 1. Auditoria de Counts
add_to_report(">>> 1. AUDITORIA DA MATRIZ DE COUNTS")
counts_df <- fread(cmd = paste("zcat", counts_file), data.table = FALSE)
sample_cols <- colnames(counts_df)[-1]
gene_ids <- counts_df[, 1]

add_to_report(" - Nome da coluna de ID: ", colnames(counts_df)[1])
add_to_report(" - Total de genes (linhas): ", nrow(counts_df))
add_to_report(" - Total de amostras (colunas): ", length(sample_cols))
add_to_report(" - GeneIDs duplicados: ", sum(duplicated(gene_ids)))

counts_mat <- as.matrix(counts_df[, -1])
add_to_report(" - Valores ausentes (NA): ", sum(is.na(counts_mat)))
add_to_report(" - Contém valores negativos: ", any(counts_mat < 0, na.rm = TRUE))

lib_sizes <- colSums(counts_mat, na.rm = TRUE)
add_to_report(" - Mediana do tamanho de biblioteca: ", format(median(lib_sizes), big.mark="."))
add_to_report("\n------------------------------------------------------------------\n")

# 2. Auditoria de Anotação
add_to_report(">>> 2. AUDITORIA DA ANOTAÇÃO GÊNICA")
annot_df <- fread(cmd = paste("zcat", annot_file), data.table = FALSE)
annot_id_col <- colnames(annot_df)[1]

add_to_report(" - Total de genes anotados: ", nrow(annot_df))
add_to_report(" - GeneIDs duplicados na anotação: ", sum(duplicated(annot_df[[annot_id_col]])))

common_genes <- intersect(as.character(gene_ids), as.character(annot_df[[annot_id_col]]))
pct_mapped <- (length(common_genes) / length(gene_ids)) * 100
add_to_report(" - Mapeamento com a matriz: ", length(common_genes), sprintf(" (%.2f%%)", pct_mapped))
add_to_report("\n------------------------------------------------------------------\n")

# 3. Auditoria do SOFT
add_to_report(">>> 3. AUDITORIA DOS METADADOS (SOFT FILE)")
gds <- getGEO(filename = soft_file)
gsm_list <- GSMList(gds)
soft_gsms <- names(gsm_list)

add_to_report(" - Amostras no SOFT: ", length(soft_gsms))
add_to_report(" - GSMs pareadas (SOFT e Counts): ", length(intersect(sample_cols, soft_gsms)))

missing_in_counts <- setdiff(soft_gsms, sample_cols)
add_to_report(" - GSMs no SOFT ausentes nos Counts: ", length(missing_in_counts))
if (length(missing_in_counts) > 0) {
  add_to_report("   [Ausentes]: ", paste(missing_in_counts, collapse = ", "))
}

writeLines(report_lines, con = output_txt)
cat("[SUCCESS] Relatório de auditoria salvo em: ", output_txt, "\n")
