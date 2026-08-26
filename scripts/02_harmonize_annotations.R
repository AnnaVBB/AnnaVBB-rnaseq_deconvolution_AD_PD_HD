#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Uso: Rscript 02_harmonize_annotations.R <counts.tsv.gz> <annot.tsv.gz> <out_rds>")
}

counts_file <- args[1]
annot_file  <- args[2]
output_rds  <- args[3]

# 1. Carregar contagens e anotação
counts_df <- fread(cmd = paste("zcat", counts_file), data.table = FALSE)
annot_df  <- fread(cmd = paste("zcat", annot_file), data.table = FALSE)

id_col <- colnames(counts_df)[1]
colnames(counts_df)[1] <- "GeneID"
colnames(annot_df)[1]  <- "GeneID"

# 2. Mesclar matriz com anotação HGNC/Ensembl
counts_annotated <- merge(annot_df, counts_df, by = "GeneID")

# 3. Filtrar genes sem Symbol ou Ensembl válido
counts_clean <- counts_annotated %>%
  filter(!is.na(Symbol) & Symbol != "" & Symbol != "-") %>%
  filter(!is.na(EnsemblGeneID) & EnsemblGeneID != "" & EnsemblGeneID != "-")

# Tratar duplicidades de Symbol mantendo o gene com maior contagem total
sample_cols <- setdiff(colnames(counts_df), "GeneID")
mat_only <- as.matrix(counts_clean[, sample_cols])
rowSums_val <- rowSums(mat_only)

counts_clean$row_sum <- rowSums_val
counts_dedup <- counts_clean %>%
  arrange(desc(row_sum)) %>%
  distinct(Symbol, .keep_all = TRUE) %>%
  select(-row_sum)

# 4. Filtragem por baixa expressão (cpm > 0.5 em pelo menos 20% das amostras)
mat_final <- as.matrix(counts_dedup[, sample_cols])
rownames(mat_final) <- counts_dedup$Symbol

keep <- rowSums(cpm(mat_final) > 0.5) >= (0.20 * ncol(mat_final))
mat_filtered <- mat_final[keep, ]
gene_info    <- counts_dedup[keep, c("GeneID", "Symbol", "EnsemblGeneID", "GeneType")]

# 5. Salvar objeto R serializado contendo matriz e mapa de genes
res_list <- list(
  counts = mat_filtered,
  annotation = gene_info
)

saveRDS(res_list, file = output_rds)
cat("[SUCCESS] Matriz harmonizada salva em: ", output_rds, "\n")
cat("[INFO] Genes retidos pós-filtro: ", nrow(mat_filtered), "\n")
