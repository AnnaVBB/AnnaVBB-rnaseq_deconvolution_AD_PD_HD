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

# 1. Carregar contagens e anotação (leitura direta sem dependência de zcat)
counts_df <- fread(counts_file, data.table = FALSE)
annot_df  <- fread(annot_file, data.table = FALSE)

# Normalizar o nome da primeira coluna para GeneID
colnames(counts_df)[1] <- "GeneID"
colnames(annot_df)[1]  <- "GeneID"

# 2. Converter IDs para character para evitar falhas de tipo no merge
counts_df$GeneID <- as.character(counts_df$GeneID)
annot_df$GeneID  <- as.character(annot_df$GeneID)

# Mesclar matriz com anotação
counts_annotated <- merge(annot_df, counts_df, by = "GeneID")

# Identificar a coluna com os símbolos dos genes (Symbol ou GeneSymbol)
symbol_col <- intersect(c("Symbol", "GeneSymbol", "symbol"), colnames(annot_df))[1]
if (is.na(symbol_col)) {
  stop("Não foi possível encontrar a coluna de Símbolo do Gene no arquivo de anotação.")
}

# 3. Limpeza de genes inválidos
counts_clean <- counts_annotated %>%
  filter(!is.na(.data[[symbol_col]]) & .data[[symbol_col]] != "" & .data[[symbol_col]] != "-")

# Resolver duplicidades mantendo a linha de maior contagem total
sample_cols <- setdiff(colnames(counts_df), "GeneID")
mat_only <- as.matrix(counts_clean[, sample_cols])
rowSums_val <- rowSums(mat_only, na.rm = TRUE)

counts_clean$row_sum <- rowSums_val
counts_dedup <- counts_clean %>%
  arrange(desc(row_sum)) %>%
  distinct(.data[[symbol_col]], .keep_all = TRUE) %>%
  select(-row_sum)

# 4. Filtragem por baixa expressão (cpm > 0.5 em pelo menos 20% das amostras)
mat_final <- as.matrix(counts_dedup[, sample_cols])
rownames(mat_final) <- counts_dedup[[symbol_col]]

# Garantir dados numéricos
class(mat_final) <- "numeric"

keep <- rowSums(cpm(mat_final) > 0.5, na.rm = TRUE) >= (0.20 * ncol(mat_final))
mat_filtered <- mat_final[keep, ]

# Salvar lista final
res_list <- list(
  counts = mat_filtered,
  annotation = counts_dedup[keep, setdiff(colnames(annot_df), sample_cols)]
)

saveRDS(res_list, file = output_rds)
cat("[SUCCESS] Matriz harmonizada salva em: ", output_rds, "\n")
cat("[INFO] Genes retidos pós-filtro: ", nrow(mat_filtered), "\n")