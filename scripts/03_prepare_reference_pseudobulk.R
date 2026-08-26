#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Uso: Rscript 03_prepare_reference_pseudobulk.R <sc_matrix.rds/tsv> <sc_meta.csv> <out_ref_rds>")
}

sc_matrix_file <- args[1]
sc_meta_file   <- args[2]
output_rds     <- args[3]

cat("[INFO] Carregando matriz e metadados de single-cell (Mathys 2019)...\n")

# Se for RDS ou TSV
if (endsWith(sc_matrix_file, ".rds")) {
  sc_counts <- readRDS(sc_matrix_file)
} else {
  sc_counts <- fread(sc_matrix_file, data.table = FALSE)
  rownames(sc_counts) <- sc_counts[, 1]
  sc_counts <- as.matrix(sc_counts[, -1])
}

sc_meta <- read_csv(sc_meta_file, show_col_types = FALSE)

# Garantir tipos celulares padrão no Córtex
valid_cell_types <- c("Ast", "End", "Ex", "In", "Mic", "Oli", "OPC")
sc_meta_clean <- sc_meta %>% 
  filter(broad.cell.type %in% valid_cell_types)

common_cells <- intersect(colnames(sc_counts), sc_meta_clean$cell_id)
sc_counts_sub <- sc_counts[, common_cells]

cat("[INFO] Matriz Single-Cell filtrada: ", ncol(sc_counts_sub), " células e ", nrow(sc_counts_sub), " genes.\n")

# Construir perfil médio por tipo celular (Signature Matrix)
cell_types <- sc_meta_clean$broad.cell.type[match(common_cells, sc_meta_clean$cell_id)]
unique_types <- unique(cell_types)

sig_matrix <- matrix(0, nrow = nrow(sc_counts_sub), ncol = length(unique_types))
rownames(sig_matrix) <- rownames(sc_counts_sub)
colnames(sig_matrix) <- unique_types

for (ct in unique_types) {
  cells_in_type <- which(cell_types == ct)
  if (length(cells_in_type) > 1) {
    sig_matrix[, ct] <- rowMeans(sc_counts_sub[, cells_in_type], na.rm = TRUE)
  } else {
    sig_matrix[, ct] <- sc_counts_sub[, cells_in_type]
  }
}

ref_data <- list(
  signature_matrix = sig_matrix,
  raw_sc_counts    = sc_counts_sub,
  sc_metadata      = sc_meta_clean
)

saveRDS(ref_data, file = output_rds)
cat("[SUCCESS] Matriz de Referência Single-Cell salva em: ", output_rds, "\n")
