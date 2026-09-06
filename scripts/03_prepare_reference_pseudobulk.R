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

sc_matrix_file <- args[1] #guarda o caminho de expressão single-cell/single-nucleus
sc_meta_file   <- args[2] 
output_rds     <- args[3]

cat("[INFO] Carregando matriz e metadados de single-cell (Mathys 2019)...\n")

# Se for RDS ou TSV
if (endsWith(sc_matrix_file, ".rds")) {
  sc_counts <- readRDS(sc_matrix_file)
} else {
  sc_counts <- fread(sc_matrix_file, data.table = FALSE)
  rownames(sc_counts) <- sc_counts[, 1]
  sc_counts <- as.matrix(sc_counts[, -1]) #todas as colunas, exceto a primeira. Os genes estão nas linhas e as células nas colunas
}

sc_meta <- read_csv(sc_meta_file, show_col_types = FALSE)

# Garantir tipos celulares padrão no Córtex
#Astrocytes, Endothelial cells, Excitatory neurons, inhibitory neurons, microglia, oligodendrocytes, oligodendrocyte precursor cells
valid_cell_types <- c("Ast", "End", "Ex", "In", "Mic", "Oli", "OPC")
sc_meta_clean <- sc_meta %>% 
  filter(broad.cell.type %in% valid_cell_types)

common_cells <- intersect(colnames(sc_counts), sc_meta_clean$cell_id)
sc_counts_sub <- sc_counts[, common_cells] #selecionamos da matriz apenas as colunas presentes em common_cells

cat("[INFO] Matriz Single-Cell filtrada: ", ncol(sc_counts_sub), " células e ", nrow(sc_counts_sub), " genes.\n")

# Construir perfil médio por tipo celular (Signature Matrix)
cell_types <- sc_meta_clean$broad.cell.type[match(common_cells, sc_meta_clean$cell_id)] #"Para cada célula de common_cells, em qual posição ela aparece dentro de sc_meta_clean$cell_id?"
unique_types <- unique(cell_types) #remove os valores repetidos

sig_matrix <- matrix(0, nrow = nrow(sc_counts_sub), ncol = length(unique_types)) #Matriz inicialmente preenchida por zeros, o número de linhas seá o número de genes e o número de colunas o número de tipos celulares
rownames(sig_matrix) <- rownames(sc_counts_sub) #os genes da matriz original passam a ser os genes da signature matrix
colnames(sig_matrix) <- unique_types #As colunas passam a ter o nome dos tipos celulares

#Passamos por cada tipo celular, encontrando as células daquele tipo
for (ct in unique_types) {
  cells_in_type <- which(cell_types == ct)
  #Se houver várias células daquele tipo entramos nesse bloco
  if (length(cells_in_type) > 1) {
    #Calculamos a expressão média de cada gene entre todas as células daquele tipo
    sig_matrix[, ct] <- rowMeans(sc_counts_sub[, cells_in_type], na.rm = TRUE)
  } else {
    #No caso de não existir uma média a ser calculada, usamos a expressão daquela única célula
    sig_matrix[, ct] <- sc_counts_sub[, cells_in_type]
  }
}

#guardamos em uma lista a signature_matrix, as contagens originais das células que foram utilizadas e os metadados das células selecionadas
ref_data <- list(
  signature_matrix = sig_matrix,
  raw_sc_counts    = sc_counts_sub,
  sc_metadata      = sc_meta_clean
)

#Salvar toda a referência em um .rds
saveRDS(ref_data, file = output_rds)
cat("[SUCCESS] Matriz de Referência Single-Cell salva em: ", output_rds, "\n")
