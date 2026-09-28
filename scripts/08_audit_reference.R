#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

# ============================================================
# 06 — AUDITORIA DA REFERÊNCIA snRNA-seq MATHYS 2019
# ============================================================
#
# Objetivos:
#   1. Confirmar a integridade estrutural da matriz filtrada
#   2. Verificar correspondência genes × matriz × núcleos
#   3. Inspecionar as colunas reais do metadata
#   4. Inspecionar os arquivos auxiliares de metadata
#
# IMPORTANTE:
#   Este script NÃO normaliza, filtra ou modifica a referência.
# ============================================================

ref_dir <- "data/reference/mathys2019/raw"

matrix_file <- file.path(ref_dir, "filtered_count_matrix.mtx")
genes_file  <- file.path(ref_dir, "filtered_gene_row_names.txt")
cells_file  <- file.path(ref_dir, "filtered_column_metadata.txt")

metadata_files <- c(
  sample_key  = "snRNAseqPFC_BA10_Sample_key.csv",
  assay       = "snRNAseqPFC_BA10_assay_scRNAseq_metadata.csv",
  biospecimen = "snRNAseqPFC_BA10_biospecimen_metadata.csv",
  id_mapping  = "snRNAseqPFC_BA10_id_mapping.csv"
)

# ============================================================
# 1. VERIFICAÇÃO DOS ARQUIVOS
# ============================================================

cat("\n============================================\n")
cat("1. ARQUIVOS DA REFERÊNCIA\n")
cat("============================================\n\n")

all_files <- c(
  matrix_file,
  genes_file,
  cells_file,
  file.path(ref_dir, metadata_files)
)

missing_files <- all_files[!file.exists(all_files)]

if (length(missing_files) > 0) {
  cat("Arquivos ausentes:\n")
  print(missing_files)
  stop("Auditoria interrompida: existem arquivos ausentes.")
}

cat("Todos os arquivos esperados foram encontrados.\n\n")

file_info <- data.frame(
  file = basename(all_files),
  size_MB = round(file.info(all_files)$size / 1024^2, 2),
  stringsAsFactors = FALSE
)

print(file_info, row.names = FALSE)

# ============================================================
# 2. MATRIZ DE CONTAGENS
# ============================================================

cat("\n============================================\n")
cat("2. MATRIZ DE CONTAGENS\n")
cat("============================================\n\n")

counts <- Matrix::readMM(matrix_file)

cat("Classe:", class(counts)[1], "\n")
cat("Número de linhas :", nrow(counts), "\n")
cat("Número de colunas:", ncol(counts), "\n")
cat("Valores não-zero:", length(counts@x), "\n")

total_entries <- as.double(nrow(counts)) * as.double(ncol(counts))
sparsity <- 1 - (length(counts@x) / total_entries)

cat(
  "Proporção de zeros:",
  round(sparsity * 100, 2),
  "%\n"
)

cat(
  "Menor valor não-zero:",
  min(counts@x),
  "\n"
)

cat(
  "Maior valor:",
  max(counts@x),
  "\n"
)

# ============================================================
# 3. ARQUIVO DE GENES
# ============================================================

cat("\n============================================\n")
cat("3. GENES\n")
cat("============================================\n\n")

genes <- data.table::fread(
  genes_file,
  header = FALSE,
  data.table = FALSE
)

cat("Número de linhas :", nrow(genes), "\n")
cat("Número de colunas:", ncol(genes), "\n")

cat("\nPrimeiras linhas:\n")
print(head(genes))

# ============================================================
# 4. METADATA DOS NÚCLEOS
# ============================================================

cat("\n============================================\n")
cat("4. METADATA DOS NÚCLEOS\n")
cat("============================================\n\n")

cell_meta <- data.table::fread(
  cells_file,
  data.table = FALSE
)

cat("Número de linhas :", nrow(cell_meta), "\n")
cat("Número de colunas:", ncol(cell_meta), "\n")

cat("\nNomes das colunas:\n")
print(colnames(cell_meta))

cat("\nPrimeiras linhas:\n")
print(head(cell_meta))

# ============================================================
# 5. CONSISTÊNCIA ESTRUTURAL
# ============================================================

cat("\n============================================\n")
cat("5. CONSISTÊNCIA ESTRUTURAL\n")
cat("============================================\n\n")

genes_match <- nrow(counts) == nrow(genes)
cells_match <- ncol(counts) == nrow(cell_meta)

cat(
  "Linhas da matriz == genes:",
  genes_match,
  "\n"
)

cat(
  "Colunas da matriz == linhas do metadata:",
  cells_match,
  "\n"
)

if (!genes_match) {
  warning("Número de genes não corresponde ao número de linhas da matriz.")
}

if (!cells_match) {
  warning("Número de núcleos não corresponde ao número de colunas da matriz.")
}

# ============================================================
# 6. AUDITORIA DOS METADADOS AUXILIARES
# ============================================================

cat("\n============================================\n")
cat("6. METADADOS AUXILIARES\n")
cat("============================================\n")

aux_metadata <- list()

for (nm in names(metadata_files)) {

  f <- file.path(ref_dir, metadata_files[[nm]])

  cat("\n--------------------------------------------\n")
  cat("Arquivo:", metadata_files[[nm]], "\n")
  cat("Objeto :", nm, "\n")
  cat("--------------------------------------------\n")

  dat <- data.table::fread(
    f,
    data.table = FALSE
  )

  aux_metadata[[nm]] <- dat

  cat("Linhas :", nrow(dat), "\n")
  cat("Colunas:", ncol(dat), "\n")

  cat("\nNomes das colunas:\n")
  print(colnames(dat))

  cat("\nPrimeiras 3 linhas:\n")
  print(head(dat, 3))
}

# ============================================================
# 7. RESUMO
# ============================================================

cat("\n============================================\n")
cat("7. RESUMO DA AUDITORIA INICIAL\n")
cat("============================================\n\n")

cat("Genes na matriz :", nrow(counts), "\n")
cat("Núcleos         :", ncol(counts), "\n")
cat(
  "Sparsidade     :",
  round(sparsity * 100, 2),
  "% zeros\n"
)

cat(
  "Genes compatíveis com matriz:",
  genes_match,
  "\n"
)

cat(
  "Metadata compatível com matriz:",
  cells_match,
  "\n"
)

cat("\nArquivos auxiliares carregados:\n")

for (nm in names(aux_metadata)) {
  cat(
    " -", nm, ":",
    nrow(aux_metadata[[nm]]), "x",
    ncol(aux_metadata[[nm]]), "\n"
  )
}

cat("\nAuditoria estrutural inicial concluída.\n")
cat("Nenhuma transformação foi aplicada aos dados.\n")

# ============================================================
# 8. AUDITORIA BIOLÓGICA DA REFERÊNCIA
# ============================================================

cat("\n============================================\n")
cat("8. AUDITORIA BIOLÓGICA DA REFERÊNCIA\n")
cat("============================================\n\n")

# ------------------------------------------------------------
# 8.1 Tipos celulares amplos
# ------------------------------------------------------------

cat("=== TIPOS CELULARES (broad.cell.type) ===\n\n")

broad_counts <- sort(
  table(cell_meta$broad.cell.type, useNA = "ifany"),
  decreasing = TRUE
)

print(broad_counts)

cat(
  "\nNúmero de tipos celulares:",
  length(unique(cell_meta$broad.cell.type)),
  "\n"
)

# ------------------------------------------------------------
# 8.2 Subclusters
# ------------------------------------------------------------

cat("\n=== SUBCLUSTERS ===\n\n")

subcluster_counts <- sort(
  table(cell_meta$Subcluster, useNA = "ifany"),
  decreasing = TRUE
)

print(subcluster_counts)

cat(
  "\nNúmero de subclusters:",
  length(unique(cell_meta$Subcluster)),
  "\n"
)

# ------------------------------------------------------------
# 8.3 Indivíduos / projid
# ------------------------------------------------------------

cat("\n=== INDIVÍDUOS REPRESENTADOS ===\n\n")

n_projids <- length(unique(cell_meta$projid))

cat(
  "Número de projid na matriz:",
  n_projids,
  "\n"
)

nuclei_per_projid <- sort(
  table(cell_meta$projid),
  decreasing = TRUE
)

cat("\nNúcleos por projid:\n")
print(nuclei_per_projid)

# ------------------------------------------------------------
# 8.4 Missing values
# ------------------------------------------------------------

cat("\n=== MISSING VALUES ===\n\n")

important_columns <- c(
  "TAG",
  "projid",
  "broad.cell.type",
  "Subcluster"
)

for (col in important_columns) {
  cat(
    col, ":",
    sum(is.na(cell_meta[[col]])),
    "NA\n"
  )
}

# ------------------------------------------------------------
# 8.5 Genes duplicados
# ------------------------------------------------------------

cat("\n=== GENES DUPLICADOS ===\n\n")

gene_names <- genes[[1]]

cat(
  "Genes totais:",
  length(gene_names),
  "\n"
)

cat(
  "Genes únicos:",
  length(unique(gene_names)),
  "\n"
)

cat(
  "Genes duplicados:",
  sum(duplicated(gene_names)),
  "\n"
)

if (anyDuplicated(gene_names) > 0) {

  duplicated_genes <- unique(
    gene_names[duplicated(gene_names)]
  )

  cat("\nPrimeiros genes duplicados:\n")
  print(head(duplicated_genes, 20))
}

# ------------------------------------------------------------
# 8.6 Compatibilidade dos projid
# ------------------------------------------------------------

cat("\n=== COMPATIBILIDADE DOS PROJID ===\n\n")

biospecimen <- aux_metadata$biospecimen
sample_key  <- aux_metadata$sample_key

projid_cells <- unique(as.character(cell_meta$projid))
projid_bio   <- unique(as.character(biospecimen$projid))
projid_key   <- unique(as.character(sample_key$projid))

cat(
  "projid únicos - cell metadata:",
  length(projid_cells),
  "\n"
)

cat(
  "projid únicos - biospecimen:",
  length(projid_bio),
  "\n"
)

cat(
  "projid únicos - sample key:",
  length(projid_key),
  "\n"
)

cat("\nprojid da matriz ausentes no biospecimen:\n")
print(setdiff(projid_cells, projid_bio))

cat("\nprojid da matriz ausentes no sample key:\n")
print(setdiff(projid_cells, projid_key))

# ------------------------------------------------------------
# 8.7 Relação projid -> individualID
# ------------------------------------------------------------

cat("\n=== PROJID -> INDIVIDUAL ID ===\n\n")

donor_map <- unique(
  biospecimen[, c("projid", "individualID")]
)

cat(
  "Número de pares projid/individualID:",
  nrow(donor_map),
  "\n"
)

cat(
  "Número de individualID únicos:",
  length(unique(donor_map$individualID)),
  "\n"
)

duplicated_projid <- donor_map$projid[
  duplicated(donor_map$projid)
]

if (length(duplicated_projid) == 0) {
  cat("Cada projid possui um único individualID: TRUE\n")
} else {
  cat("Cada projid possui um único individualID: FALSE\n")
  print(unique(duplicated_projid))
}

# ------------------------------------------------------------
# 8.8 Adicionar individualID ao metadata em memória
# ------------------------------------------------------------

cell_meta_audit <- merge(
  cell_meta,
  donor_map,
  by = "projid",
  all.x = TRUE,
  sort = FALSE
)

cat(
  "\nNúcleos sem individualID após merge:",
  sum(is.na(cell_meta_audit$individualID)),
  "\n"
)

cat(
  "Número de individualID representados:",
  length(unique(cell_meta_audit$individualID)),
  "\n"
)

# ------------------------------------------------------------
# 8.9 Tabela donor × cell type
# ------------------------------------------------------------

cat("\n=== NÚCLEOS POR INDIVÍDUO × TIPO CELULAR ===\n\n")

donor_cell_table <- table(
  cell_meta_audit$individualID,
  cell_meta_audit$broad.cell.type
)

print(donor_cell_table)

cat("\nNúmero de combinações donor × cell type com zero núcleos:\n")

cat(
  sum(donor_cell_table == 0),
  "de",
  length(donor_cell_table),
  "\n"
)

cat("\nAuditoria biológica concluída.\n")
