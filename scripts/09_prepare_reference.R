#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

# ============================================================
# 07_prepare_reference.R
#
# Prepara a referência snRNA-seq Mathys et al. 2019
# para deconvolução celular.
#
# IMPORTANTE:
# - dados Mathys/ROSMAP são de acesso controlado;
# - outputs derivados permanecem locais;
# - nenhuma normalização é aplicada à matriz original;
# - identidade dos doadores é preservada;
# - matriz permanece esparsa.
# ============================================================


# ============================================================
# 1. CAMINHOS
# ============================================================

ref_dir <- "data/reference/mathys2019/raw"
out_dir <- "data/reference/mathys2019/prepared"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(
  ref_dir,
  "filtered_count_matrix.mtx"
)

gene_file <- file.path(
  ref_dir,
  "filtered_gene_row_names.txt"
)

cell_file <- file.path(
  ref_dir,
  "filtered_column_metadata.txt"
)

biospecimen_file <- file.path(
  ref_dir,
  "snRNAseqPFC_BA10_biospecimen_metadata.csv"
)


# ============================================================
# 2. VERIFICAR ARQUIVOS
# ============================================================

cat("\n============================================\n")
cat("1. VERIFICAÇÃO DOS ARQUIVOS\n")
cat("============================================\n\n")

required_files <- c(
  matrix_file,
  gene_file,
  cell_file,
  biospecimen_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {

  cat("Arquivos ausentes:\n")
  print(missing_files)

  stop(
    "Preparação interrompida: existem arquivos obrigatórios ausentes."
  )
}

cat("Todos os arquivos obrigatórios foram encontrados.\n")


# ============================================================
# 3. CARREGAR MATRIZ E METADADOS
# ============================================================

cat("\n============================================\n")
cat("2. CARREGAMENTO DOS DADOS\n")
cat("============================================\n\n")

cat("Lendo matriz de contagens...\n")

counts <- Matrix::readMM(matrix_file)

# Converter dgTMatrix -> dgCMatrix.
# O formato dgCMatrix é mais conveniente para operações posteriores.
counts <- as(counts, "CsparseMatrix")

cat(
  "Matriz:",
  nrow(counts), "genes x",
  ncol(counts), "núcleos\n"
)

cat("Lendo genes...\n")

genes <- data.table::fread(
  gene_file,
  header = FALSE
)

setnames(genes, "V1", "gene_symbol")

cat("Genes:", nrow(genes), "\n")

cat("Lendo metadata dos núcleos...\n")

cell_meta <- data.table::fread(
  cell_file
)

cat("Núcleos:", nrow(cell_meta), "\n")

cat("Lendo biospecimen metadata...\n")

biospecimen <- data.table::fread(
  biospecimen_file
)

cat("Biospecimens:", nrow(biospecimen), "\n")


# ============================================================
# 4. VALIDAÇÃO ESTRUTURAL
# ============================================================

cat("\n============================================\n")
cat("3. VALIDAÇÃO ESTRUTURAL\n")
cat("============================================\n\n")

if (nrow(counts) != nrow(genes)) {
  stop("Número de genes incompatível com as linhas da matriz.")
}

if (ncol(counts) != nrow(cell_meta)) {
  stop("Número de núcleos incompatível com as colunas da matriz.")
}

if (anyDuplicated(genes$gene_symbol) > 0) {
  stop("Foram encontrados genes duplicados.")
}

required_cell_columns <- c(
  "TAG",
  "projid",
  "broad.cell.type",
  "Subcluster"
)

missing_cell_columns <- setdiff(
  required_cell_columns,
  names(cell_meta)
)

if (length(missing_cell_columns) > 0) {

  stop(
    paste(
      "Colunas ausentes no metadata celular:",
      paste(missing_cell_columns, collapse = ", ")
    )
  )
}

required_bio_columns <- c(
  "projid",
  "individualID"
)

missing_bio_columns <- setdiff(
  required_bio_columns,
  names(biospecimen)
)

if (length(missing_bio_columns) > 0) {

  stop(
    paste(
      "Colunas ausentes no biospecimen:",
      paste(missing_bio_columns, collapse = ", ")
    )
  )
}

cat("Estrutura matriz/genes/células: OK\n")


# ============================================================
# 5. MAPEAR PROJID -> INDIVIDUALID
# ============================================================

cat("\n============================================\n")
cat("4. MAPEAMENTO DOS DOADORES\n")
cat("============================================\n\n")

donor_map <- unique(
  biospecimen[, .(
    projid,
    individualID
  )]
)

if (anyDuplicated(donor_map$projid) > 0) {
  stop("projid não possui relação 1:1 com individualID.")
}

# match() preserva rigorosamente a ordem das células.
# Isso é fundamental porque a linha i de cell_meta deve continuar
# correspondendo à coluna i da matriz counts.

idx <- match(
  as.character(cell_meta$projid),
  as.character(donor_map$projid)
)

cell_meta[, individualID := donor_map$individualID[idx]]

if (anyNA(cell_meta$individualID)) {
  stop("Existem núcleos sem individualID após o mapeamento.")
}

cat(
  "Doadores representados:",
  uniqueN(cell_meta$individualID),
  "\n"
)


# ============================================================
# 6. HARMONIZAR TIPOS CELULARES
# ============================================================

cat("\n============================================\n")
cat("5. HARMONIZAÇÃO DOS TIPOS CELULARES\n")
cat("============================================\n\n")

celltype_map <- c(
  "Ex"  = "Excitatory",
  "In"  = "Inhibitory",
  "Ast" = "Astrocyte",
  "Oli" = "Oligodendrocyte",
  "Opc" = "OPC",
  "Mic" = "Microglia",
  "End" = "Endothelial",
  "Per" = "Pericyte"
)

observed_types <- sort(
  unique(cell_meta$broad.cell.type)
)

unknown_types <- setdiff(
  observed_types,
  names(celltype_map)
)

if (length(unknown_types) > 0) {

  stop(
    paste(
      "Tipos celulares não reconhecidos:",
      paste(unknown_types, collapse = ", ")
    )
  )
}

cell_meta[, cell_type := unname(
  celltype_map[broad.cell.type]
)]

cat("Mapeamento utilizado:\n")

print(
  data.table(
    original = names(celltype_map),
    harmonized = unname(celltype_map)
  )
)

cat("\nContagem de núcleos:\n")

print(
  sort(
    table(cell_meta$cell_type),
    decreasing = TRUE
  )
)


# ============================================================
# 7. IDENTIFICADORES INTERNOS ANÔNIMOS
# ============================================================

cat("\n============================================\n")
cat("6. IDENTIFICADORES INTERNOS\n")
cat("============================================\n\n")

# Criamos IDs internos para outputs derivados.
# individualID original continua apenas no ambiente/dados locais.

donors <- sort(
  unique(cell_meta$individualID)
)

anonymous_map <- data.table(
  individualID = donors,
  donor_id = sprintf(
    "Donor_%02d",
    seq_along(donors)
  )
)

idx_donor <- match(
  cell_meta$individualID,
  anonymous_map$individualID
)

cell_meta[, donor_id := anonymous_map$donor_id[idx_donor]]

cat(
  "IDs internos criados:",
  uniqueN(cell_meta$donor_id),
  "\n"
)


# ============================================================
# 8. QC DONOR × CELL TYPE
# ============================================================

cat("\n============================================\n")
cat("7. QC DONOR x CELL TYPE\n")
cat("============================================\n\n")

donor_cell_qc <- cell_meta[
  ,
  .(
    n_nuclei = .N
  ),
  by = .(
    donor_id,
    cell_type
  )
]

all_combinations <- CJ(
  donor_id = sort(unique(cell_meta$donor_id)),
  cell_type = sort(unique(cell_meta$cell_type)),
  unique = TRUE
)

donor_cell_qc <- merge(
  all_combinations,
  donor_cell_qc,
  by = c("donor_id", "cell_type"),
  all.x = TRUE,
  sort = TRUE
)

donor_cell_qc[
  is.na(n_nuclei),
  n_nuclei := 0L
]

cat(
  "Combinações donor x cell type:",
  nrow(donor_cell_qc),
  "\n"
)

cat(
  "Combinações com zero núcleos:",
  sum(donor_cell_qc$n_nuclei == 0),
  "\n"
)

cat("\nResumo por tipo celular:\n")

qc_summary <- donor_cell_qc[
  ,
  .(
    donors_total = .N,
    donors_present = sum(n_nuclei > 0),
    donors_absent = sum(n_nuclei == 0),
    min_nuclei = min(n_nuclei),
    median_nuclei = median(n_nuclei),
    max_nuclei = max(n_nuclei),
    total_nuclei = sum(n_nuclei)
  ),
  by = cell_type
][order(-total_nuclei)]

print(qc_summary)


# ============================================================
# 9. PSEUDOBULK DONOR × CELL TYPE
# ============================================================

cat("\n============================================\n")
cat("8. PSEUDOBULK DONOR x CELL TYPE\n")
cat("============================================\n\n")

# Cada coluna da matriz original corresponde a um núcleo.
# Agrupamos núcleos do mesmo donor e cell type por SOMA.
#
# A soma preserva a natureza de contagens e gera perfis
# donor-específicos por tipo celular.

group_id <- paste(
  cell_meta$donor_id,
  cell_meta$cell_type,
  sep = "__"
)

group_levels <- sort(unique(group_id))

group_factor <- factor(
  group_id,
  levels = group_levels
)

# Matriz esparsa:
# núcleos × grupos donor-celltype

aggregation_matrix <- sparse.model.matrix(
  ~ 0 + group_factor
)

# counts:
# genes × núcleos
#
# resultado:
# genes × donor-celltype

pseudobulk <- counts %*% aggregation_matrix

colnames(pseudobulk) <- group_levels
rownames(pseudobulk) <- genes$gene_symbol

cat(
  "Pseudobulk:",
  nrow(pseudobulk),
  "genes x",
  ncol(pseudobulk),
  "combinações observadas donor-celltype\n"
)

cat(
  "Soma total da matriz original:",
  sum(counts),
  "\n"
)

cat(
  "Soma total do pseudobulk:",
  sum(pseudobulk),
  "\n"
)

if (sum(counts) != sum(pseudobulk)) {
  stop("A soma das contagens mudou durante a agregação.")
}

cat("Conservação das contagens: OK\n")


# ============================================================
# 10. METADATA DO PSEUDOBULK
# ============================================================

cat("\n============================================\n")
cat("9. METADATA DO PSEUDOBULK\n")
cat("============================================\n\n")

pseudobulk_meta <- data.table(
  sample_id = group_levels
)

split_group <- tstrsplit(
  pseudobulk_meta$sample_id,
  "__",
  fixed = TRUE
)

pseudobulk_meta[, donor_id := split_group[[1]]]
pseudobulk_meta[, cell_type := split_group[[2]]]

nuclei_lookup <- donor_cell_qc[
  n_nuclei > 0,
  .(
    donor_id,
    cell_type,
    n_nuclei
  )
]

pseudobulk_meta <- merge(
  pseudobulk_meta,
  nuclei_lookup,
  by = c("donor_id", "cell_type"),
  all.x = TRUE,
  sort = FALSE
)

# Restaurar exatamente a ordem das colunas do pseudobulk.

pseudobulk_meta <- pseudobulk_meta[
  match(
    colnames(pseudobulk),
    sample_id
  )
]

if (!identical(
  pseudobulk_meta$sample_id,
  colnames(pseudobulk)
)) {
  stop("Ordem do pseudobulk metadata incompatível com a matriz.")
}

cat(
  "Perfis pseudobulk:",
  nrow(pseudobulk_meta),
  "\n"
)


# ============================================================
# 11. SALVAR OBJETOS LOCAIS
# ============================================================

cat("\n============================================\n")
cat("10. SALVANDO REFERÊNCIA PREPARADA\n")
cat("============================================\n\n")

# Matriz original esparsa.
saveRDS(
  counts,
  file.path(
    out_dir,
    "mathys2019_counts_sparse.rds"
  ),
  compress = FALSE
)

# Genes.
data.table::fwrite(
  genes,
  file.path(
    out_dir,
    "mathys2019_genes.tsv"
  ),
  sep = "\t"
)

# Metadata nucleus-level.
#
# ATENÇÃO:
# contém identificadores originais e permanece em diretório
# protegido pelo .gitignore.
data.table::fwrite(
  cell_meta,
  file.path(
    out_dir,
    "mathys2019_cell_metadata.tsv"
  ),
  sep = "\t"
)

# Mapeamento original -> ID interno.
#
# Também permanece privado/local.
data.table::fwrite(
  anonymous_map,
  file.path(
    out_dir,
    "mathys2019_donor_private_map.tsv"
  ),
  sep = "\t"
)

# QC com IDs internos.
data.table::fwrite(
  donor_cell_qc,
  file.path(
    out_dir,
    "mathys2019_donor_celltype_qc.tsv"
  ),
  sep = "\t"
)

data.table::fwrite(
  qc_summary,
  file.path(
    out_dir,
    "mathys2019_celltype_qc_summary.tsv"
  ),
  sep = "\t"
)

# Pseudobulk.
saveRDS(
  pseudobulk,
  file.path(
    out_dir,
    "mathys2019_pseudobulk_counts.rds"
  ),
  compress = FALSE
)

data.table::fwrite(
  pseudobulk_meta,
  file.path(
    out_dir,
    "mathys2019_pseudobulk_metadata.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 12. VALIDAÇÃO FINAL
# ============================================================

cat("\n============================================\n")
cat("11. VALIDAÇÃO FINAL\n")
cat("============================================\n\n")

cat(
  "Genes:",
  nrow(counts),
  "\n"
)

cat(
  "Núcleos:",
  ncol(counts),
  "\n"
)

cat(
  "Doadores:",
  uniqueN(cell_meta$donor_id),
  "\n"
)

cat(
  "Tipos celulares:",
  uniqueN(cell_meta$cell_type),
  "\n"
)

cat(
  "Perfis pseudobulk:",
  ncol(pseudobulk),
  "\n"
)

cat(
  "Combinações donor-celltype ausentes:",
  sum(donor_cell_qc$n_nuclei == 0),
  "\n"
)

cat("\nArquivos gerados em:\n")
cat(out_dir, "\n")

cat("\nPreparação da referência concluída.\n")
cat("Nenhuma normalização foi aplicada às contagens originais.\n")
cat("Nenhuma célula ou tipo celular foi removido nesta etapa.\n")
