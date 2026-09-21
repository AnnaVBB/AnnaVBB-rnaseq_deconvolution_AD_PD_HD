#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(edgeR)
})

# Este script harmoniza a matriz de contagens com a tabela de anotação utilizando GeneID como identificador principal.
# IMPORTANTE: GeneID é mantido como identificador primário durante as análises estatísticas. O símbolo gênico (Symbol) é utilizado posteriormente para interpretação biológica.
# O script também realiza a filtragem de baixa expressão: CPM >= 1 em pelo menos 20% das amostras.

# ----------------------------------------------------------
# 0. Receber argumentos do terminal
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Uso: Rscript 02_harmonize_annotations.R ",
    "<counts.tsv.gz> <annot.tsv.gz> <out_rds>"
  )
}
counts_file <- args[1]
annot_file  <- args[2]
output_rds  <- args[3]


# ------------------------------------------------------------
# 1. Carregar contagens e anotação
# ------------------------------------------------------------
# fread consegue ler diretamente arquivos .gz, portanto não
# precisamos utilizar zcat para realizar a descompressão.

counts_df <- fread(
  counts_file,
  data.table = FALSE
)

annot_df <- fread(
  annot_file,
  data.table = FALSE
)

# ------------------------------------------------------------
# 2. Padronizar a coluna de identificação gênica
# ------------------------------------------------------------
# A primeira coluna das duas tabelas corresponde ao identificador
# do gene. Padronizamos seu nome para GeneID para facilitar
# as operações seguintes.
colnames(counts_df)[1] <- "GeneID"
colnames(annot_df)[1]  <- "GeneID"


# Converter GeneID para texto evita problemas de compatibilidade
# durante a comparação entre counts e anotação.
counts_df$GeneID <- as.character(counts_df$GeneID)
annot_df$GeneID  <- as.character(annot_df$GeneID)


# ------------------------------------------------------------
# 3. Verificações básicas de integridade
# ------------------------------------------------------------
# GeneIDs duplicados na matriz de counts poderiam causar
# ambiguidades durante as análises.
if (anyDuplicated(counts_df$GeneID) > 0) {
  stop("Foram encontrados GeneIDs duplicados na matriz de counts.")
}


# A anotação também deve possuir uma única linha por GeneID.
if (anyDuplicated(annot_df$GeneID) > 0) {
  warning(
    "Foram encontrados GeneIDs duplicados na anotação. ",
    "Será mantida a primeira ocorrência de cada GeneID."
  )

  annot_df <- annot_df %>%
    distinct(GeneID, .keep_all = TRUE)
}

# ------------------------------------------------------------
# 4. Identificar as colunas das amostras
# ------------------------------------------------------------
# Todas as colunas depois de GeneID correspondem às amostras GSM.
sample_cols <- setdiff(
  colnames(counts_df),
  "GeneID"
)

if (length(sample_cols) == 0) {
  stop("Nenhuma coluna de amostra foi encontrada na matriz de counts.")
}

# ------------------------------------------------------------
# 5. Criar matriz numérica de counts
# ------------------------------------------------------------
# Retiramos temporariamente GeneID e mantemos apenas os valores
# de expressão para realizar os cálculos de CPM.
counts_mat <- as.matrix(
  counts_df[, sample_cols, drop = FALSE]
)

# Garantir que os valores da matriz sejam numéricos.
storage.mode(counts_mat) <- "numeric"

# GeneID é colocado como nome das linhas para preservar a
# identidade de cada gene durante a filtragem.
rownames(counts_mat) <- counts_df$GeneID

# ------------------------------------------------------------
# 6. Verificar valores inválidos
# ------------------------------------------------------------
# Matrizes de contagens não devem conter valores ausentes
# nem valores negativos.

if (anyNA(counts_mat)) {
  stop("A matriz de counts contém valores ausentes (NA).")
}

if (any(counts_mat < 0)) {
  stop("A matriz de counts contém valores negativos.")
}

# ------------------------------------------------------------
# 7. Filtragem de genes com baixa expressão
# ------------------------------------------------------------
# Criamos um objeto DGEList para calcular CPM utilizando edgeR.
dge <- edgeR::DGEList(
  counts = counts_mat
)
cpm_mat <- edgeR::cpm(dge)

# O critério definido para o projeto é:
# CPM >= 1 em pelo menos 20% das amostras.
# ceiling() transforma 20% do número total de amostras no menor número inteiro de amostras que satisfaz esse percentual.

min_samples <- ceiling(
  0.20 * ncol(dge)
)

keep <- rowSums(
  cpm_mat >= 1,
  na.rm = TRUE
) >= min_samples


# Aplicar a filtragem à matriz de counts.
mat_filtered <- counts_mat[
  keep,
  ,
  drop = FALSE
]

cat(
  "[INFO] Genes antes da filtragem:",
  nrow(counts_mat),
  "\n"
)

cat(
  "[INFO] Critério de filtragem: CPM >= 1 em pelo menos",
  min_samples,
  "amostras\n"
)

cat(
  "[INFO] Genes mantidos após filtragem:",
  nrow(mat_filtered),
  "\n"
)

cat(
  "[INFO] Genes removidos:",
  sum(!keep),
  "\n"
)

# ------------------------------------------------------------
# 8. Harmonizar GeneID com a anotação
# ------------------------------------------------------------
# Selecionamos somente os genes que permaneceram após a
# filtragem de baixa expressão.

retained_gene_ids <- rownames(
  mat_filtered
)

# A anotação é filtrada para os GeneIDs retidos. Não removemos genes apenas porque Symbol está ausente.
# O GeneID continua sendo o identificador primário da análise.

annotation_filtered <- annot_df %>%
  filter(
    GeneID %in% retained_gene_ids
  )


# Reordenamos a anotação na mesma ordem dos genes presentes na matriz de expressão.
annotation_filtered <- annotation_filtered[
  match(
    retained_gene_ids,
    annotation_filtered$GeneID
  ),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 9. Verificar genes sem correspondência na anotação
# ------------------------------------------------------------
n_missing_annotation <- sum(
  is.na(annotation_filtered$GeneID)
)

cat(
  "[INFO] Genes sem correspondência na anotação:",
  n_missing_annotation,
  "\n"
)

# ------------------------------------------------------------
# 10. Salvar objeto harmonizado
# ------------------------------------------------------------
# O arquivo RDS contém dois componentes:
# counts: matriz filtrada de expressão, com GeneID nas linhas.
# annotation: tabela de anotação correspondente aos mesmos GeneIDs.
# Não realizamos normalização TMM neste script. A normalização será realizada posteriormente no modelo estatístico de expressão diferencial.

res_list <- list(
  counts = mat_filtered,
  annotation = annotation_filtered
)


# Criar o diretório de saída caso ainda não exista.
dir.create(
  dirname(output_rds),
  recursive = TRUE,
  showWarnings = FALSE
)


saveRDS(
  res_list,
  file = output_rds
)


cat(
  "[SUCCESS] Matriz harmonizada salva em:",
  output_rds,
  "\n"
)

cat(
  "[INFO] Genes retidos pós-filtro:",
  nrow(mat_filtered),
  "\n"
)