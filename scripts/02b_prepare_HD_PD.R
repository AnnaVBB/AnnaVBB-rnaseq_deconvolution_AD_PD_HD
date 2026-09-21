#!/usr/bin/env Rscript

# ============================================================
# 02b_prepare_HD_PD.R
#
# Prepara as matrizes de expressão processadas pelos autores
# para:
#   GSE64810 = Huntington's disease (HD)
#   GSE68719 = Parkinson's disease (PD)
#
# Entrada:
#   - DESeq2 normalized counts fornecidos pelos autores
#   - anotação NCBI para recuperação/complementação de símbolos
#
# Estratégia:
#   1. Preservar o universo de genes disponibilizado pelo estudo
#   2. Não aplicar TMM ou nova normalização
#   3. Não aplicar cutoff adicional de baixa expressão
#   4. Aplicar log2(normalized counts + 1)
#   5. Remover versão dos Ensembl IDs
#   6. Preservar identificador original
#   7. Produzir estrutura comum:
#
#      gene_id | original_id | symbol | sample1 | sample2 | ...
#
# Uso:
#
# Rscript scripts/02b_prepare_HD_PD.R \
#   <dataset> \
#   <expression_file> \
#   <annotation_file> \
#   <output_file>
#
# dataset deve ser:
#   GSE64810
#   ou
#   GSE68719
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

if (length(args) != 4) {
  stop(
    paste0(
      "\nUso:\n\n",
      "Rscript scripts/02b_prepare_HD_PD.R ",
      "<dataset> ",
      "<expression_file> ",
      "<annotation_file> ",
      "<output_file>\n"
    )
  )
}

dataset <- args[1]
expression_file <- args[2]
annotation_file <- args[3]
output_file <- args[4]


# ------------------------------------------------------------
# 3. Validar dataset
# ------------------------------------------------------------

allowed_datasets <- c(
  "GSE64810",
  "GSE68719"
)

if (!dataset %in% allowed_datasets) {
  stop(
    "Dataset inválido: ",
    dataset,
    "\nUse GSE64810 ou GSE68719."
  )
}


# ------------------------------------------------------------
# 4. Verificar arquivos
# ------------------------------------------------------------

if (!file.exists(expression_file)) {
  stop(
    "Arquivo de expressão não encontrado: ",
    expression_file
  )
}

if (!file.exists(annotation_file)) {
  stop(
    "Arquivo de anotação não encontrado: ",
    annotation_file
  )
}

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 5. Carregar expressão
# ------------------------------------------------------------

message(
  "Carregando matriz de expressão: ",
  dataset
)

expr <- fread(
  expression_file,
  data.table = FALSE
)

n_genes_initial <- nrow(expr)


# ------------------------------------------------------------
# 6. Identificar estrutura específica do dataset
# ------------------------------------------------------------

if (dataset == "GSE64810") {

  # HD:
  # primeira coluna = Ensembl ID com versão
  # demais colunas = 69 amostras

  original_id <- as.character(
    expr[[1]]
  )

  sample_cols <- colnames(expr)[-1]

  expected_samples <- 69

  study_symbol <- rep(
    NA_character_,
    nrow(expr)
  )

} else if (dataset == "GSE68719") {

  # PD:
  # EnsemblID
  # symbol
  # + 73 amostras

  required_columns <- c(
    "EnsemblID",
    "symbol"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(expr)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Colunas obrigatórias ausentes no PD: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  original_id <- as.character(
    expr$EnsemblID
  )

  study_symbol <- as.character(
    expr$symbol
  )

  sample_cols <- setdiff(
    colnames(expr),
    c("EnsemblID", "symbol")
  )

  expected_samples <- 73
}


# ------------------------------------------------------------
# 7. Validar número de amostras
# ------------------------------------------------------------

if (length(sample_cols) != expected_samples) {
  stop(
    "Esperadas ",
    expected_samples,
    " amostras para ",
    dataset,
    ", mas foram encontradas ",
    length(sample_cols),
    "."
  )
}

message(
  "Amostras encontradas: ",
  length(sample_cols)
)


# ------------------------------------------------------------
# 8. Construir matriz numérica
# ------------------------------------------------------------

expression_matrix <- as.matrix(
  expr[, sample_cols, drop = FALSE]
)

storage.mode(expression_matrix) <- "numeric"

if (anyNA(expression_matrix)) {
  stop(
    "Foram encontrados valores NA na matriz de expressão."
  )
}

if (any(expression_matrix < 0)) {
  stop(
    "Foram encontrados valores negativos na matriz."
  )
}


# ------------------------------------------------------------
# 9. Auditoria de genes all-zero
#
# Não removemos automaticamente.
# Apenas registramos, pois a estratégia da Fase 1 é preservar
# o universo de genes disponibilizado pelos autores.
# ------------------------------------------------------------

all_zero <- rowSums(
  expression_matrix > 0
) == 0

n_all_zero <- sum(all_zero)


# ------------------------------------------------------------
# 10. Transformação
#
# Os valores já são normalized counts produzidos pelo
# processamento dos autores.
#
# Portanto:
#   - NÃO aplicar TMM
#   - NÃO aplicar DESeq2 novamente
#   - NÃO aplicar voom nesta etapa
# ------------------------------------------------------------

log_expr <- log2(
  expression_matrix + 1
)


# ------------------------------------------------------------
# 11. Harmonizar Ensembl
# Exemplo: ENSG00000141510.15 --> ENSG00000141510
# ------------------------------------------------------------
gene_id <- sub(
  "\\..*$",
  "",
  original_id
)

gene_id[
  gene_id == ""
] <- NA_character_


# ------------------------------------------------------------
# 12. Carregar anotação NCBI
# ------------------------------------------------------------
message("Carregando anotação...")

annotation <- fread(
  annotation_file,
  data.table = FALSE
)

required_annotation <- c(
  "Symbol",
  "EnsemblGeneID"
)

missing_annotation <- setdiff(
  required_annotation,
  colnames(annotation)
)

if (length(missing_annotation) > 0) {
  stop(
    "Colunas ausentes na anotação: ",
    paste(missing_annotation, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 13. Preparar mapa Ensembl -> Symbol
#
# Um mesmo Ensembl pode estar associado a múltiplos símbolos
# na anotação NCBI.
#
# Regra:
#   - 1 símbolo único  -> preservar
#   - >1 símbolos      -> NA (mapeamento ambíguo)
#
# O gene_id da matriz de expressão nunca é duplicado por causa
# de uma anotação auxiliar.
# ------------------------------------------------------------

annotation_map_raw <- annotation %>%
  transmute(
    gene_id = sub(
      "\\..*$",
      "",
      as.character(EnsemblGeneID)
    ),
    annotation_symbol = as.character(Symbol)
  ) %>%
  filter(
    !is.na(gene_id),
    gene_id != ""
  ) %>%
  distinct()


annotation_map <- annotation_map_raw %>%
  group_by(gene_id) %>%
  summarise(
    n_symbols = n_distinct(
      annotation_symbol[
        !is.na(annotation_symbol) &
        annotation_symbol != ""
      ]
    ),

    annotation_symbol = {
      valid_symbols <- unique(
        annotation_symbol[
          !is.na(annotation_symbol) &
          annotation_symbol != ""
        ]
      )

      if (length(valid_symbols) == 1) {
        valid_symbols
      } else {
        NA_character_
      }
    },

    .groups = "drop"
  )


n_ambiguous_annotation <- sum(
  annotation_map$n_symbols > 1
)

# ------------------------------------------------------------
# 14. Construir informações dos genes
# ------------------------------------------------------------
gene_info <- data.frame(
  gene_id = gene_id,
  original_id = original_id,
  study_symbol = study_symbol,
  stringsAsFactors = FALSE
)

gene_info <- gene_info %>%
  left_join(
    annotation_map,
    by = "gene_id"
  )

n_ambiguous_in_dataset <- sum(
  gene_info$n_symbols > 1,
  na.rm = TRUE
)

# ------------------------------------------------------------
# 15. Escolher símbolo
# PD: prioridade para symbol fornecido pelo estudo
# HD: study_symbol é NA, então usamos anotação NCBI
# ------------------------------------------------------------

gene_info <- gene_info %>%
  mutate(
    symbol = case_when(
      !is.na(study_symbol) &
        study_symbol != "" ~ study_symbol,

      !is.na(annotation_symbol) &
        annotation_symbol != "" ~ annotation_symbol,

      TRUE ~ NA_character_
    )
  )


# ------------------------------------------------------------
# 16. Verificar se a junção aumentou o número de linhas
# Isso detecta mapeamentos Ensembl -> múltiplos símbolos antes de combinarmos a expressão.
# ------------------------------------------------------------

if (nrow(gene_info) != n_genes_initial) {
  stop(
    paste0(
      "A junção com a anotação alterou o número de linhas: ",
      n_genes_initial,
      " -> ",
      nrow(gene_info),
      ".\n",
      "Provável mapeamento Ensembl para múltiplas anotações. ",
      "Investigue antes de prosseguir."
    )
  )
}


# ------------------------------------------------------------
# 17. Montar tabela preparada
# ------------------------------------------------------------

expression_df <- as.data.frame(
  log_expr,
  check.names = FALSE
)

prepared <- bind_cols(
  gene_info %>%
    select(
      gene_id,
      original_id,
      symbol
    ),
  expression_df
)


# ------------------------------------------------------------
# 18. Auditoria dos identificadores
# ------------------------------------------------------------
n_missing_gene_id <- sum(
  is.na(prepared$gene_id)
)

n_missing_symbol <- sum(
  is.na(prepared$symbol) |
    prepared$symbol == ""
)

n_duplicated_original <- sum(
  duplicated(prepared$original_id)
)

n_duplicated_gene_id <- sum(
  duplicated(
    prepared$gene_id[
      !is.na(prepared$gene_id)
    ]
  )
)


# ------------------------------------------------------------
# 19. Variabilidade
# Genes com variância zero não são removidos aqui.
# Apenas registramos para decidir no modelo estatístico.
# ------------------------------------------------------------
gene_variance <- apply(
  log_expr,
  1,
  var
)

n_zero_variance <- sum(
  gene_variance == 0 |
    is.na(gene_variance)
)


# ------------------------------------------------------------
# 20. Resumo
# ------------------------------------------------------------
cat(
  "\n",
  "============================================\n",
  "PREPARAÇÃO DA EXPRESSÃO — ",
  dataset,
  "\n",
  "============================================\n",
  "Genes no arquivo original: ",
  n_genes_initial,
  "\n",
  "Genes preservados: ",
  nrow(prepared),
  "\n",
  "Amostras: ",
  ncol(log_expr),
  "\n",
  "Genes all-zero: ",
  n_all_zero,
  "\n",
  "Genes com variância zero: ",
  n_zero_variance,
  "\n",
  "Original IDs duplicados: ",
  n_duplicated_original,
  "\n",
  "Ensembl ausentes: ",
  n_missing_gene_id,
  "\n",
  "Ensembl duplicados: ",
  n_duplicated_gene_id,
  "\n",
  "Símbolos ausentes: ",
  n_missing_symbol,
  "\n",
  "Ensembl ambíguos na anotação NCBI: ",
  n_ambiguous_annotation,
  "\n",
  "Ensembl ambíguos presentes no dataset: ",
  n_ambiguous_in_dataset,
  "\n",
  "Transformação: log2(normalized counts + 1)\n",
  "Nova normalização aplicada: NÃO\n",
  "Filtro adicional de baixa expressão: NÃO\n",
  "============================================\n",
  sep = ""
)


# ------------------------------------------------------------
# 21. Salvar
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