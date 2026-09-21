#!/usr/bin/env Rscript

# ============================================================
# 00_auditoria.R
# Auditoria dos arquivos de expressão processados disponibilizados pelos estudos originais.
#
# Datasets:
#   GSE53697 -> Alzheimer (AD)
#   GSE64810 -> Huntington (HD)
#   GSE68719 -> Parkinson (PD)
#
# Objetivo: verificar a integridade e a estrutura dos arquivos antes de qualquer processamento analítico.
#
# Uso:
# Rscript scripts/00_auditoria.R \
#   <dataset> \
#   <expression_file> \
#   <output_report>
# ============================================================

# ------------------------------------------------------------
# 1. Pacotes
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# 2. Argumentos
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "Uso:\n",
      "Rscript scripts/00_auditoria.R ",
      "<dataset> <expression_file> <output_report>\n\n",
      "Exemplo:\n",
      "Rscript scripts/00_auditoria.R ",
      "GSE53697 ",
      "data/raw/GSE53697/GSE53697_RNAseq_AD.txt.gz ",
      "results/00_audit/GSE53697_audit.txt"
    )
  )
}

dataset <- args[1]
expression_file <- args[2]
output_report <- args[3]


# ------------------------------------------------------------
# 3. Verificar dataset
# ------------------------------------------------------------

valid_datasets <- c(
  "GSE53697",
  "GSE64810",
  "GSE68719"
)

if (!dataset %in% valid_datasets) {
  stop(
    "Dataset não reconhecido: ", dataset,
    "\nDatasets permitidos: ",
    paste(valid_datasets, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 4. Verificar arquivo
# ------------------------------------------------------------
if (!file.exists(expression_file)) {
  stop(
    "Arquivo de expressão não encontrado:\n",
    expression_file
  )
}


# ------------------------------------------------------------
# 5. Criar diretório de saída
# ------------------------------------------------------------
output_dir <- dirname(output_report)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------
# 6. Carregar arquivo
# ------------------------------------------------------------
cat("Carregando:", expression_file, "\n")

expr <- fread(
  expression_file,
  check.names = FALSE
)

cat(
  "Arquivo carregado:",
  nrow(expr), "linhas x",
  ncol(expr), "colunas\n"
)


# ------------------------------------------------------------
# 7. Configuração específica por dataset
# ------------------------------------------------------------
# Cada estudo disponibilizou seus dados processados em um formato diferente.
# Identificamos a coluna de gene, as colunas de expressão, os grupos e o número esperado de amostras

if (dataset == "GSE53697") {
  disease <- "AD"
  gene_id_col <- "GeneID"
  control_cols <- grep(
    "^C[0-9]+_raw$",
    colnames(expr),
    value = TRUE
  )

  disease_cols <- grep(
    "^A[0-9]+_raw$",
    colnames(expr),
    value = TRUE
  )

  expression_cols <- c(
    control_cols,
    disease_cols
  )

  expected_control <- 8
  expected_disease <- 9
  expected_total <- 17

  matrix_description <- paste0(
    "Raw-expression columns disponibilizadas pelos autores ",
    "(sufixo _raw); o arquivo também contém RPKM."
  )


} else if (dataset == "GSE64810") {
  disease <- "HD"
  gene_id_col <- "V1"
  control_cols <- grep(
    "^C_",
    colnames(expr),
    value = TRUE
  )

  disease_cols <- grep(
    "^H_",
    colnames(expr),
    value = TRUE
  )

  expression_cols <- c(
    control_cols,
    disease_cols
  )

  expected_control <- 49
  expected_disease <- 20
  expected_total <- 69

  matrix_description <- paste0(
    "Counts normalizados pelo DESeq2 ",
    "disponibilizados pelos autores."
  )


} else if (dataset == "GSE68719") {
  disease <- "PD"
  gene_id_col <- "EnsemblID"
  control_cols <- grep(
    "^C_",
    colnames(expr),
    value = TRUE
  )

  disease_cols <- grep(
    "^P_",
    colnames(expr),
    value = TRUE
  )

  expression_cols <- c(
    control_cols,
    disease_cols
  )

  expected_control <- 44
  expected_disease <- 29
  expected_total <- 73
  matrix_description <- paste0(
    "Counts normalizados pelo DESeq2 ",
    "disponibilizados pelos autores."
  )
}


# ------------------------------------------------------------
# 8. Verificar coluna de identificador gênico
# ------------------------------------------------------------
if (!gene_id_col %in% colnames(expr)) {
  stop(
    "Coluna de identificador gênico esperada não encontrada: ",
    gene_id_col
  )
}


# ------------------------------------------------------------
# 9. Verificar colunas de expressão
# ------------------------------------------------------------
if (length(expression_cols) == 0) {
  stop(
    "Nenhuma coluna de expressão foi identificada para ",
    dataset
  )
}

expr_matrix <- as.matrix(
  expr[, ..expression_cols]
)

storage.mode(expr_matrix) <- "numeric"

# ------------------------------------------------------------
# 10. Estatísticas básicas
# ------------------------------------------------------------
n_genes <- nrow(expr)
n_control <- length(control_cols)
n_disease <- length(disease_cols)
n_samples <- length(expression_cols)

duplicated_gene_ids <- sum(
  duplicated(expr[[gene_id_col]])
)

missing_gene_ids <- sum(
  is.na(expr[[gene_id_col]]) |
    trimws(as.character(expr[[gene_id_col]])) == ""
)

n_na <- sum(is.na(expr_matrix))

n_negative <- sum(
  expr_matrix < 0,
  na.rm = TRUE
)

n_zero <- sum(
  expr_matrix == 0,
  na.rm = TRUE
)

n_values <- sum(!is.na(expr_matrix))

zero_percentage <- if (n_values > 0) {
  100 * n_zero / n_values
} else {
  NA_real_
}


# ------------------------------------------------------------
# 11. Verificar valores inteiros
# ------------------------------------------------------------
finite_values <- expr_matrix[
  is.finite(expr_matrix)
]

if (length(finite_values) > 0) {
  integer_values <- abs(
    finite_values - round(finite_values)
  ) < 1e-8
  all_integer <- all(integer_values)
  non_integer_percentage <- 100 * mean(
    !integer_values
  )

} else {
  all_integer <- NA
  non_integer_percentage <- NA_real_
}


# ------------------------------------------------------------
# 12. Distribuição da expressão
# ------------------------------------------------------------
expression_summary <- summary(
  as.numeric(expr_matrix)
)

# ------------------------------------------------------------
# 13. Verificar número esperado de amostras
# ------------------------------------------------------------
control_ok <- n_control == expected_control
disease_ok <- n_disease == expected_disease
total_ok <- n_samples == expected_total

# ------------------------------------------------------------
# 14. Status geral
# ------------------------------------------------------------
critical_checks <- c(
  duplicated_gene_ids == 0,
  missing_gene_ids == 0,
  n_na == 0,
  n_negative == 0,
  control_ok,
  disease_ok,
  total_ok
)

audit_status <- if (all(critical_checks)) {
  "PASS"
} else {
  "CHECK"
}

# ------------------------------------------------------------
# 15. Construir relatório
# ------------------------------------------------------------
report <- c(

  "============================================================",
  "AUDITORIA DOS DADOS DE EXPRESSÃO",
  "============================================================",
  "",

  paste("Dataset:", dataset),
  paste("Doença:", disease),
  paste("Arquivo:", expression_file),
  paste("Status geral:", audit_status),

  "",
  "------------------------------------------------------------",
  "TIPO DE MATRIZ",
  "------------------------------------------------------------",
  "",

  matrix_description,

  "",
  "IMPORTANTE:",
  paste0(
    "Esta auditoria descreve a matriz como disponibilizada ",
    "pelos autores. Nenhuma normalização, transformação log2 ",
    "ou filtragem foi realizada."
  ),

  "",
  "------------------------------------------------------------",
  "DIMENSÕES",
  "------------------------------------------------------------",
  "",

  paste("Genes:", n_genes),
  paste("Amostras de expressão:", n_samples),
  paste("Controles:", n_control),
  paste(disease, ":", n_disease),

  "",
  paste(
    "Controles esperados:",
    expected_control,
    "->",
    ifelse(control_ok, "OK", "CHECK")
  ),

  paste(
    paste0(disease, " esperados:"),
    expected_disease,
    "->",
    ifelse(disease_ok, "OK", "CHECK")
  ),

  paste(
    "Total esperado:",
    expected_total,
    "->",
    ifelse(total_ok, "OK", "CHECK")
  ),

  "",
  "------------------------------------------------------------",
  "IDENTIFICADORES GÊNICOS",
  "------------------------------------------------------------",
  "",

  paste("Coluna de ID:", gene_id_col),
  paste("IDs duplicados:", duplicated_gene_ids),
  paste("IDs ausentes/vazios:", missing_gene_ids),

  "",
  "Primeiros IDs:",
  paste(
    head(as.character(expr[[gene_id_col]]), 10),
    collapse = ", "
  ),

  "",
  "------------------------------------------------------------",
  "INTEGRIDADE DA MATRIZ",
  "------------------------------------------------------------",
  "",

  paste("Valores NA:", n_na),
  paste("Valores negativos:", n_negative),
  paste("Valores iguais a zero:", n_zero),

  paste0(
    "Percentual de zeros: ",
    round(zero_percentage, 2),
    "%"
  ),

  paste(
    "Todos os valores são inteiros:",
    all_integer
  ),

  paste0(
    "Percentual de valores não inteiros: ",
    round(non_integer_percentage, 2),
    "%"
  ),

  "",
  "------------------------------------------------------------",
  "DISTRIBUIÇÃO DOS VALORES DE EXPRESSÃO",
  "------------------------------------------------------------",
  "",

  capture.output(
    print(expression_summary)
  ),

  "",
  "------------------------------------------------------------",
  "AMOSTRAS",
  "------------------------------------------------------------",
  "",

  "Controles:",
  paste(control_cols, collapse = ", "),

  "",
  paste0(disease, ":"),
  paste(disease_cols, collapse = ", "),

  "",
  "============================================================",
  paste("RESULTADO FINAL:", audit_status),
  "============================================================"
)


# ------------------------------------------------------------
# 16. Salvar relatório
# ------------------------------------------------------------
writeLines(
  report,
  con = output_report
)

# ------------------------------------------------------------
# 17. Mostrar resumo no terminal
# ------------------------------------------------------------
cat("\n")
cat("============================================\n")
cat("AUDITORIA CONCLUÍDA\n")
cat("============================================\n")
cat("Dataset:", dataset, "\n")
cat("Genes:", n_genes, "\n")
cat("Amostras:", n_samples, "\n")
cat("Control:", n_control, "\n")
cat(disease, ":", n_disease, "\n")
cat("IDs duplicados:", duplicated_gene_ids, "\n")
cat("NA:", n_na, "\n")
cat("Negativos:", n_negative, "\n")
cat("Todos inteiros:", all_integer, "\n")
cat("Status:", audit_status, "\n")
cat("\nRelatório salvo em:\n")
cat(output_report, "\n")