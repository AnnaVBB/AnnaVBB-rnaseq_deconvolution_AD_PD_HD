#!/usr/bin/env Rscript

# ============================================================
# 03_exploratory_deg.R
#
# Análise exploratória e expressão diferencial — Fase 1
#
# Datasets:
#   GSE53697 = Alzheimer (AD)
#   GSE64810 = Huntington (HD)
#   GSE68719 = Parkinson (PD)
#
# Entrada:
#   - matriz preparada por 02a/02b
#   - metadata preparado por 01_parse_metadata.R
#
# Estratégia principal:
#   - análise separada para cada dataset
#   - modelo comum: ~ group
#   - Control como referência
#   - contraste: Disease - Control
#   - limma sobre expressão log2 transformada
#   - eBayes(trend = TRUE)
#
# Critério principal:
#   adj.P.Val < 0.05
#   |logFC| > 0.58
#
# Critério de sensibilidade / artigo comparativo:
#   P.Value < 0.05
#   |logFC| > 1
#
# Uso:
#
# Rscript scripts/03_exploratory_deg.R \
#   <dataset> \
#   <expression_file> \
#   <metadata_file> \
#   <output_dir>
#
# Exemplo:
#
# Rscript scripts/03_exploratory_deg.R \
#   GSE53697 \
#   data/processed/GSE53697_expression_prepared.csv \
#   metadata/GSE53697_metadata.csv \
#   results/01_core_original/GSE53697
#
# ============================================================


# ------------------------------------------------------------
# 1. Pacotes
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(limma)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})


# ------------------------------------------------------------
# 2. Argumentos
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    paste0(
      "\nUso:\n\n",
      "Rscript scripts/03_exploratory_deg.R ",
      "<dataset> ",
      "<expression_file> ",
      "<metadata_file> ",
      "<output_dir>\n"
    )
  )
}

dataset <- args[1]
expression_file <- args[2]
metadata_file <- args[3]
output_dir <- args[4]


# ------------------------------------------------------------
# 3. Configuração dos datasets
# ------------------------------------------------------------

dataset_config <- list(

  GSE53697 = list(
    disease = "AD",
    disease_label = "Alzheimer"
  ),

  GSE64810 = list(
    disease = "HD",
    disease_label = "Huntington"
  ),

  GSE68719 = list(
    disease = "PD",
    disease_label = "Parkinson"
  )

)

if (!dataset %in% names(dataset_config)) {
  stop(
    "Dataset inválido: ",
    dataset,
    "\nUse GSE53697, GSE64810 ou GSE68719."
  )
}

disease <- dataset_config[[dataset]]$disease
disease_label <- dataset_config[[dataset]]$disease_label


# ------------------------------------------------------------
# 4. Verificar arquivos
# ------------------------------------------------------------

if (!file.exists(expression_file)) {
  stop(
    "Arquivo de expressão não encontrado: ",
    expression_file
  )
}

if (!file.exists(metadata_file)) {
  stop(
    "Arquivo de metadata não encontrado: ",
    metadata_file
  )
}


# ------------------------------------------------------------
# 5. Criar diretórios
# ------------------------------------------------------------

deg_dir <- file.path(
  output_dir,
  "deg"
)

plot_dir <- file.path(
  output_dir,
  "plots"
)

dir.create(
  deg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 6. Carregar dados
# ------------------------------------------------------------

message(
  "Carregando dados — ",
  dataset
)

expr <- read_csv(
  expression_file,
  show_col_types = FALSE
)

meta <- read_csv(
  metadata_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 7. Validar estrutura da expressão
# ------------------------------------------------------------

required_expr <- c(
  "gene_id",
  "original_id",
  "symbol"
)

missing_expr <- setdiff(
  required_expr,
  colnames(expr)
)

if (length(missing_expr) > 0) {
  stop(
    "Colunas ausentes na matriz de expressão: ",
    paste(missing_expr, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 8. Validar metadata
# ------------------------------------------------------------

required_meta <- c(
  "sample_id",
  "group"
)

missing_meta <- setdiff(
  required_meta,
  colnames(meta)
)

if (length(missing_meta) > 0) {
  stop(
    "Colunas ausentes no metadata: ",
    paste(missing_meta, collapse = ", ")
  )
}

if (anyDuplicated(meta$sample_id) > 0) {
  stop(
    "Existem sample_id duplicados no metadata."
  )
}


# ------------------------------------------------------------
# 9. Identificar amostras
# ------------------------------------------------------------

sample_cols <- setdiff(
  colnames(expr),
  required_expr
)

if (length(sample_cols) == 0) {
  stop(
    "Nenhuma coluna de amostra encontrada."
  )
}


# ------------------------------------------------------------
# 10. Validar expressão x metadata
# ------------------------------------------------------------

expression_without_metadata <- setdiff(
  sample_cols,
  meta$sample_id
)

metadata_without_expression <- setdiff(
  meta$sample_id,
  sample_cols
)

if (
  length(expression_without_metadata) > 0 ||
  length(metadata_without_expression) > 0
) {

  stop(
    paste0(
      "Expressão e metadata possuem amostras diferentes.\n",
      "Expressão sem metadata: ",
      paste(
        expression_without_metadata,
        collapse = ", "
      ),
      "\n",
      "Metadata sem expressão: ",
      paste(
        metadata_without_expression,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 11. Reordenar metadata explicitamente
#
# Mesmo que já esteja na mesma ordem, fazemos isso como
# proteção contra alterações futuras.
# ------------------------------------------------------------

meta <- meta[
  match(
    sample_cols,
    meta$sample_id
  ),
  ,
  drop = FALSE
]

if (!identical(
  sample_cols,
  meta$sample_id
)) {
  stop(
    "Falha ao alinhar expressão e metadata."
  )
}


# ------------------------------------------------------------
# 12. Construir matriz de expressão
# ------------------------------------------------------------

expression_matrix <- as.matrix(
  expr[
    ,
    sample_cols,
    drop = FALSE
  ]
)

storage.mode(
  expression_matrix
) <- "numeric"

rownames(
  expression_matrix
) <- seq_len(
  nrow(expression_matrix)
)


# ------------------------------------------------------------
# 13. Validação numérica
# ------------------------------------------------------------

if (anyNA(expression_matrix)) {
  stop(
    "Existem valores NA na matriz de expressão."
  )
}

if (any(!is.finite(expression_matrix))) {
  stop(
    "Existem valores não finitos na matriz."
  )
}

if (any(expression_matrix < 0)) {
  stop(
    paste0(
      "Foram encontrados valores negativos. ",
      "A matriz preparada deveria estar em escala ",
      "log2(x + 1)."
    )
  )
}


# ------------------------------------------------------------
# 14. Verificar variância
# ------------------------------------------------------------

gene_variance <- apply(
  expression_matrix,
  1,
  var
)

zero_variance <- (
  is.na(gene_variance) |
  gene_variance == 0
)

n_zero_variance <- sum(
  zero_variance
)

if (n_zero_variance > 0) {

  message(
    "Removendo ",
    n_zero_variance,
    " genes com variância zero."
  )

  expression_matrix <- expression_matrix[
    !zero_variance,
    ,
    drop = FALSE
  ]

  expr <- expr[
    !zero_variance,
    ,
    drop = FALSE
  ]
}


# ------------------------------------------------------------
# 15. Validar grupos
# ------------------------------------------------------------

observed_groups <- unique(
  meta$group
)

if (!"Control" %in% observed_groups) {
  stop(
    "Grupo Control não encontrado."
  )
}

if (!disease %in% observed_groups) {
  stop(
    "Grupo ",
    disease,
    " não encontrado."
  )
}

unexpected_groups <- setdiff(
  observed_groups,
  c(
    "Control",
    disease
  )
)

if (length(unexpected_groups) > 0) {
  stop(
    "Grupos inesperados encontrados: ",
    paste(
      unexpected_groups,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 16. Definir fator
#
# Control é explicitamente a referência.
# ------------------------------------------------------------

meta$group <- factor(
  meta$group,
  levels = c(
    "Control",
    disease
  )
)

cat(
  "\nGrupos utilizados:\n"
)

print(
  table(meta$group)
)


# ============================================================
# PARTE A — QC EXPLORATÓRIO
# ============================================================


# ------------------------------------------------------------
# 17. PCA
#
# prcomp espera:
#   linhas   = amostras
#   colunas  = genes
#
# scale. = FALSE porque os genes já estão em escala log2 e
# não queremos dar peso igual artificialmente a genes de
# baixa e alta variabilidade.
# ------------------------------------------------------------

message(
  "Calculando PCA..."
)

pca <- prcomp(
  t(expression_matrix),
  center = TRUE,
  scale. = FALSE
)

variance_explained <- (
  pca$sdev^2 /
  sum(pca$sdev^2)
) * 100


# ------------------------------------------------------------
# 18. Tabela PCA
# ------------------------------------------------------------

pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = meta$group,
  stringsAsFactors = FALSE
)

write_csv(
  pca_df,
  file.path(
    plot_dir,
    "pca_coordinates.csv"
  )
)


# ------------------------------------------------------------
# 19. Plot PCA
# ------------------------------------------------------------
group_colors <- c(
  "Control" = "blue",
  "AD" = "red",
  "HD" = "red",
  "PD" = "red"
)

pca_plot <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = group,
    shape = group
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  labs(
    title = paste0(
      dataset,
      " — ",
      disease_label,
      " vs Control"
    ),
    subtitle = "PCA da matriz de expressão preparada",
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 1),
      "%)"
    ),
    color = "Grupo",
    shape = "Grupo"
  ) +
  theme_bw()

ggsave(
  filename = file.path(
    plot_dir,
    "pca.png"
  ),
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)
if (!file.exists(file.path(plot_dir, "pca.png"))) {
  stop("Falha ao produzir pca.png")
}

message(
  "PCA salvo em: ",
  file.path(plot_dir, "pca.png")
)

# ------------------------------------------------------------
# 20. Distribuição por amostra
#
# Boxplot em formato longo.
# ------------------------------------------------------------

distribution_df <- data.frame(
  sample_id = rep(
    colnames(expression_matrix),
    each = nrow(expression_matrix)
  ),
  expression = as.vector(
    expression_matrix
  ),
  stringsAsFactors = FALSE
)

distribution_df$group <- meta$group[
  match(
    distribution_df$sample_id,
    meta$sample_id
  )
]

distribution_plot <- ggplot(
  distribution_df,
  aes(
    x = sample_id,
    y = expression,
    group = sample_id
  )
) +
  geom_boxplot(
    outlier.size = 0.2
  ) +
  labs(
    title = paste0(
      dataset,
      " — distribuição da expressão"
    ),
    x = "Amostra",
    y = "Expressão log2"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6
    )
  )

ggsave(
  filename = file.path(
    plot_dir,
    "expression_distribution.png"
  ),
  plot = distribution_plot,
  width = 14,
  height = 6,
  dpi = 300
)


# ============================================================
# PARTE B — EXPRESSÃO DIFERENCIAL
# ============================================================


# ------------------------------------------------------------
# 21. Matriz de design
#
# Sem intercepto:
#
# Control
# Disease
#
# Isso torna o contraste explícito:
#
# Disease - Control
# ------------------------------------------------------------

design <- model.matrix(
  ~ 0 + group,
  data = meta
)

colnames(design) <- levels(
  meta$group
)

rownames(design) <- meta$sample_id


# ------------------------------------------------------------
# 22. Contraste
# ------------------------------------------------------------

contrast_string <- paste0(
  disease,
  "-Control"
)

contrast_matrix <- makeContrasts(
  contrasts = contrast_string,
  levels = design
)


# ------------------------------------------------------------
# 23. limma
#
# Não usamos:
#   - DGEList
#   - calcNormFactors
#   - voom
#
# porque a entrada desta fase já é expressão processada pelos
# estudos e transformada para log2.
#
# trend = TRUE permite tendência média-variância no eBayes.
# ------------------------------------------------------------

message(
  "Executando limma: ",
  contrast_string
)

fit <- lmFit(
  expression_matrix,
  design
)

fit <- contrasts.fit(
  fit,
  contrast_matrix
)

fit <- eBayes(
  fit,
  trend = TRUE
)


# ------------------------------------------------------------
# 24. Extrair tabela completa
# ------------------------------------------------------------

deg <- topTable(
  fit,
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)


# ------------------------------------------------------------
# 25. Verificar correspondência de linhas
# ------------------------------------------------------------

if (nrow(deg) != nrow(expr)) {
  stop(
    "Número de genes no resultado limma não corresponde ",
    "à matriz de expressão."
  )
}


# ------------------------------------------------------------
# 26. Acrescentar identificadores
# ------------------------------------------------------------

deg <- data.frame(
  gene_id = expr$gene_id,
  original_id = expr$original_id,
  symbol = expr$symbol,
  deg,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 27. Classificação dos resultados
#
# Critério principal:
#   FDR < 0.05
#   |logFC| > 0.58
#
# Critério artigo:
#   p < 0.05
#   |logFC| > 1
# ------------------------------------------------------------

deg <- deg %>%
  mutate(

    significant_primary = (
      adj.P.Val < 0.05 &
      abs(logFC) > 0.58
    ),

    significant_article = (
      P.Value < 0.05 &
      abs(logFC) > 1
    ),

    direction_primary = case_when(
      significant_primary &
        logFC > 0 ~ "Up",

      significant_primary &
        logFC < 0 ~ "Down",

      TRUE ~ "NS"
    ),

    direction_article = case_when(
      significant_article &
        logFC > 0 ~ "Up",

      significant_article &
        logFC < 0 ~ "Down",

      TRUE ~ "NS"
    )
  )


# ------------------------------------------------------------
# 28. Ordenar por FDR
# ------------------------------------------------------------

deg_sorted <- deg %>%
  arrange(
    adj.P.Val,
    P.Value
  )


# ------------------------------------------------------------
# 29. Salvar tabela completa
# ------------------------------------------------------------

write_csv(
  deg_sorted,
  file.path(
    deg_dir,
    "deg_all.csv"
  ),
  na = "NA"
)


# ------------------------------------------------------------
# 30. Salvar DEGs principais
# ------------------------------------------------------------

deg_primary <- deg_sorted %>%
  filter(
    significant_primary
  )

write_csv(
  deg_primary,
  file.path(
    deg_dir,
    "deg_significant_primary.csv"
  ),
  na = "NA"
)


# ------------------------------------------------------------
# 31. Salvar DEGs segundo critério do artigo
# ------------------------------------------------------------

deg_article <- deg_sorted %>%
  filter(
    significant_article
  )

write_csv(
  deg_article,
  file.path(
    deg_dir,
    "deg_significant_article.csv"
  ),
  na = "NA"
)


# ============================================================
# PARTE C — VISUALIZAÇÃO DOS DEGs
# ============================================================


# ------------------------------------------------------------
# 32. Volcano
# Classificação:
#   Up   = FDR < 0.05 e logFC > 0.58
#   Down = FDR < 0.05 e logFC < -0.58
#   NS   = não significativo
# ------------------------------------------------------------
volcano_df <- deg %>%
  mutate(
    minus_log10_fdr = -log10(
      pmax(
        adj.P.Val,
        .Machine$double.xmin
      )
    ),
    volcano_class = factor(
      direction_primary,
      levels = c("Down", "NS", "Up")
    )
  )

volcano_colors <- c(
  "Down" = "blue",
  "NS"   = "grey70",
  "Up"   = "red"
)

volcano_plot <- ggplot(
  volcano_df,
  aes(
    x = logFC,
    y = minus_log10_fdr,
    color = volcano_class
  )
) +
  geom_point(
    alpha = 0.65,
    size = 1.5
  ) +
  scale_color_manual(
    values = volcano_colors,
    drop = FALSE
  ) +
  geom_vline(
    xintercept = c(-0.58, 0.58),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = paste0(
      dataset,
      " — ",
      disease_label,
      " vs Control"
    ),
    subtitle = "Critério principal: FDR < 0.05 e |logFC| > 0.58",
    x = "logFC",
    y = "-log10(FDR)",
    color = "Classificação"
  ) +
  theme_bw()

ggsave(
  filename = file.path(
    plot_dir,
    "volcano_primary.png"
  ),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 32b. Volcano — critério do artigo
#
# Critério:
#   p < 0.05
#   |logFC| > 1
# ------------------------------------------------------------

volcano_article_df <- deg %>%
  mutate(
    minus_log10_p = -log10(
      pmax(
        P.Value,
        .Machine$double.xmin
      )
    ),
    volcano_class = factor(
      direction_article,
      levels = c("Down", "NS", "Up")
    )
  )

volcano_article_plot <- ggplot(
  volcano_article_df,
  aes(
    x = logFC,
    y = minus_log10_p,
    color = volcano_class
  )
) +
  geom_point(
    alpha = 0.65,
    size = 1.5
  ) +
  scale_color_manual(
    values = volcano_colors,
    drop = FALSE
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = paste0(
      dataset,
      " — ",
      disease_label,
      " vs Control"
    ),
    subtitle = "Critério comparativo: p < 0.05 e |logFC| > 1",
    x = "logFC",
    y = "-log10(p-value)",
    color = "Classificação"
  ) +
  theme_bw()

ggsave(
  filename = file.path(
    plot_dir,
    "volcano_article.png"
  ),
  plot = volcano_article_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 33. Heatmap exploratório — top 50 genes mais variáveis
#
# Este heatmap NÃO representa DEGs.
# É utilizado apenas para exploração/QC das amostras.
# ------------------------------------------------------------

gene_variance_final <- apply(
  expression_matrix,
  1,
  var
)

top_variable_indices <- order(
  gene_variance_final,
  decreasing = TRUE
)[
  seq_len(
    min(
      50,
      length(gene_variance_final)
    )
  )
]

variable_matrix <- expression_matrix[
  top_variable_indices,
  ,
  drop = FALSE
]

variable_info <- expr[
  top_variable_indices,
  ,
  drop = FALSE
]

variable_labels <- ifelse(
  !is.na(variable_info$symbol) &
    variable_info$symbol != "",
  variable_info$symbol,
  ifelse(
    !is.na(variable_info$gene_id) &
      variable_info$gene_id != "",
    variable_info$gene_id,
    as.character(variable_info$original_id)
  )
)

rownames(variable_matrix) <- make.unique(
  as.character(variable_labels)
)

annotation_col <- data.frame(
  Group = meta$group
)

rownames(annotation_col) <- meta$sample_id

pheatmap(
  variable_matrix,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 7,
  border_color = NA,
  filename = file.path(
    plot_dir,
    "heatmap_top50_variable.png"
  ),
  width = 10,
  height = 10
)

# ------------------------------------------------------------
# 34. Heatmap dos DEGs principais 
#
# Para evitar heatmaps enormes:
#   - selecionar no máximo 50 DEGs
#   - priorizar menor FDR
#
# Se houver menos de 2 DEGs, não produzir heatmap.
# ------------------------------------------------------------

if (nrow(deg_primary) >= 2) {
  n_heatmap <- min(
    50,
    nrow(deg_primary)
  )

  heatmap_genes <- deg_primary %>%
    arrange(
      adj.P.Val,
      desc(abs(logFC))
    ) %>%
    slice_head(
      n = n_heatmap
    )

  heatmap_indices <- match(
    heatmap_genes$original_id,
    expr$original_id
  )

  heatmap_matrix <- expression_matrix[
    heatmap_indices,
    ,
    drop = FALSE
  ]

  heatmap_labels <- ifelse(
    !is.na(heatmap_genes$symbol) &
      heatmap_genes$symbol != "",
    heatmap_genes$symbol,
    heatmap_genes$gene_id
  )

  # Evitar rownames repetidos no heatmap
  heatmap_labels <- make.unique(
    as.character(heatmap_labels)
  )

  rownames(
    heatmap_matrix
  ) <- heatmap_labels

  annotation_col <- data.frame(
    Group = meta$group
  )

  rownames(
    annotation_col
  ) <- meta$sample_id

  pheatmap(
    heatmap_matrix,
    scale = "row",
    annotation_col = annotation_col,
    show_colnames = FALSE,
    fontsize_row = 7,
    border_color = NA,
    filename = file.path(
      plot_dir,
      "heatmap_top50_primary.png"
    ),
    width = 10,
    height = 10
  )

} else {

  message(
    "Heatmap não produzido: menos de 2 DEGs ",
    "pelo critério principal."
  )
}


# ============================================================
# PARTE D — RESUMO
# ============================================================
# ------------------------------------------------------------
# 34. Contagens
# ------------------------------------------------------------

n_primary <- sum(
  deg$significant_primary
)

n_primary_up <- sum(
  deg$direction_primary == "Up"
)

n_primary_down <- sum(
  deg$direction_primary == "Down"
)

n_article <- sum(
  deg$significant_article
)

n_article_up <- sum(
  deg$direction_article == "Up"
)

n_article_down <- sum(
  deg$direction_article == "Down"
)


# ------------------------------------------------------------
# 35. Menores valores estatísticos
# ------------------------------------------------------------

min_p <- min(
  deg$P.Value,
  na.rm = TRUE
)

min_fdr <- min(
  deg$adj.P.Val,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 36. Salvar resumo
# ------------------------------------------------------------

summary_df <- data.frame(

  dataset = dataset,

  disease = disease,

  n_samples = ncol(
    expression_matrix
  ),

  n_control = sum(
    meta$group == "Control"
  ),

  n_disease = sum(
    meta$group == disease
  ),

  n_genes_tested = nrow(
    expression_matrix
  ),

  primary_total = n_primary,

  primary_up = n_primary_up,

  primary_down = n_primary_down,

  article_total = n_article,

  article_up = n_article_up,

  article_down = n_article_down,

  min_p_value = min_p,

  min_fdr = min_fdr,

  stringsAsFactors = FALSE
)

write_csv(
  summary_df,
  file.path(
    deg_dir,
    "deg_summary.csv"
  )
)


# ------------------------------------------------------------
# 37. Resumo no terminal
# ------------------------------------------------------------

cat(
  "\n",
  "============================================\n",
  "ANÁLISE DIFERENCIAL — ",
  dataset,
  "\n",
  "============================================\n",
  "Doença: ",
  disease_label,
  "\n",
  "Contraste: ",
  disease,
  " - Control\n",
  "Amostras: ",
  ncol(expression_matrix),
  "\n",
  "  Control: ",
  sum(meta$group == "Control"),
  "\n",
  "  ",
  disease,
  ": ",
  sum(meta$group == disease),
  "\n",
  "Genes testados: ",
  nrow(expression_matrix),
  "\n",
  "\n",
  "Critério principal:\n",
  "  FDR < 0.05 e |logFC| > 0.58\n",
  "  Total: ",
  n_primary,
  "\n",
  "  Up: ",
  n_primary_up,
  "\n",
  "  Down: ",
  n_primary_down,
  "\n",
  "\n",
  "Critério artigo:\n",
  "  p < 0.05 e |logFC| > 1\n",
  "  Total: ",
  n_article,
  "\n",
  "  Up: ",
  n_article_up,
  "\n",
  "  Down: ",
  n_article_down,
  "\n",
  "\n",
  "Menor p-value: ",
  format(
    min_p,
    scientific = TRUE
  ),
  "\n",
  "Menor FDR: ",
  format(
    min_fdr,
    scientific = TRUE
  ),
  "\n",
  "============================================\n",
  sep = ""
)

message(
  "\nResultados salvos em:\n",
  output_dir
)