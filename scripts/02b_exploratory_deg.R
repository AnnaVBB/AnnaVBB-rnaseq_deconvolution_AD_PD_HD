#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tibble)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# ------------------------------------------------------------
# 0. Argumentos
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Uso: Rscript 02b_exploratory_deg.R ",
    "<counts.tsv.gz> <metadata.csv> <out_prefix>"
  )
}

tsv_file   <- args[1]
meta_file  <- args[2]
out_prefix <- args[3]

cat(
  "[INFO] Processando dataset:",
  tsv_file,
  "\n"
)

# Criar pasta de saída caso não exista
dir.create(
  dirname(out_prefix),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Carregar dados brutos e metadata
# ------------------------------------------------------------

counts <- read.delim(
  gzfile(tsv_file),
  row.names = 1,
  check.names = FALSE
)

metadata <- readr::read_csv(
  meta_file,
  show_col_types = FALSE
)

# Garantir que GeneIDs sejam tratados como texto
rownames(counts) <- as.character(rownames(counts))


# ------------------------------------------------------------
# 1.1 Verificações básicas
# ------------------------------------------------------------

if (!"sample_id" %in% colnames(metadata)) {
  stop("A coluna 'sample_id' não foi encontrada no metadata.")
}

if (!"group" %in% colnames(metadata)) {
  stop("A coluna 'group' não foi encontrada no metadata.")
}

if (anyDuplicated(rownames(counts)) > 0) {
  stop("Foram encontrados GeneIDs duplicados na matriz de contagens.")
}

if (anyDuplicated(metadata$sample_id) > 0) {
  stop("Foram encontrados sample_id duplicados no metadata.")
}

if (anyNA(counts)) {
  stop("A matriz de contagens contém valores NA.")
}

if (any(counts < 0)) {
  stop("A matriz de contagens contém valores negativos.")
}


# ------------------------------------------------------------
# 1.2 Alinhar amostras entre counts e metadata
# ------------------------------------------------------------

common_samples <- intersect(
  colnames(counts),
  metadata$sample_id
)

if (length(common_samples) < 2) {
  stop(
    "Menos de 2 amostras em comum entre counts e metadata."
  )
}

counts <- counts[
  ,
  common_samples,
  drop = FALSE
]

metadata <- metadata %>%
  dplyr::filter(sample_id %in% common_samples) %>%
  dplyr::slice(
    match(common_samples, sample_id)
  )

if (!identical(
  colnames(counts),
  metadata$sample_id
)) {
  stop(
    "A ordem das amostras em counts e metadata não coincide."
  )
}


# ------------------------------------------------------------
# 1.3 Definir comparação Control vs Disease
# ------------------------------------------------------------

disease_groups <- unique(
  metadata$group[
    !is.na(metadata$group) &
      metadata$group != "Control"
  ]
)

if (length(disease_groups) != 1) {
  stop(
    "Esperado exatamente um grupo de doença além de Control. ",
    "Encontrado: ",
    paste(
      disease_groups,
      collapse = ", "
    )
  )
}

disease_label <- disease_groups[1]

metadata$group <- factor(
  metadata$group,
  levels = c(
    "Control",
    disease_label
  )
)

if (anyNA(metadata$group)) {
  stop(
    "Existem amostras com grupo inválido após a criação do fator."
  )
}

group_counts <- table(metadata$group)

if (any(group_counts < 2)) {
  stop(
    "Cada grupo deve possuir pelo menos 2 amostras. ",
    "Distribuição encontrada: ",
    paste(
      names(group_counts),
      group_counts,
      sep = "=",
      collapse = ", "
    )
  )
}

cat(
  "[INFO] Amostras:",
  ncol(counts),
  "| Grupo Controle vs",
  disease_label,
  "\n"
)

cat(
  "[INFO] Distribuição dos grupos:",
  paste(
    names(group_counts),
    group_counts,
    sep = "=",
    collapse = ", "
  ),
  "\n"
)


# ------------------------------------------------------------
# 2. Criar objeto DGEList
# ------------------------------------------------------------

dge <- edgeR::DGEList(
  counts = counts
)

cat(
  "[INFO] Genes antes da filtragem:",
  nrow(dge),
  "\n"
)


# ------------------------------------------------------------
# 2.1 Filtragem por expressão
#
# Critério metodológico:
# CPM >= 1 em pelo menos 20% das amostras
# ------------------------------------------------------------

cpm_mat <- edgeR::cpm(dge)

min_samples <- ceiling(
  0.20 * ncol(dge)
)

keep <- rowSums(
  cpm_mat >= 1,
  na.rm = TRUE
) >= min_samples

cat(
  "[INFO] Critério de filtragem: CPM >= 1 em pelo menos",
  min_samples,
  "amostras\n"
)

cat(
  "[INFO] Genes mantidos após filtragem:",
  sum(keep),
  "\n"
)

cat(
  "[INFO] Genes removidos:",
  sum(!keep),
  "\n"
)

if (sum(keep) < 2) {
  stop(
    "Menos de 2 genes permaneceram após a filtragem por CPM."
  )
}

dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]


# ------------------------------------------------------------
# 2.2 Normalização TMM
# ------------------------------------------------------------

dge <- edgeR::calcNormFactors(
  dge,
  method = "TMM"
)

cat(
  "[INFO] Normalização TMM concluída.\n"
)


# ------------------------------------------------------------
# 3. Design experimental
# ------------------------------------------------------------

design <- model.matrix(
  ~ group,
  data = metadata
)

colnames(design) <- c(
  "Intercept",
  "Disease"
)

cat(
  "[INFO] Design experimental criado.\n"
)


# ------------------------------------------------------------
# 4. limma-voom + expressão diferencial
# ------------------------------------------------------------

v <- limma::voom(
  dge,
  design,
  plot = FALSE
)

fit <- limma::lmFit(
  v,
  design
)

fit <- limma::eBayes(
  fit
)

res_table <- limma::topTable(
  fit,
  coef = "Disease",
  number = Inf,
  sort.by = "P"
) %>%
  tibble::rownames_to_column(
    "gene_id"
  ) %>%
  tibble::as_tibble()

# Garantir GeneID como character para posterior join
res_table$gene_id <- as.character(
  res_table$gene_id
)

cat(
  "[INFO] Genes testados na análise diferencial:",
  nrow(res_table),
  "\n"
)


# ------------------------------------------------------------
# 5. Adicionar símbolo do gene
# ------------------------------------------------------------

annotation_file <-
  "data/raw/Human.GRCh38.p13.annot.tsv.gz"

if (file.exists(annotation_file)) {

  annotation <- read.delim(
    gzfile(annotation_file),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(annotation)[1] <- "GeneID"

  annotation$GeneID <- as.character(
    annotation$GeneID
  )

  if ("Symbol" %in% colnames(annotation)) {

    annotation_small <- annotation %>%
      dplyr::select(
        GeneID,
        Symbol
      ) %>%
      dplyr::rename(
        gene_id = GeneID,
        gene_symbol = Symbol
      ) %>%
      dplyr::distinct(
        gene_id,
        .keep_all = TRUE
      )

    res_table <- res_table %>%
      dplyr::left_join(
        annotation_small,
        by = "gene_id"
      )

  } else {

    warning(
      "A coluna 'Symbol' não foi encontrada no arquivo de anotação."
    )

    res_table$gene_symbol <- NA_character_
  }

} else {

  warning(
    "Arquivo de anotação não encontrado: ",
    annotation_file
  )

  res_table$gene_symbol <- NA_character_
}


# GeneID como fallback visual caso Symbol não esteja disponível
res_table <- res_table %>%
  dplyr::mutate(
    gene_symbol_display = dplyr::if_else(
      is.na(gene_symbol) |
        gene_symbol == "",
      gene_id,
      gene_symbol
    )
  )


# ------------------------------------------------------------
# 6. Classificar genes significativos
#
# Critérios:
# FDR < 0.05
# |log2FC| > 0.58
# ------------------------------------------------------------

res_table <- res_table %>%
  dplyr::mutate(
    significant =
      !is.na(adj.P.Val) &
      adj.P.Val < 0.05 &
      abs(logFC) > 0.58
  )

n_significant <- sum(
  res_table$significant,
  na.rm = TRUE
)

n_up <- sum(
  res_table$significant &
    res_table$logFC > 0,
  na.rm = TRUE
)

n_down <- sum(
  res_table$significant &
    res_table$logFC < 0,
  na.rm = TRUE
)

cat(
  "[INFO] DEGs significativos (FDR < 0.05 e |log2FC| > 0.58):",
  n_significant,
  "\n"
)

cat(
  "[INFO] Upregulated:",
  n_up,
  "| Downregulated:",
  n_down,
  "\n"
)

cat(
  "[INFO] Menor P-value:",
  min(
    res_table$P.Value,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "[INFO] Menor FDR:",
  min(
    res_table$adj.P.Val,
    na.rm = TRUE
  ),
  "\n"
)


# ------------------------------------------------------------
# 7. Salvar resultados DEG
# ------------------------------------------------------------

deg_out_csv <- paste0(
  out_prefix,
  "_deg_results.csv"
)

readr::write_csv(
  res_table,
  deg_out_csv
)

cat(
  "[SUCCESS] Resultados DEG salvos em:",
  deg_out_csv,
  "\n"
)


# ------------------------------------------------------------
# 8. PCA
#
# Utiliza a matriz de expressão transformada pelo voom
# após a filtragem por expressão.
# ------------------------------------------------------------

pca <- prcomp(
  t(v$E),
  scale. = TRUE
)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = metadata$group,
  Sample = metadata$sample_id
)

var_explained <- round(
  100 *
    (
      pca$sdev^2 /
        sum(pca$sdev^2)
    ),
  1
)

p_pca <- ggplot2::ggplot(
  pca_df,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    color = Group,
    label = Sample
  )
) +
  ggplot2::geom_point(
    size = 3
  ) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = paste(
      "PCA Plot -",
      disease_label,
      "vs Control"
    ),
    x = paste0(
      "PC1 (",
      var_explained[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      var_explained[2],
      "%)"
    )
  )

ggplot2::ggsave(
  paste0(
    out_prefix,
    "_pca.png"
  ),
  p_pca,
  width = 6,
  height = 5,
  dpi = 300
)

cat(
  "[SUCCESS] PCA salvo.\n"
)


# ------------------------------------------------------------
# 9. Volcano plot
# ------------------------------------------------------------

p_volcano <- ggplot2::ggplot(
  res_table,
  ggplot2::aes(
    x = logFC,
    y = -log10(P.Value)
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      color = significant
    ),
    alpha = 0.6,
    size = 1.5
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "FALSE" = "grey60",
      "TRUE" = "red3"
    )
  ) +
  ggplot2::geom_vline(
    xintercept = c(
      -0.58,
      0.58
    ),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = paste(
      "Volcano Plot:",
      disease_label,
      "vs Control"
    ),
    x = "Log2 Fold Change",
    y = "-Log10 P-Value",
    color = paste0(
      "Significant\n",
      "(FDR < 0.05 and |log2FC| > 0.58)"
    )
  )

ggplot2::ggsave(
  paste0(
    out_prefix,
    "_volcano.png"
  ),
  p_volcano,
  width = 7,
  height = 6,
  dpi = 300
)

cat(
  "[SUCCESS] Volcano plot salvo.\n"
)


# ------------------------------------------------------------
# 10. Heatmap
#
# Top 50 genes ranqueados pelo menor FDR.
# Eles não são necessariamente DEGs significativos.
# ------------------------------------------------------------

top_genes <- res_table %>%
  dplyr::filter(
    !is.na(adj.P.Val)
  ) %>%
  dplyr::arrange(
    adj.P.Val,
    P.Value
  ) %>%
  dplyr::slice_head(
    n = 50
  ) %>%
  dplyr::pull(
    gene_id
  )

top_genes <- intersect(
  top_genes,
  rownames(v$E)
)

mat_heatmap <- v$E[
  top_genes,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 10.1 Usar símbolos como labels quando disponíveis
# ------------------------------------------------------------

heatmap_labels <- res_table %>%
  dplyr::filter(
    gene_id %in% top_genes
  ) %>%
  dplyr::select(
    gene_id,
    gene_symbol_display
  )

heatmap_labels <- heatmap_labels[
  match(
    rownames(mat_heatmap),
    heatmap_labels$gene_id
  ),
  ,
  drop = FALSE
]

rownames(mat_heatmap) <-
  make.unique(
    heatmap_labels$gene_symbol_display
  )


# ------------------------------------------------------------
# 10.2 Gerar heatmap
# ------------------------------------------------------------

if (nrow(mat_heatmap) > 1) {

  annotation_col <- data.frame(
    Group = metadata$group
  )

  rownames(annotation_col) <-
    metadata$sample_id

  pheatmap::pheatmap(
    mat_heatmap,
    scale = "row",
    annotation_col = annotation_col,
    show_colnames = TRUE,
    show_rownames = TRUE,
    main = paste(
      "Heatmap - Top 50 genes by FDR -",
      disease_label
    ),
    filename = paste0(
      out_prefix,
      "_heatmap.png"
    ),
    width = 8,
    height = 9
  )

  cat(
    "[SUCCESS] Heatmap salvo.\n"
  )

} else {

  warning(
    "Não foi possível gerar o heatmap: ",
    "menos de 2 genes disponíveis."
  )
}


# ------------------------------------------------------------
# 11. Finalização
# ------------------------------------------------------------

cat(
  "[SUCCESS] Análise concluída para",
  disease_label,
  "vs Control.\n"
)

