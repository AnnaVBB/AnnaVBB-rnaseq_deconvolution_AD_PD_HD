#!/usr/bin/env Rscript

# ==============================================================================
# 05_deseq2.R
#
# OBJETIVOS
# ------------------------------------------------------------------------------
# 1. Executar expressão diferencial com DESeq2 separadamente para AD, HD e PD.
# 2. Utilizar Control como grupo de referência.
# 3. Aplicar o teste padrão do DESeq2 para H0: log2FC = 0.
# 4. Classificar DEGs pelo critério principal:
#
#       padj < 0.05 e |log2FoldChange| > 0.58
#
# 5. Realizar análise de sensibilidade:
#
#       padj < 0.05 e |log2FoldChange| > 1
#
# 6. Anotar os resultados usando gene_master.tsv.
# 7. Gerar:
#       - tabela completa do DESeq2;
#       - tabela de DEGs do critério principal;
#       - MA plot;
#       - volcano plot;
#       - heatmap dos até 50 DEGs mais significativos.
#
# IMPORTANTE
# ------------------------------------------------------------------------------
# - Cada doença é analisada independentemente.
# - O modelo principal é ~ group.
# - Nenhuma covariável adicional é incluída nesta análise principal.
# - O limiar de log2FC NÃO é passado para results().
# - Independent filtering do DESeq2 permanece ativado.
# - O log2FoldChange não-shrunken é utilizado para a classificação principal.
# - VST é utilizado apenas para visualização do heatmap.
# - Nenhuma amostra é removida nesta etapa.
#
# ENTRADAS
# ------------------------------------------------------------------------------
# results/qc/counts_AD_DE_filtered.tsv
# results/qc/counts_HD_DE_filtered.tsv
# results/qc/counts_PD_DE_filtered.tsv
#
# results/metadata/metadata_AD.tsv
# results/metadata/metadata_HD.tsv
# results/metadata/metadata_PD.tsv
#
# results/harmonization/gene_master.tsv
#
# SAÍDAS
# ------------------------------------------------------------------------------
# results/deseq2/
#   AD_deseq2_results.tsv
#   AD_DEGs.tsv
#   HD_deseq2_results.tsv
#   HD_DEGs.tsv
#   PD_deseq2_results.tsv
#   PD_DEGs.tsv
#   deseq2_summary.tsv
#
# results/plots/AD/
#   MA_plot.png
#   volcano_plot.png
#   heatmap_top50_DEGs.png
#
# results/plots/HD/
#   ...
#
# results/plots/PD/
#   ...
#
# ==============================================================================


# ==============================================================================
# 1. PACOTES
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})


# ==============================================================================
# 2. PARÂMETROS
# ==============================================================================

alpha <- 0.05

# Critério principal
lfc_threshold <- 0.58

# Análise de sensibilidade
lfc_sensitivity <- 1

# Número máximo de DEGs no heatmap
n_heatmap <- 50

# Número máximo de genes rotulados no volcano
n_volcano_labels <- 10


# ==============================================================================
# 3. CORES
# ==============================================================================

# Convenção fixa para MA plot e volcano:
#
# Up      = vermelho
# Down    = azul
# Not_DEG = cinza

deg_colors <- c(
  "Down" = "#2166AC",
  "Not_DEG" = "#BDBDBD",
  "Up" = "#B2182B"
)


# Cores dos grupos utilizadas na anotação do heatmap.
# Control permanece igual nos três datasets.

group_colors <- list(
  AD = c(
    "Control" = "#4D4D4D",
    "AD" = "#D55E00"
  ),
  HD = c(
    "Control" = "#4D4D4D",
    "HD" = "#0072B2"
  ),
  PD = c(
    "Control" = "#4D4D4D",
    "PD" = "#009E73"
  )
)


# ==============================================================================
# 4. CAMINHOS
# ==============================================================================

counts_files <- c(
  AD = "results/qc/counts_AD_DE_filtered.tsv",
  HD = "results/qc/counts_HD_DE_filtered.tsv",
  PD = "results/qc/counts_PD_DE_filtered.tsv"
)


metadata_files <- c(
  AD = "results/metadata/metadata_AD.tsv",
  HD = "results/metadata/metadata_HD.tsv",
  PD = "results/metadata/metadata_PD.tsv"
)


gene_master_file <- "results/harmonization/gene_master.tsv"


output_dir <- "results/deseq2"

plots_dir <- "results/plots"


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


for (dataset in c("AD", "HD", "PD")) {

  dir.create(
    file.path(
      plots_dir,
      dataset
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 5. VERIFICAR ARQUIVOS DE ENTRADA
# ==============================================================================

input_files <- c(
  counts_files,
  metadata_files,
  gene_master_file
)


missing_files <- input_files[
  !file.exists(input_files)
]


if (length(missing_files) > 0) {

  stop(
    paste0(
      "Arquivos necessários não encontrados:\n",
      paste(
        missing_files,
        collapse = "\n"
      ),
      "\n\nO 05_deseq2.R deve ser executado somente após ",
      "a conclusão das etapas 03 e 04."
    )
  )
}


# ==============================================================================
# 6. CARREGAR E VALIDAR GENE MASTER
# ==============================================================================

gene_master <- fread(
  gene_master_file,
  colClasses = list(
    character = "gene_id"
  )
)


if (!"gene_id" %in% names(gene_master)) {

  stop(
    "gene_master.tsv não contém a coluna gene_id."
  )
}


if (anyNA(gene_master$gene_id)) {

  stop(
    "gene_master.tsv contém gene_id ausente."
  )
}


if (anyDuplicated(gene_master$gene_id)) {

  stop(
    "gene_master.tsv contém gene_id duplicado."
  )
}


# Detectar os nomes reais das colunas de anotação.
# O script tolera pequenas diferenças de nomenclatura no gene_master.

find_column <- function(candidates, available_names) {

  hit <- candidates[
    candidates %in% available_names
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


symbol_col <- find_column(
  c(
    "symbol",
    "Symbol",
    "final_symbol",
    "gene_symbol"
  ),
  names(gene_master)
)


ensembl_col <- find_column(
  c(
    "ensembl_id",
    "EnsemblGeneID",
    "ensembl_gene_id"
  ),
  names(gene_master)
)


gene_type_col <- find_column(
  c(
    "gene_type",
    "GeneType",
    "geneType"
  ),
  names(gene_master)
)


annotation_dt <- data.table(
  gene_id = as.character(
    gene_master$gene_id
  )
)


if (!is.na(symbol_col)) {

  annotation_dt[
    ,
    symbol := as.character(
      gene_master[[symbol_col]]
    )
  ]

} else {

  annotation_dt[
    ,
    symbol := NA_character_
  ]
}


if (!is.na(ensembl_col)) {

  annotation_dt[
    ,
    ensembl_id := as.character(
      gene_master[[ensembl_col]]
    )
  ]

} else {

  annotation_dt[
    ,
    ensembl_id := NA_character_
  ]
}


if (!is.na(gene_type_col)) {

  annotation_dt[
    ,
    gene_type := as.character(
      gene_master[[gene_type_col]]
    )
  ]

} else {

  annotation_dt[
    ,
    gene_type := NA_character_
  ]
}


cat("\nAnotação carregada do gene_master.tsv\n")
cat(
  "Genes no gene_master:",
  nrow(annotation_dt),
  "\n"
)
cat(
  "Genes com símbolo:",
  sum(
    !is.na(annotation_dt$symbol) &
      annotation_dt$symbol != ""
  ),
  "\n"
)


# ==============================================================================
# 7. FUNÇÃO PRINCIPAL
# ==============================================================================

run_deseq2 <- function(
  dataset_name,
  disease_level,
  counts_file,
  metadata_file
) {

  cat("\n")
  cat("============================================\n")
  cat(
    "DESeq2 -",
    dataset_name,
    "\n"
  )
  cat("============================================\n")


  # ============================================================================
  # 7.1. Ler contagens
  # ============================================================================

  counts_dt <- fread(
    counts_file,
    check.names = FALSE
  )


  if (!"GeneID" %in% names(counts_dt)) {

    stop(
      dataset_name,
      ": coluna GeneID não encontrada."
    )
  }


  gene_ids <- as.character(
    counts_dt$GeneID
  )


  if (anyNA(gene_ids)) {

    stop(
      dataset_name,
      ": existem GeneIDs ausentes."
    )
  }


  if (anyDuplicated(gene_ids)) {

    stop(
      dataset_name,
      ": existem GeneIDs duplicados."
    )
  }


  sample_columns <- setdiff(
    names(counts_dt),
    "GeneID"
  )


  count_matrix <- as.matrix(
    counts_dt[
      ,
      ..sample_columns
    ]
  )


  storage.mode(count_matrix) <- "numeric"

  rownames(count_matrix) <- gene_ids


  # ============================================================================
  # 7.2. Validar contagens
  # ============================================================================

  if (anyNA(count_matrix)) {

    stop(
      dataset_name,
      ": existem NA na matriz de contagens."
    )
  }


  if (any(count_matrix < 0)) {

    stop(
      dataset_name,
      ": existem contagens negativas."
    )
  }


  n_fractional <- sum(
    abs(
      count_matrix -
        round(count_matrix)
    ) > 1e-8
  )


  if (n_fractional > 0) {

    stop(
      dataset_name,
      ": existem ",
      n_fractional,
      " valores não inteiros."
    )
  }


  storage.mode(count_matrix) <- "integer"


  # ============================================================================
  # 7.3. Ler metadata
  # ============================================================================

  metadata <- fread(
    metadata_file,
    na.strings = c(
      "NA",
      ""
    )
  )


  required_metadata <- c(
    "sample_id",
    "group"
  )


  missing_metadata <- setdiff(
    required_metadata,
    names(metadata)
  )


  if (length(missing_metadata) > 0) {

    stop(
      dataset_name,
      ": campos ausentes no metadata: ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  }


  if (anyNA(metadata$sample_id)) {

    stop(
      dataset_name,
      ": existem sample_id ausentes."
    )
  }


  if (anyDuplicated(metadata$sample_id)) {

    stop(
      dataset_name,
      ": existem sample_id duplicados."
    )
  }


  if (anyNA(metadata$group)) {

    stop(
      dataset_name,
      ": existem grupos ausentes."
    )
  }


  # ============================================================================
  # 7.4. Validar counts <-> metadata
  # ============================================================================

  if (
    !setequal(
      colnames(count_matrix),
      metadata$sample_id
    )
  ) {

    stop(
      dataset_name,
      ": amostras de counts e metadata não correspondem."
    )
  }


  metadata <- metadata[
    match(
      colnames(count_matrix),
      sample_id
    )
  ]


  if (
    !identical(
      metadata$sample_id,
      colnames(count_matrix)
    )
  ) {

    stop(
      dataset_name,
      ": falha ao alinhar counts e metadata."
    )
  }


  # ============================================================================
  # 7.5. Definir grupos
  # ============================================================================

  metadata[
    ,
    group := factor(
      group,
      levels = c(
        "Control",
        disease_level
      )
    )
  ]


  if (anyNA(metadata$group)) {

    stop(
      dataset_name,
      ": grupo inesperado detectado."
    )
  }


  cat(
    "Genes recebidos do QC:",
    nrow(count_matrix),
    "\n"
  )


  cat(
    "Amostras:",
    ncol(count_matrix),
    "\n"
  )


  cat("\nGrupos:\n")

  print(
    table(
      metadata$group
    )
  )


  # ============================================================================
  # 7.6. Criar DESeqDataSet
  # ============================================================================

  coldata <- as.data.frame(
    metadata
  )


  rownames(coldata) <- coldata$sample_id


  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = ~ group
  )


  # ============================================================================
  # 7.7. Executar DESeq2
  # ============================================================================

  cat("\nExecutando DESeq2...\n")


  dds <- DESeq(
    dds
  )


  # ============================================================================
  # 7.8. Verificar contraste
  # ============================================================================

  cat("\nCoeficientes disponíveis:\n")

  print(
    resultsNames(dds)
  )


  # ============================================================================
  # 7.9. Extrair resultado
  #
  # O teste padrão permanece H0: log2FC = 0.
  #
  # NÃO usamos lfcThreshold aqui.
  #
  # O limiar de efeito é aplicado posteriormente apenas para classificação.
  # ============================================================================

  res <- results(
    dds,
    contrast = c(
      "group",
      disease_level,
      "Control"
    ),
    alpha = alpha,
    independentFiltering = TRUE
  )


  # ============================================================================
  # 7.10. Transformar resultado em tabela
  # ============================================================================

  res_dt <- as.data.table(
    as.data.frame(res),
    keep.rownames = "gene_id"
  )


  res_dt[
    ,
    gene_id := as.character(
      gene_id
    )
  ]


  # ============================================================================
  # 7.11. Adicionar anotação
  # ============================================================================

  res_dt <- merge(
    res_dt,
    annotation_dt,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
  )


  # Recuperar explicitamente a ordem original do resultado do DESeq2.

  result_order <- match(
    rownames(res),
    res_dt$gene_id
  )


  res_dt <- res_dt[
    result_order
  ]


  if (
    !identical(
      res_dt$gene_id,
      rownames(res)
    )
  ) {

    stop(
      dataset_name,
      ": falha ao preservar a ordem dos resultados após anotação."
    )
  }


  # Label para figuras.
  # Símbolo é preferido; GeneID é usado quando símbolo não está disponível.

  res_dt[
    ,
    plot_label := fifelse(
      !is.na(symbol) &
        symbol != "",
      symbol,
      gene_id
    )
  ]


  # ============================================================================
  # 7.12. Classificação principal
  # ============================================================================

  res_dt[
    ,
    significant_FDR := (
      !is.na(padj) &
        padj < alpha
    )
  ]


  res_dt[
    ,
    DEG := (
      !is.na(padj) &
        padj < alpha &
        !is.na(log2FoldChange) &
        abs(log2FoldChange) > lfc_threshold
    )
  ]


  res_dt[
    ,
    direction := fcase(
      DEG &
        log2FoldChange > lfc_threshold,
      "Up",

      DEG &
        log2FoldChange < -lfc_threshold,
      "Down",

      default = "Not_DEG"
    )
  ]


  res_dt[
    ,
    direction := factor(
      direction,
      levels = c(
        "Down",
        "Not_DEG",
        "Up"
      )
    )
  ]


  # ============================================================================
  # 7.13. Análise de sensibilidade
  #
  # Mesmo controle de FDR:
  #
  #   padj < 0.05
  #
  # com limiar de efeito mais restritivo:
  #
  #   |log2FC| > 1
  #
  # ============================================================================

  res_dt[
    ,
    DEG_sensitivity := (
      !is.na(padj) &
        padj < alpha &
        !is.na(log2FoldChange) &
        abs(log2FoldChange) > lfc_sensitivity
    )
  ]


  res_dt[
    ,
    direction_sensitivity := fcase(
      DEG_sensitivity &
        log2FoldChange > lfc_sensitivity,
      "Up",

      DEG_sensitivity &
        log2FoldChange < -lfc_sensitivity,
      "Down",

      default = "Not_DEG"
    )
  ]


  # ============================================================================
  # 7.14. Estatísticas
  # ============================================================================

  n_total <- nrow(
    res_dt
  )


  n_pvalue <- sum(
    !is.na(res_dt$pvalue)
  )


  n_padj <- sum(
    !is.na(res_dt$padj)
  )


  n_fdr <- sum(
    res_dt$significant_FDR
  )


  n_deg <- sum(
    res_dt$DEG
  )


  n_up <- sum(
    res_dt$direction == "Up"
  )


  n_down <- sum(
    res_dt$direction == "Down"
  )


  n_deg_sensitivity <- sum(
    res_dt$DEG_sensitivity
  )


  n_up_sensitivity <- sum(
    res_dt$direction_sensitivity == "Up"
  )


  n_down_sensitivity <- sum(
    res_dt$direction_sensitivity == "Down"
  )


  cat("\n")
  cat("--------------------------------------------\n")
  cat("RESUMO DA EXPRESSÃO DIFERENCIAL\n")
  cat("--------------------------------------------\n")


  cat(
    "Genes recebidos:",
    n_total,
    "\n"
  )


  cat(
    "Genes com p-value:",
    n_pvalue,
    "\n"
  )


  cat(
    "Genes com padj:",
    n_padj,
    "\n"
  )


  cat(
    "FDR <",
    alpha,
    ":",
    n_fdr,
    "\n"
  )


  cat("\nCRITÉRIO PRIMÁRIO\n")


  cat(
    "DEGs (padj <",
    alpha,
    "e |log2FC| >",
    lfc_threshold,
    "):",
    n_deg,
    "\n"
  )


  cat(
    "Up:",
    n_up,
    "\n"
  )


  cat(
    "Down:",
    n_down,
    "\n"
  )


  cat("\nANÁLISE DE SENSIBILIDADE\n")


  cat(
    "DEGs (padj <",
    alpha,
    "e |log2FC| >",
    lfc_sensitivity,
    "):",
    n_deg_sensitivity,
    "\n"
  )


  cat(
    "Up:",
    n_up_sensitivity,
    "\n"
  )


  cat(
    "Down:",
    n_down_sensitivity,
    "\n"
  )


  # ============================================================================
  # 7.15. Ordenar resultados para exportação
  # ============================================================================

  res_dt[
    ,
    abs_log2FoldChange := abs(
      log2FoldChange
    )
  ]


  setorder(
    res_dt,
    padj,
    pvalue,
    -abs_log2FoldChange,
    na.last = TRUE
  )


  deg_dt <- res_dt[
    DEG == TRUE
  ]


  # ============================================================================
  # 7.16. Salvar tabelas
  # ============================================================================

  full_output <- file.path(
    output_dir,
    paste0(
      dataset_name,
      "_deseq2_results.tsv"
    )
  )


  deg_output <- file.path(
    output_dir,
    paste0(
      dataset_name,
      "_DEGs.tsv"
    )
  )


  fwrite(
    res_dt,
    full_output,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )


  fwrite(
    deg_dt,
    deg_output,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )


  cat("\nTabelas salvas:\n")
  cat(
    " -",
    full_output,
    "\n"
  )
  cat(
    " -",
    deg_output,
    "\n"
  )


  # ============================================================================
  # 7.17. MA PLOT
  #
  # x = média normalizada das contagens (baseMean)
  # y = log2FoldChange
  #
  # Up      = vermelho
  # Down    = azul
  # Not_DEG = cinza
  # ============================================================================

  ma_dt <- res_dt[
    !is.na(baseMean) &
      baseMean > 0 &
      !is.na(log2FoldChange)
  ]


  ma_plot <- ggplot(
    ma_dt,
    aes(
      x = baseMean,
      y = log2FoldChange,
      color = direction
    )
  ) +
    geom_point(
      alpha = 0.55,
      size = 1.3
    ) +
    geom_hline(
      yintercept = c(
        -lfc_threshold,
        lfc_threshold
      ),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    scale_x_log10() +
    scale_color_manual(
      values = deg_colors,
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "MA plot - ",
        dataset_name
      ),
      subtitle = paste0(
        "Primary DEG criterion: FDR < ",
        alpha,
        " and |log2FC| > ",
        lfc_threshold
      ),
      x = "Mean normalized expression (baseMean)",
      y = "log2 fold change",
      color = "Classification"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(
        face = "bold"
      )
    )


  ma_output <- file.path(
    plots_dir,
    dataset_name,
    "MA_plot.png"
  )


  ggsave(
    filename = ma_output,
    plot = ma_plot,
    width = 8,
    height = 6,
    dpi = 300
  )


  cat(
    "MA plot:",
    ma_output,
    "\n"
  )


  # ============================================================================
  # 7.18. VOLCANO PLOT
  #
  # x = log2FoldChange
  # y = -log10(padj)
  #
  # Linhas:
  #   FDR = 0.05
  #   log2FC = +/- 0.58
  #
  # Os 10 DEGs com menor padj são rotulados.
  # ============================================================================

  volcano_dt <- copy(
    res_dt[
      !is.na(log2FoldChange) &
        !is.na(padj)
    ]
  )


  # Evitar -log10(0) = Inf exclusivamente para visualização.
  #
  # Os valores originais de padj NÃO são modificados na tabela do DESeq2.

  positive_padj <- volcano_dt[
    padj > 0,
    padj
  ]


  if (length(positive_padj) > 0) {

    plot_padj_floor <- min(
      positive_padj
    ) / 10

  } else {

    plot_padj_floor <- .Machine$double.xmin
  }


  plot_padj_floor <- max(
    plot_padj_floor,
    .Machine$double.xmin
  )


  volcano_dt[
    ,
    padj_plot := pmax(
      padj,
      plot_padj_floor
    )
  ]


  volcano_dt[
    ,
    neg_log10_padj := -log10(
      padj_plot
    )
  ]


  # Top genes para labels:
  # primeiro menor padj, depois maior |log2FC|.

  volcano_labels <- volcano_dt[
    DEG == TRUE
  ]


  setorder(
    volcano_labels,
    padj,
    -abs_log2FoldChange,
    na.last = TRUE
  )


  volcano_labels <- head(
    volcano_labels,
    n_volcano_labels
  )


  volcano_plot <- ggplot(
    volcano_dt,
    aes(
      x = log2FoldChange,
      y = neg_log10_padj,
      color = direction
    )
  ) +
    geom_point(
      alpha = 0.60,
      size = 1.4
    ) +
    geom_vline(
      xintercept = c(
        -lfc_threshold,
        lfc_threshold
      ),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    geom_hline(
      yintercept = -log10(
        alpha
      ),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    scale_color_manual(
      values = deg_colors,
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "Volcano plot - ",
        dataset_name
      ),
      subtitle = paste0(
        "Primary DEG criterion: FDR < ",
        alpha,
        " and |log2FC| > ",
        lfc_threshold
      ),
      x = "log2 fold change",
      y = expression(
        -log[10]("adjusted p-value")
      ),
      color = "Classification"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(
        face = "bold"
      )
    )


  if (nrow(volcano_labels) > 0) {

    volcano_plot <- volcano_plot +
      geom_text_repel(
        data = volcano_labels,
        aes(
          label = plot_label
        ),
        size = 3.2,
        box.padding = 0.45,
        point.padding = 0.25,
        min.segment.length = 0,
        max.overlaps = Inf,
        show.legend = FALSE,
        seed = 123
      )
  }


  volcano_output <- file.path(
    plots_dir,
    dataset_name,
    "volcano_plot.png"
  )


  ggsave(
    filename = volcano_output,
    plot = volcano_plot,
    width = 8.5,
    height = 6.5,
    dpi = 300
  )


  cat(
    "Volcano plot:",
    volcano_output,
    "\n"
  )


  # ============================================================================
  # 7.19. HEATMAP DOS TOP 50 DEGs
  #
  # Seleção:
  #   1. somente DEGs do critério principal;
  #   2. menor padj;
  #   3. maior |log2FC| como segundo critério;
  #   4. máximo de 50 genes.
  #
  # Expressão:
  #   VST -> Z-score por gene.
  #
  # IMPORTANTE:
  #   VST é usado somente para visualização.
  # ============================================================================

  if (nrow(deg_dt) == 0) {

    cat(
      "Heatmap:",
      dataset_name,
      "não possui DEGs pelo critério principal; figura não gerada.\n"
    )

  } else {

    top_degs <- copy(
      deg_dt
    )


    setorder(
      top_degs,
      padj,
      -abs_log2FoldChange,
      na.last = TRUE
    )


    top_degs <- head(
      top_degs,
      n_heatmap
    )


    top_gene_ids <- top_degs$gene_id


    cat(
      "Heatmap:",
      length(top_gene_ids),
      "DEGs selecionados.\n"
    )


    # VST sobre o DESeqDataSet já estimado.
    #
    # blind = FALSE:
    # neste ponto o objetivo é visualização dos DEGs após o modelo,
    # preservando a estrutura do design experimental.

    vsd <- varianceStabilizingTransformation(
      dds,
      blind = FALSE
    )


    vst_matrix <- assay(
      vsd
    )


    heatmap_matrix <- vst_matrix[
      top_gene_ids,
      ,
      drop = FALSE
    ]


    # Z-score por gene.

    heatmap_z <- t(
      scale(
        t(
          heatmap_matrix
        )
      )
    )


    # Remover genes com variância zero caso exista algum.
    # Isso evita NA após scale().

    valid_heatmap_rows <- apply(
      heatmap_z,
      1,
      function(x) {
        all(
          is.finite(x)
        )
      }
    )


    heatmap_z <- heatmap_z[
      valid_heatmap_rows,
      ,
      drop = FALSE
    ]


    if (nrow(heatmap_z) == 0) {

      cat(
        "Heatmap:",
        dataset_name,
        "não possui genes com variância suficiente para gerar a figura.\n"
      )

    } else {

      # Labels dos genes.

      heatmap_annotation <- top_degs[
        match(
          rownames(heatmap_z),
          gene_id
        )
      ]


      gene_labels <- heatmap_annotation$plot_label


      # Evitar labels duplicados quando dois GeneIDs compartilham símbolo.

      duplicated_labels <- duplicated(
        gene_labels
      ) |
        duplicated(
          gene_labels,
          fromLast = TRUE
        )


      gene_labels[
        duplicated_labels
      ] <- paste0(
        gene_labels[
          duplicated_labels
        ],
        " [",
        heatmap_annotation$gene_id[
          duplicated_labels
        ],
        "]"
      )


      rownames(
        heatmap_z
      ) <- gene_labels


      # Anotação das amostras.

      annotation_col <- data.frame(
        Group = metadata$group,
        row.names = metadata$sample_id
      )


      annotation_colors <- list(
        Group = group_colors[[dataset_name]]
      )


      heatmap_output <- file.path(
        plots_dir,
        dataset_name,
        "heatmap_top50_DEGs.png"
      )


      # Ajustar altura conforme o número de genes.

      heatmap_height <- max(
        7,
        min(
          14,
          4 +
            0.18 *
              nrow(heatmap_z)
        )
      )


      png(
        filename = heatmap_output,
        width = 2600,
        height = round(
          300 *
            heatmap_height
        ),
        res = 300
      )


      pheatmap(
        heatmap_z,
        annotation_col = annotation_col,
        annotation_colors = annotation_colors,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        show_colnames = FALSE,
        fontsize_row = if (
          nrow(heatmap_z) > 35
        ) {
          6
        } else {
          8
        },
        main = paste0(
          "Top ",
          nrow(heatmap_z),
          " DEGs - ",
          dataset_name
        ),
        border_color = NA
      )


      dev.off()


      cat(
        "Heatmap:",
        heatmap_output,
        "\n"
      )
    }
  }


  # ============================================================================
  # 7.20. Retornar resumo
  # ============================================================================

  return(
    data.table(
      dataset = dataset_name,
      disease = disease_level,
      n_samples = ncol(count_matrix),
      n_genes_input = n_total,
      n_pvalue = n_pvalue,
      n_padj = n_padj,
      n_FDR_005 = n_fdr,

      # Critério principal:
      # padj < 0.05 e |log2FC| > 0.58

      n_DEG_primary = n_deg,
      n_up_primary = n_up,
      n_down_primary = n_down,

      # Sensibilidade:
      # padj < 0.05 e |log2FC| > 1

      n_DEG_sensitivity = n_deg_sensitivity,
      n_up_sensitivity = n_up_sensitivity,
      n_down_sensitivity = n_down_sensitivity
    )
  )
}


# ==============================================================================
# 8. EXECUTAR AD
# ==============================================================================

summary_AD <- run_deseq2(
  dataset_name = "AD",
  disease_level = "AD",
  counts_file = counts_files[["AD"]],
  metadata_file = metadata_files[["AD"]]
)


# ==============================================================================
# 9. EXECUTAR HD
# ==============================================================================

summary_HD <- run_deseq2(
  dataset_name = "HD",
  disease_level = "HD",
  counts_file = counts_files[["HD"]],
  metadata_file = metadata_files[["HD"]]
)


# ==============================================================================
# 10. EXECUTAR PD
# ==============================================================================

summary_PD <- run_deseq2(
  dataset_name = "PD",
  disease_level = "PD",
  counts_file = counts_files[["PD"]],
  metadata_file = metadata_files[["PD"]]
)


# ==============================================================================
# 11. RESUMO FINAL
# ==============================================================================

summary_all <- rbindlist(
  list(
    summary_AD,
    summary_HD,
    summary_PD
  )
)


summary_output <- file.path(
  output_dir,
  "deseq2_summary.tsv"
)


fwrite(
  summary_all,
  summary_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


cat("\n")
cat("============================================\n")
cat("RESUMO FINAL DESeq2\n")
cat("============================================\n")


print(
  summary_all
)


cat(
  "\nResumo salvo em:",
  summary_output,
  "\n"
)


# ==============================================================================
# 12. FINALIZAÇÃO
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("05_deseq2.R CONCLUÍDO\n")
cat("============================================\n")


cat(
  "Critério primário de DEG: padj <",
  alpha,
  "e |log2FoldChange| >",
  lfc_threshold,
  "\n"
)


cat(
  "Sensibilidade: padj <",
  alpha,
  "e |log2FoldChange| >",
  lfc_sensitivity,
  "\n"
)


cat(
  "Modelo: ~ group\n"
)


cat(
  "Referência: Control\n"
)


cat(
  "Independent filtering do DESeq2: ativado\n"
)


cat(
  "Nenhuma covariável adicional foi incluída no modelo principal.\n"
)


cat(
  "Cada doença foi analisada independentemente.\n"
)


cat(
  "MA plot e volcano utilizam a classificação principal de DEGs.\n"
)


cat(
  "Volcano: Up = vermelho; Down = azul; Not_DEG = cinza.\n"
)


cat(
  "Heatmap: até 50 DEGs com menor padj; desempate por maior |log2FC|.\n"
)


cat(
  "VST no heatmap foi utilizada somente para visualização.\n"
)