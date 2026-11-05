#!/usr/bin/env Rscript

# ==============================================================================
# 06_pathway_enrichment.R
#
# OBJETIVOS
# ------------------------------------------------------------------------------
# 1. Realizar enriquecimento funcional separadamente para AD, HD e PD.
#
# 2. ORA (Over-Representation Analysis):
#       - GO Biological Process
#       - KEGG
#
#    Utilizando os DEGs do critério principal:
#
#       FDR < 0.05 e |log2FoldChange| > 0.58
#
#    Conjuntos:
#       - All
#       - Up
#       - Down
#
# 3. GSEA:
#       - GO Biological Process
#       - KEGG
#
#    Utilizando todos os genes com estatística de Wald válida.
#
# 4. Utilizar NCBI GeneID / ENTREZID.
#
# 5. Na ORA, utilizar como universo apenas os genes efetivamente
#    testados pelo DESeq2 em cada dataset.
#
# INTERPRETAÇÃO DO GSEA
# ------------------------------------------------------------------------------
# NES > 0:
#   conjunto enriquecido no topo do ranking, isto é, genes associados
#   à expressão aumentada na doença.
#
# NES < 0:
#   conjunto enriquecido na parte inferior do ranking, isto é, genes
#   associados à expressão reduzida na doença / maior expressão em controles.
#
# ENTRADAS
# ------------------------------------------------------------------------------
# results/deseq2/AD_deseq2_results.tsv
# results/deseq2/HD_deseq2_results.tsv
# results/deseq2/PD_deseq2_results.tsv
#
# SAÍDAS
# ------------------------------------------------------------------------------
# results/enrichment/
#   ORA_GO_BP.tsv
#   ORA_KEGG.tsv
#   GSEA_GO_BP.tsv
#   GSEA_KEGG.tsv
#   enrichment_summary.tsv
#
# results/plots/AD/
# results/plots/HD/
# results/plots/PD/
#
# ==============================================================================


# ==============================================================================
# 1. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table",
  "ggplot2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "fgsea"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (length(missing_packages) > 0) {

  stop(
    paste0(
      "Pacotes necessários não instalados:\n",
      paste(
        missing_packages,
        collapse = "\n"
      )
    )
  )
}


suppressPackageStartupMessages({

  library(data.table)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)

})


# ==============================================================================
# 2. PARÂMETROS
# ==============================================================================

alpha <- 0.05

min_gs_size <- 10
max_gs_size <- 500

# Número máximo de termos ORA mostrados por conjunto
# All / Up / Down.
n_ora_terms_per_set <- 7

# Número máximo de termos mostrados em cada gráfico GSEA.
n_gsea_terms <- 15


# ==============================================================================
# 3. CORES
# ==============================================================================

direction_colors <- c(
  "Down" = "#2166AC",
  "Up" = "#B2182B"
)


ora_colors <- c(
  "All" = "#666666",
  "Down" = "#2166AC",
  "Up" = "#B2182B"
)


# ==============================================================================
# 4. DATASETS E CAMINHOS
# ==============================================================================

datasets <- c(
  "AD",
  "HD",
  "PD"
)


deseq_files <- setNames(
  file.path(
    "results/deseq2",
    paste0(
      datasets,
      "_deseq2_results.tsv"
    )
  ),
  datasets
)


output_dir <- "results/enrichment"

plots_dir <- "results/plots"


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


for (dataset_name in datasets) {

  dir.create(
    file.path(
      plots_dir,
      dataset_name
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 5. VALIDAR ARQUIVOS DE ENTRADA
# ==============================================================================

missing_files <- deseq_files[
  !file.exists(deseq_files)
]


if (length(missing_files) > 0) {

  stop(
    paste0(
      "Arquivos necessários não encontrados:\n",
      paste(
        missing_files,
        collapse = "\n"
      ),
      "\n\nO 06_pathway_enrichment.R deve ser executado após ",
      "o 05_deseq2.R."
    )
  )
}


# ==============================================================================
# 6. FUNÇÕES AUXILIARES
# ==============================================================================


# ------------------------------------------------------------------------------
# 6.1. Quebrar descrições longas para os gráficos
# ------------------------------------------------------------------------------

wrap_label <- function(
  x,
  width = 50
) {

  vapply(
    x,
    function(label) {

      paste(
        strwrap(
          label,
          width = width
        ),
        collapse = "\n"
      )

    },
    character(1)
  )
}


# ------------------------------------------------------------------------------
# 6.2. Converter GeneIDs separados por "/" em símbolos
# ------------------------------------------------------------------------------

ids_to_symbols <- function(
  id_string,
  symbol_map
) {

  if (
    is.na(id_string) ||
    id_string == ""
  ) {

    return(
      NA_character_
    )
  }


  ids <- strsplit(
    id_string,
    "/",
    fixed = TRUE
  )[[1]]


  symbols <- unname(
    symbol_map[
      ids
    ]
  )


  missing_symbol <- (
    is.na(symbols) |
      symbols == ""
  )


  symbols[
    missing_symbol
  ] <- ids[
    missing_symbol
  ]


  paste(
    unique(symbols),
    collapse = "/"
  )
}


# ------------------------------------------------------------------------------
# 6.3. ORA GO Biological Process
# ------------------------------------------------------------------------------

run_ora_go <- function(
  gene_ids,
  background_ids
) {

  if (length(gene_ids) == 0) {

    return(NULL)
  }


  tryCatch(

    enrichGO(
      gene = gene_ids,
      universe = background_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",

      # Recuperar todos os termos possíveis.
      # A classificação por FDR será realizada depois.
      pvalueCutoff = 1,
      qvalueCutoff = 1,

      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,

      readable = FALSE
    ),

    error = function(e) {

      message(
        "Erro no enrichGO: ",
        conditionMessage(e)
      )

      NULL
    }
  )
}


# ------------------------------------------------------------------------------
# 6.4. ORA KEGG
# ------------------------------------------------------------------------------

run_ora_kegg <- function(
  gene_ids,
  background_ids
) {

  if (length(gene_ids) == 0) {

    return(NULL)
  }


  tryCatch(

    enrichKEGG(
      gene = gene_ids,
      universe = background_ids,

      organism = "hsa",
      keyType = "ncbi-geneid",

      pAdjustMethod = "BH",

      pvalueCutoff = 1,
      qvalueCutoff = 1,

      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,

      use_internal_data = FALSE
    ),

    error = function(e) {

      message(
        "Erro no enrichKEGG: ",
        conditionMessage(e)
      )

      NULL
    }
  )
}


# ------------------------------------------------------------------------------
# 6.5. GSEA GO Biological Process
# ------------------------------------------------------------------------------

run_gsea_go <- function(
  gene_list
) {

  tryCatch(

    gseGO(
      geneList = gene_list,

      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",

      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,

      pvalueCutoff = 1,
      pAdjustMethod = "BH",

      eps = 0,

      verbose = FALSE
    ),

    error = function(e) {

      message(
        "Erro no gseGO: ",
        conditionMessage(e)
      )

      NULL
    }
  )
}


# ------------------------------------------------------------------------------
# 6.6. GSEA KEGG
# ------------------------------------------------------------------------------

run_gsea_kegg <- function(
  gene_list
) {

  tryCatch(

    gseKEGG(
      geneList = gene_list,

      organism = "hsa",
      keyType = "ncbi-geneid",

      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,

      pvalueCutoff = 1,
      pAdjustMethod = "BH",

      eps = 0,

      verbose = FALSE
    ),

    error = function(e) {

      message(
        "Erro no gseKEGG: ",
        conditionMessage(e)
      )

      NULL
    }
  )
}


# ------------------------------------------------------------------------------
# 6.7. Formatar resultado ORA
# ------------------------------------------------------------------------------

format_ora_result <- function(
  enrichment_object,
  dataset_name,
  database_name,
  gene_set_name,
  symbol_map
) {

  if (is.null(enrichment_object)) {

    return(NULL)
  }


  enrichment_df <- as.data.frame(
    enrichment_object
  )


  if (nrow(enrichment_df) == 0) {

    return(NULL)
  }


  enrichment_dt <- as.data.table(
    enrichment_df
  )


  enrichment_dt[
    ,
    dataset := dataset_name
  ]


  enrichment_dt[
    ,
    database := database_name
  ]


  enrichment_dt[
    ,
    gene_set := gene_set_name
  ]


  enrichment_dt[
    ,
    significant := (
      !is.na(p.adjust) &
        p.adjust < alpha
    )
  ]


  if ("geneID" %in% names(enrichment_dt)) {

    enrichment_dt[
      ,
      gene_symbols := vapply(
        geneID,
        ids_to_symbols,
        character(1),
        symbol_map = symbol_map
      )
    ]
  }


  preferred_columns <- c(
    "dataset",
    "database",
    "gene_set",
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "pvalue",
    "p.adjust",
    "qvalue",
    "significant",
    "Count",
    "geneID",
    "gene_symbols"
  )


  preferred_columns <- preferred_columns[
    preferred_columns %in%
      names(enrichment_dt)
  ]


  setcolorder(
    enrichment_dt,
    c(
      preferred_columns,
      setdiff(
        names(enrichment_dt),
        preferred_columns
      )
    )
  )


  return(
    enrichment_dt
  )
}


# ------------------------------------------------------------------------------
# 6.8. Formatar resultado GSEA
# ------------------------------------------------------------------------------

format_gsea_result <- function(
  enrichment_object,
  dataset_name,
  database_name,
  symbol_map
) {

  if (is.null(enrichment_object)) {

    return(NULL)
  }


  enrichment_df <- as.data.frame(
    enrichment_object
  )


  if (nrow(enrichment_df) == 0) {

    return(NULL)
  }


  enrichment_dt <- as.data.table(
    enrichment_df
  )


  enrichment_dt[
    ,
    dataset := dataset_name
  ]


  enrichment_dt[
    ,
    database := database_name
  ]


  enrichment_dt[
    ,
    direction := fcase(
      NES > 0,
      "Up",

      NES < 0,
      "Down",

      default = NA_character_
    )
  ]


  enrichment_dt[
    ,
    significant := (
      !is.na(p.adjust) &
        p.adjust < alpha
    )
  ]


  if (
    "core_enrichment" %in%
      names(enrichment_dt)
  ) {

    enrichment_dt[
      ,
      core_symbols := vapply(
        core_enrichment,
        ids_to_symbols,
        character(1),
        symbol_map = symbol_map
      )
    ]
  }


  preferred_columns <- c(
    "dataset",
    "database",
    "ID",
    "Description",
    "setSize",
    "enrichmentScore",
    "NES",
    "direction",
    "pvalue",
    "p.adjust",
    "qvalue",
    "significant",
    "rank",
    "leading_edge",
    "core_enrichment",
    "core_symbols"
  )


  preferred_columns <- preferred_columns[
    preferred_columns %in%
      names(enrichment_dt)
  ]


  setcolorder(
    enrichment_dt,
    c(
      preferred_columns,
      setdiff(
        names(enrichment_dt),
        preferred_columns
      )
    )
  )


  return(
    enrichment_dt
  )
}


# ------------------------------------------------------------------------------
# 6.9. Criar gráfico ORA
# ------------------------------------------------------------------------------

plot_ora <- function(
  ora_dt,
  dataset_name,
  database_name
) {

  if (
    is.null(ora_dt) ||
    nrow(ora_dt) == 0
  ) {

    return(
      invisible(NULL)
    )
  }


  plot_dt <- ora_dt[
    dataset == dataset_name &
      database == database_name &
      significant == TRUE
  ]


  if (nrow(plot_dt) == 0) {

    message(
      dataset_name,
      " - ",
      database_name,
      ": nenhum termo ORA com FDR < ",
      alpha,
      "."
    )

    return(
      invisible(NULL)
    )
  }


  # Selecionar os termos mais significativos
  # dentro de cada conjunto All / Up / Down.

  plot_dt <- plot_dt[
    order(
      p.adjust,
      -Count
    ),
    head(
      .SD,
      n_ora_terms_per_set
    ),
    by = gene_set
  ]


  plot_dt[
    ,
    p_adjust_plot := pmax(
      p.adjust,
      .Machine$double.xmin
    )
  ]


  plot_dt[
    ,
    neg_log10_fdr := -log10(
      p_adjust_plot
    )
  ]


  plot_dt[
    ,
    plot_label := paste0(
      wrap_label(
        Description,
        width = 45
      ),
      " [",
      gene_set,
      "]"
    )
  ]


  plot_dt <- plot_dt[
    order(
      neg_log10_fdr
    )
  ]


  plot_dt[
    ,
    plot_label := factor(
      plot_label,
      levels = unique(
        plot_label
      )
    )
  ]


  p <- ggplot(
    plot_dt,
    aes(
      x = neg_log10_fdr,
      y = plot_label,
      color = gene_set,
      size = Count
    )
  ) +
    geom_point(
      alpha = 0.85
    ) +
    scale_color_manual(
      values = ora_colors,
      drop = FALSE
    ) +
    labs(
      title = paste0(
        database_name,
        " ORA - ",
        dataset_name
      ),
      subtitle = paste0(
        "Significant terms: BH-FDR < ",
        alpha
      ),
      x = expression(
        -log[10]("adjusted p-value")
      ),
      y = NULL,
      color = "Gene set",
      size = "Genes"
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 8
      ),
      legend.position = "right"
    )


  output_file <- file.path(
    plots_dir,
    dataset_name,
    paste0(
      "enrichment_ORA_",
      database_name,
      ".png"
    )
  )


  ggsave(
    filename = output_file,
    plot = p,
    width = 10,
    height = 7,
    dpi = 300
  )


  message(
    "Gráfico salvo: ",
    output_file
  )


  invisible(NULL)
}


# ------------------------------------------------------------------------------
# 6.10. Criar gráfico GSEA
# ------------------------------------------------------------------------------

plot_gsea <- function(
  gsea_dt,
  dataset_name,
  database_name
) {

  if (
    is.null(gsea_dt) ||
    nrow(gsea_dt) == 0
  ) {

    return(
      invisible(NULL)
    )
  }


  plot_dt <- gsea_dt[
    dataset == dataset_name &
      database == database_name &
      significant == TRUE
  ]


  if (nrow(plot_dt) == 0) {

    message(
      dataset_name,
      " - ",
      database_name,
      ": nenhum termo GSEA com FDR < ",
      alpha,
      "."
    )

    return(
      invisible(NULL)
    )
  }


  plot_dt[
    ,
    abs_NES := abs(
      NES
    )
  ]


  setorder(
    plot_dt,
    p.adjust,
    -abs_NES
  )


  plot_dt <- head(
    plot_dt,
    n_gsea_terms
  )


  plot_dt[
    ,
    plot_label := wrap_label(
      Description,
      width = 48
    )
  ]


  plot_dt <- plot_dt[
    order(
      NES
    )
  ]


  plot_dt[
    ,
    plot_label := factor(
      plot_label,
      levels = unique(
        plot_label
      )
    )
  ]


  p <- ggplot(
    plot_dt,
    aes(
      x = NES,
      y = plot_label,
      color = direction,
      size = setSize
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_point(
      alpha = 0.9
    ) +
    scale_color_manual(
      values = direction_colors,
      drop = FALSE
    ) +
    labs(
      title = paste0(
        database_name,
        " GSEA - ",
        dataset_name
      ),
      subtitle = paste0(
        "Significant terms: BH-FDR < ",
        alpha
      ),
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      color = "Direction",
      size = "Gene set size"
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 8
      ),
      legend.position = "right"
    )


  output_file <- file.path(
    plots_dir,
    dataset_name,
    paste0(
      "enrichment_GSEA_",
      database_name,
      ".png"
    )
  )


  ggsave(
    filename = output_file,
    plot = p,
    width = 10,
    height = 7,
    dpi = 300
  )


  message(
    "Gráfico salvo: ",
    output_file
  )


  invisible(NULL)
}


# ==============================================================================
# 7. LISTAS PARA ACUMULAR RESULTADOS
# ==============================================================================

ora_go_results <- list()
ora_kegg_results <- list()

gsea_go_results <- list()
gsea_kegg_results <- list()

summary_results <- list()


# ==============================================================================
# 8. EXECUTAR DATASET POR DATASET
# ==============================================================================

for (dataset_name in datasets) {

  cat("\n")
  cat("============================================================\n")

  cat(
    "ENRIQUECIMENTO FUNCIONAL - ",
    dataset_name,
    "\n",
    sep = ""
  )

  cat("============================================================\n")


  # ============================================================================
  # 8.1. Ler resultado do DESeq2
  # ============================================================================

  res_dt <- fread(
    deseq_files[[dataset_name]]
  )


  required_columns <- c(
    "gene_id",
    "symbol",
    "stat",
    "pvalue",
    "padj",
    "log2FoldChange",
    "DEG",
    "direction"
  )


  missing_columns <- setdiff(
    required_columns,
    names(res_dt)
  )


  if (length(missing_columns) > 0) {

    stop(
      dataset_name,
      ": colunas ausentes no resultado DESeq2: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }


  res_dt[
    ,
    gene_id := as.character(
      gene_id
    )
  ]


  res_dt[
    ,
    symbol := as.character(
      symbol
    )
  ]


  if (anyDuplicated(res_dt$gene_id)) {

    stop(
      dataset_name,
      ": existem GeneIDs duplicados no resultado DESeq2."
    )
  }


  # ============================================================================
  # 8.2. Criar GeneID -> Symbol
  # ============================================================================

  symbol_map <- setNames(
    res_dt$symbol,
    res_dt$gene_id
  )


  # ============================================================================
  # 8.3. Definir universo da ORA
  #
  # Genes que efetivamente receberam teste estatístico.
  # ============================================================================

  background_ids <- unique(
    res_dt[
      !is.na(pvalue),
      gene_id
    ]
  )


  # ============================================================================
  # 8.4. Definir DEGs principais
  # ============================================================================

  all_deg <- unique(
    res_dt[
      DEG == TRUE,
      gene_id
    ]
  )


  up_deg <- unique(
    res_dt[
      DEG == TRUE &
        direction == "Up",
      gene_id
    ]
  )


  down_deg <- unique(
    res_dt[
      DEG == TRUE &
        direction == "Down",
      gene_id
    ]
  )


  cat(
    "Genes no background:",
    length(background_ids),
    "\n"
  )


  cat(
    "DEGs principais:",
    length(all_deg),
    "\n"
  )


  cat(
    "  Up:",
    length(up_deg),
    "\n"
  )


  cat(
    "  Down:",
    length(down_deg),
    "\n"
  )


  # ============================================================================
  # 8.5. ORA
  # ============================================================================

  deg_sets <- list(
    All = all_deg,
    Up = up_deg,
    Down = down_deg
  )


  for (gene_set_name in names(deg_sets)) {

    gene_ids <- deg_sets[[gene_set_name]]


    cat(
      "\nORA - ",
      gene_set_name,
      ": ",
      length(gene_ids),
      " genes\n",
      sep = ""
    )


    # --------------------------------------------------------------------------
    # Sem DEGs
    # --------------------------------------------------------------------------

    if (length(gene_ids) == 0) {

      cat(
        "Sem genes. ORA ignorada.\n"
      )


      summary_results[[length(summary_results) + 1]] <- NULL

      summary_results[[length(summary_results) + 1]] <- data.table(
        dataset = dataset_name,
        method = "ORA",
        database = "GO_BP",
        gene_set = gene_set_name,
        n_input = 0L,
        n_background = length(
          background_ids
        ),
        n_terms_tested = 0L,
        n_terms_FDR005 = 0L,
        n_positive = NA_integer_,
        n_negative = NA_integer_,
        status = "skipped_no_genes"
      )


      summary_results[[length(summary_results) + 1]] <- data.table(
        dataset = dataset_name,
        method = "ORA",
        database = "KEGG",
        gene_set = gene_set_name,
        n_input = 0L,
        n_background = length(
          background_ids
        ),
        n_terms_tested = 0L,
        n_terms_FDR005 = 0L,
        n_positive = NA_integer_,
        n_negative = NA_integer_,
        status = "skipped_no_genes"
      )


      next
    }


    # --------------------------------------------------------------------------
    # GO Biological Process
    # --------------------------------------------------------------------------

    go_object <- run_ora_go(
      gene_ids = gene_ids,
      background_ids = background_ids
    )


    go_dt <- format_ora_result(
      enrichment_object = go_object,
      dataset_name = dataset_name,
      database_name = "GO_BP",
      gene_set_name = gene_set_name,
      symbol_map = symbol_map
    )


    if (!is.null(go_dt)) {

      ora_go_results[[length(ora_go_results) + 1]] <- NULL

      ora_go_results[[length(ora_go_results) + 1]] <- go_dt


      n_go_tested <- nrow(
        go_dt
      )


      n_go_sig <- sum(
        go_dt$significant,
        na.rm = TRUE
      )

    } else {

      n_go_tested <- 0L
      n_go_sig <- 0L
    }


    cat(
      "GO BP significativos:",
      n_go_sig,
      "\n"
    )


    summary_results[[length(summary_results) + 1]] <- data.table(
      dataset = dataset_name,
      method = "ORA",
      database = "GO_BP",
      gene_set = gene_set_name,
      n_input = length(
        gene_ids
      ),
      n_background = length(
        background_ids
      ),
      n_terms_tested = n_go_tested,
      n_terms_FDR005 = n_go_sig,
      n_positive = NA_integer_,
      n_negative = NA_integer_,
      status = "completed"
    )


    # --------------------------------------------------------------------------
    # KEGG
    # --------------------------------------------------------------------------

    kegg_object <- run_ora_kegg(
      gene_ids = gene_ids,
      background_ids = background_ids
    )


    kegg_dt <- format_ora_result(
      enrichment_object = kegg_object,
      dataset_name = dataset_name,
      database_name = "KEGG",
      gene_set_name = gene_set_name,
      symbol_map = symbol_map
    )


    if (!is.null(kegg_dt)) {

      ora_kegg_results[[length(ora_kegg_results) + 1]] <- kegg_dt


      n_kegg_tested <- nrow(
        kegg_dt
      )


      n_kegg_sig <- sum(
        kegg_dt$significant,
        na.rm = TRUE
      )

    } else {

      n_kegg_tested <- 0L
      n_kegg_sig <- 0L
    }


    cat(
      "KEGG significativos:",
      n_kegg_sig,
      "\n"
    )


    summary_results[[length(summary_results) + 1]] <- data.table(
      dataset = dataset_name,
      method = "ORA",
      database = "KEGG",
      gene_set = gene_set_name,
      n_input = length(
        gene_ids
      ),
      n_background = length(
        background_ids
      ),
      n_terms_tested = n_kegg_tested,
      n_terms_FDR005 = n_kegg_sig,
      n_positive = NA_integer_,
      n_negative = NA_integer_,
      status = "completed"
    )
  }


  # ============================================================================
  # 8.6. Preparar ranking GSEA
  #
  # Ranking = estatística de Wald do DESeq2.
  # ============================================================================

  rank_dt <- res_dt[
    !is.na(stat) &
      is.finite(stat)
  ]


  if (nrow(rank_dt) == 0) {

    stop(
      dataset_name,
      ": nenhum gene possui estatística de Wald válida."
    )
  }


  if (anyDuplicated(rank_dt$gene_id)) {

    stop(
      dataset_name,
      ": GeneIDs duplicados no ranking GSEA."
    )
  }


  gene_list <- rank_dt$stat


  names(
    gene_list
  ) <- rank_dt$gene_id


  gene_list <- sort(
    gene_list,
    decreasing = TRUE
  )


  cat(
    "\nGenes no ranking GSEA:",
    length(gene_list),
    "\n"
  )


  cat(
    "Maior estatística:",
    max(gene_list),
    "\n"
  )


  cat(
    "Menor estatística:",
    min(gene_list),
    "\n"
  )


  # ============================================================================
  # 8.7. GSEA GO Biological Process
  # ============================================================================

  cat(
    "\nExecutando GSEA GO BP...\n"
  )


  gsea_go_object <- run_gsea_go(
    gene_list
  )


  gsea_go_dt <- format_gsea_result(
    enrichment_object = gsea_go_object,
    dataset_name = dataset_name,
    database_name = "GO_BP",
    symbol_map = symbol_map
  )


  if (!is.null(gsea_go_dt)) {

    gsea_go_results[[length(gsea_go_results) + 1]] <- gsea_go_dt


    n_gsea_go_tested <- nrow(
      gsea_go_dt
    )


    n_gsea_go_sig <- sum(
      gsea_go_dt$significant,
      na.rm = TRUE
    )


    n_gsea_go_positive <- sum(
      gsea_go_dt$significant &
        gsea_go_dt$direction == "Up",
      na.rm = TRUE
    )


    n_gsea_go_negative <- sum(
      gsea_go_dt$significant &
        gsea_go_dt$direction == "Down",
      na.rm = TRUE
    )

  } else {

    n_gsea_go_tested <- 0L
    n_gsea_go_sig <- 0L
    n_gsea_go_positive <- 0L
    n_gsea_go_negative <- 0L
  }


  cat(
    "GO BP GSEA significativos:",
    n_gsea_go_sig,
    "\n"
  )


  cat(
    "  NES positivo:",
    n_gsea_go_positive,
    "\n"
  )


  cat(
    "  NES negativo:",
    n_gsea_go_negative,
    "\n"
  )


  summary_results[[length(summary_results) + 1]] <- data.table(
    dataset = dataset_name,
    method = "GSEA",
    database = "GO_BP",
    gene_set = "Ranked",
    n_input = length(
      gene_list
    ),
    n_background = NA_integer_,
    n_terms_tested = n_gsea_go_tested,
    n_terms_FDR005 = n_gsea_go_sig,
    n_positive = n_gsea_go_positive,
    n_negative = n_gsea_go_negative,
    status = "completed"
  )


  # ============================================================================
  # 8.8. GSEA KEGG
  # ============================================================================

  cat(
    "\nExecutando GSEA KEGG...\n"
  )


  gsea_kegg_object <- run_gsea_kegg(
    gene_list
  )


  gsea_kegg_dt <- format_gsea_result(
    enrichment_object = gsea_kegg_object,
    dataset_name = dataset_name,
    database_name = "KEGG",
    symbol_map = symbol_map
  )


  if (!is.null(gsea_kegg_dt)) {

    gsea_kegg_results[[length(gsea_kegg_results) + 1]] <- gsea_kegg_dt


    n_gsea_kegg_tested <- nrow(
      gsea_kegg_dt
    )


    n_gsea_kegg_sig <- sum(
      gsea_kegg_dt$significant,
      na.rm = TRUE
    )


    n_gsea_kegg_positive <- sum(
      gsea_kegg_dt$significant &
        gsea_kegg_dt$direction == "Up",
      na.rm = TRUE
    )


    n_gsea_kegg_negative <- sum(
      gsea_kegg_dt$significant &
        gsea_kegg_dt$direction == "Down",
      na.rm = TRUE
    )

  } else {

    n_gsea_kegg_tested <- 0L
    n_gsea_kegg_sig <- 0L
    n_gsea_kegg_positive <- 0L
    n_gsea_kegg_negative <- 0L
  }


  cat(
    "KEGG GSEA significativos:",
    n_gsea_kegg_sig,
    "\n"
  )


  cat(
    "  NES positivo:",
    n_gsea_kegg_positive,
    "\n"
  )


  cat(
    "  NES negativo:",
    n_gsea_kegg_negative,
    "\n"
  )


  summary_results[[length(summary_results) + 1]] <- data.table(
    dataset = dataset_name,
    method = "GSEA",
    database = "KEGG",
    gene_set = "Ranked",
    n_input = length(
      gene_list
    ),
    n_background = NA_integer_,
    n_terms_tested = n_gsea_kegg_tested,
    n_terms_FDR005 = n_gsea_kegg_sig,
    n_positive = n_gsea_kegg_positive,
    n_negative = n_gsea_kegg_negative,
    status = "completed"
  )
}


# ==============================================================================
# 9. CONSOLIDAR RESULTADOS
# ==============================================================================

ora_go_all <- if (length(ora_go_results) > 0) {

  rbindlist(
    ora_go_results,
    fill = TRUE
  )

} else {

  data.table()
}


ora_kegg_all <- if (length(ora_kegg_results) > 0) {

  rbindlist(
    ora_kegg_results,
    fill = TRUE
  )

} else {

  data.table()
}


gsea_go_all <- if (length(gsea_go_results) > 0) {

  rbindlist(
    gsea_go_results,
    fill = TRUE
  )

} else {

  data.table()
}


gsea_kegg_all <- if (length(gsea_kegg_results) > 0) {

  rbindlist(
    gsea_kegg_results,
    fill = TRUE
  )

} else {

  data.table()
}


summary_all <- rbindlist(
  summary_results,
  fill = TRUE
)


# ==============================================================================
# 10. SALVAR TABELAS
# ==============================================================================

fwrite(
  ora_go_all,
  file.path(
    output_dir,
    "ORA_GO_BP.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  ora_kegg_all,
  file.path(
    output_dir,
    "ORA_KEGG.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  gsea_go_all,
  file.path(
    output_dir,
    "GSEA_GO_BP.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  gsea_kegg_all,
  file.path(
    output_dir,
    "GSEA_KEGG.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  summary_all,
  file.path(
    output_dir,
    "enrichment_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 11. GERAR GRÁFICOS
# ==============================================================================

for (dataset_name in datasets) {

  # ---------------------------------------------------------------------------
  # ORA
  # ---------------------------------------------------------------------------

  plot_ora(
    ora_dt = ora_go_all,
    dataset_name = dataset_name,
    database_name = "GO_BP"
  )


  plot_ora(
    ora_dt = ora_kegg_all,
    dataset_name = dataset_name,
    database_name = "KEGG"
  )


  # ---------------------------------------------------------------------------
  # GSEA
  # ---------------------------------------------------------------------------

  plot_gsea(
    gsea_dt = gsea_go_all,
    dataset_name = dataset_name,
    database_name = "GO_BP"
  )


  plot_gsea(
    gsea_dt = gsea_kegg_all,
    dataset_name = dataset_name,
    database_name = "KEGG"
  )
}


# ==============================================================================
# 12. RESUMO FINAL
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("RESUMO FINAL - ENRIQUECIMENTO FUNCIONAL\n")
cat("============================================================\n")


print(
  summary_all
)


cat("\n")
cat(
  "Resultados salvos em: ",
  output_dir,
  "\n",
  sep = ""
)


cat("\n")
cat("============================================================\n")
cat("06_pathway_enrichment.R CONCLUÍDO\n")
cat("============================================================\n")


cat(
  "ORA: GO Biological Process + KEGG\n"
)


cat(
  "GSEA: GO Biological Process + KEGG\n"
)


cat(
  "FDR para enriquecimento: ",
  alpha,
  "\n",
  sep = ""
)


cat(
  "ORA background: genes testados pelo DESeq2 de cada dataset.\n"
)


cat(
  "ORA realizada separadamente para All, Up e Down.\n"
)


cat(
  "GSEA ranqueado pela estatística de Wald do DESeq2.\n"
)


cat(
  "NES positivo = enriquecimento na direção da doença.\n"
)


cat(
  "NES negativo = enriquecimento na direção do controle.\n"
)


cat(
  "AD sem DEGs: ORA ignorada; GSEA permanece disponível.\n"
)