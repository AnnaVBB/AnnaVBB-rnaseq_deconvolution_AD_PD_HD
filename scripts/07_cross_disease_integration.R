#!/usr/bin/env Rscript

# ==============================================================================
# 07_cross_disease_integration.R
#
# INTEGRAÇÃO TRANSDIAGNÓSTICA
#
# Doenças:
#   AD = Alzheimer
#   HD = Huntington
#   PD = Parkinson
#
# OBJETIVOS
# ------------------------------------------------------------------------------
# Nível gênico:
#
# 1. Construir universo comum testável.
#
# 2. Comparar DEGs pelo critério principal:
#
#       padj < 0.05 e |log2FoldChange| > 0.58
#
# 3. Comparar DEGs pela análise de sensibilidade:
#
#       padj < 0.05 e |log2FoldChange| > 1
#
# 4. Identificar:
#       - genes exclusivos;
#       - genes compartilhados;
#       - interseção HD + PD;
#       - concordância de direção;
#       - genes com mesmo sinal nas três doenças.
#
# 5. Calcular correlações de Spearman dos log2FC.
#
# Nível funcional:
#
# 6. Identificar termos GSEA compartilhados entre doenças.
#
# 7. Avaliar concordância de NES.
#
# 8. Comparar ORA entre HD e PD.
#
# Visualizações:
#
# 9. Distribuição dos padrões de DEGs.
#
# 10. Diagramas de Venn:
#       - critério principal;
#       - análise de sensibilidade.
#
# 11. Scatter HD vs PD.
#
# 12. Heatmaps de NES compartilhados.
#
# IMPORTANTE
# ------------------------------------------------------------------------------
# - NCBI GeneID é a chave primária.
#
# - Símbolo gênico é usado apenas como anotação.
#
# - Ausência de DEG não significa ausência de efeito.
#
# - AD possui zero DEGs pelo critério principal.
#
# - A integração transdiagnóstica não depende apenas da
#   interseção estrita de listas de DEGs.
#
# - NES positivo/negativo representa direção do enriquecimento,
#   não ativação/inibição fisiológica direta.
#
# - Os diagramas de Venn são esquemáticos.
#   O tamanho dos círculos NÃO é proporcional ao número de genes.
#
# ENTRADAS
# ------------------------------------------------------------------------------
# results/deseq2/AD_deseq2_results.tsv
# results/deseq2/HD_deseq2_results.tsv
# results/deseq2/PD_deseq2_results.tsv
#
# results/enrichment/GSEA_GO_BP.tsv
# results/enrichment/GSEA_KEGG.tsv
# results/enrichment/ORA_GO_BP.tsv
# results/enrichment/ORA_KEGG.tsv
#
# SAÍDAS PRINCIPAIS
# ------------------------------------------------------------------------------
# results/integration/
#   gene_integration.tsv
#   primary_pattern_summary.tsv
#   sensitivity_pattern_summary.tsv
#   shared_primary_HD_PD.tsv
#   directionally_concordant_AD_HD_PD.tsv
#   log2FC_spearman_correlations.tsv
#   shared_GSEA_GO_BP.tsv
#   shared_GSEA_KEGG.tsv
#   shared_ORA_HD_PD_GO_BP.tsv
#   shared_ORA_HD_PD_KEGG.tsv
#   integration_summary.tsv
#
# results/plots/integration/
#   primary_deg_patterns.png
#   venn_primary_AD_HD_PD.png
#   venn_sensitivity_AD_HD_PD.png
#   log2FC_HD_vs_PD.png
#   shared_GSEA_GO_BP.png
#   shared_GSEA_KEGG.png
#
# ==============================================================================


# ==============================================================================
# 1. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table",
  "ggplot2"
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

})


# ==============================================================================
# 2. PARÂMETROS
# ==============================================================================

alpha <- 0.05


datasets <- c(
  "AD",
  "HD",
  "PD"
)


# ==============================================================================
# 3. CORES
# ==============================================================================

disease_colors <- c(

  "AD" = "#7B2CBF",

  "HD" = "#277DA1",

  "PD" = "#E75480"

)


pattern_colors <- c(

  "AD only" = "#7B2CBF",

  "HD only" = "#277DA1",

  "PD only" = "#E75480",

  "AD+HD" = "#D10404",

  "AD+PD" = "#F8A304",

  "HD+PD" = "#18C900",

  "AD+HD+PD" = "#00BFC4"

)


# ==============================================================================
# 4. CAMINHOS
# ==============================================================================

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


gsea_files <- c(

  "GO_BP" = "results/enrichment/GSEA_GO_BP.tsv",

  "KEGG" = "results/enrichment/GSEA_KEGG.tsv"

)


ora_files <- c(

  "GO_BP" = "results/enrichment/ORA_GO_BP.tsv",

  "KEGG" = "results/enrichment/ORA_KEGG.tsv"

)


output_dir <- "results/integration"

plot_dir <- "results/plots/integration"


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 5. VALIDAR ENTRADAS
# ==============================================================================

input_files <- c(
  deseq_files,
  gsea_files,
  ora_files
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
      )
    )
  )
}


# ==============================================================================
# 6. FUNÇÕES AUXILIARES
# ==============================================================================


# ------------------------------------------------------------------------------
# 6.1. Converter TRUE/FALSE de forma segura
# ------------------------------------------------------------------------------

as_logical_safe <- function(x) {

  if (is.logical(x)) {

    return(x)
  }


  x_chr <- toupper(
    as.character(x)
  )


  out <- rep(
    NA,
    length(x_chr)
  )


  out[
    x_chr %in% c(
      "TRUE",
      "T",
      "1"
    )
  ] <- TRUE


  out[
    x_chr %in% c(
      "FALSE",
      "F",
      "0"
    )
  ] <- FALSE


  bad_idx <- (
    is.na(out) &
      !is.na(x)
  )


  if (any(bad_idx)) {

    stop(
      "Valores lógicos inesperados: ",
      paste(
        unique(
          x_chr[
            bad_idx
          ]
        ),
        collapse = ", "
      )
    )
  }


  return(out)
}


# ------------------------------------------------------------------------------
# 6.2. Direção do efeito
# ------------------------------------------------------------------------------

effect_direction <- function(x) {

  fcase(

    is.na(x),
    NA_character_,

    x > 0,
    "Up",

    x < 0,
    "Down",

    default = "Zero"
  )
}


# ------------------------------------------------------------------------------
# 6.3. Primeiro símbolo válido disponível
# ------------------------------------------------------------------------------

first_nonmissing <- function(...) {

  values <- list(...)


  out <- rep(
    NA_character_,
    length(
      values[[1]]
    )
  )


  for (v in values) {

    v <- as.character(v)


    use <- (
      is.na(out) &
        !is.na(v) &
        v != ""
    )


    out[
      use
    ] <- v[
      use
    ]
  }


  return(out)
}


# ------------------------------------------------------------------------------
# 6.4. Ler resultado do DESeq2
# ------------------------------------------------------------------------------

read_deseq_result <- function(dataset_name) {

  dt <- fread(
    deseq_files[[dataset_name]]
  )


  required <- c(

    "gene_id",
    "symbol",

    "log2FoldChange",
    "stat",
    "pvalue",
    "padj",

    "DEG",
    "direction",

    "DEG_sensitivity",
    "direction_sensitivity"
  )


  missing <- setdiff(
    required,
    names(dt)
  )


  if (length(missing) > 0) {

    stop(
      dataset_name,
      ": colunas ausentes no resultado do DESeq2: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }


  dt[
    ,
    gene_id := as.character(
      gene_id
    )
  ]


  dt[
    ,
    symbol := as.character(
      symbol
    )
  ]


  dt[
    ,
    DEG := as_logical_safe(
      DEG
    )
  ]


  dt[
    ,
    DEG_sensitivity := as_logical_safe(
      DEG_sensitivity
    )
  ]


  if (anyNA(dt$gene_id)) {

    stop(
      dataset_name,
      ": existem GeneIDs ausentes."
    )
  }


  if (anyDuplicated(dt$gene_id)) {

    stop(
      dataset_name,
      ": existem GeneIDs duplicados."
    )
  }


  return(dt)
}


# ------------------------------------------------------------------------------
# 6.5. Criar resumo completo de padrões
# ------------------------------------------------------------------------------

pattern_levels <- c(

  "AD only",
  "HD only",
  "PD only",

  "AD+HD",
  "AD+PD",
  "HD+PD",

  "AD+HD+PD",

  "None"
)


make_pattern_summary <- function(
  pattern_vector,
  column_name
) {

  counts <- data.table(
    pattern = pattern_vector
  )[
    ,
    .(
      n_genes = .N
    ),
    by = pattern
  ]


  out <- data.table(
    pattern = pattern_levels
  )


  out <- merge(

    out,
    counts,

    by = "pattern",

    all.x = TRUE,

    sort = FALSE
  )


  out[
    is.na(n_genes),
    n_genes := 0L
  ]


  setnames(
    out,
    "pattern",
    column_name
  )


  return(out)
}


# ------------------------------------------------------------------------------
# 6.6. Criar diagrama de Venn customizado
#
# IMPORTANTE:
# os círculos têm o mesmo tamanho.
# Portanto, a área NÃO representa a quantidade de genes.
# ------------------------------------------------------------------------------

create_venn_plot <- function(
  pattern_dt,
  pattern_column,
  title,
  subtitle
) {

  dt <- copy(
    pattern_dt
  )


  setnames(
    dt,
    pattern_column,
    "pattern"
  )


  get_n <- function(pattern_name) {

    value <- dt[
      pattern == pattern_name,
      n_genes
    ]


    if (length(value) == 0) {

      return(0L)
    }


    return(
      as.integer(
        value[1]
      )
    )
  }


  theta <- seq(
    0,
    2 * pi,
    length.out = 600
  )


  circle_data <- rbindlist(

    list(

      data.table(
        x = -0.85 + 1.45 * cos(theta),
        y =  0.55 + 1.45 * sin(theta),
        disease = "AD"
      ),

      data.table(
        x =  0.85 + 1.45 * cos(theta),
        y =  0.55 + 1.45 * sin(theta),
        disease = "HD"
      ),

      data.table(
        x =  0.00 + 1.45 * cos(theta),
        y = -0.85 + 1.45 * sin(theta),
        disease = "PD"
      )

    )
  )


  counts <- data.table(

    pattern = c(
      "AD only",
      "HD only",
      "PD only",
      "AD+HD",
      "AD+PD",
      "HD+PD",
      "AD+HD+PD"
    ),

    n = c(
      get_n("AD only"),
      get_n("HD only"),
      get_n("PD only"),
      get_n("AD+HD"),
      get_n("AD+PD"),
      get_n("HD+PD"),
      get_n("AD+HD+PD")
    ),

    x = c(
      -1.25,
       1.25,
       0.00,
       0.00,
      -0.75,
       0.75,
       0.00
    ),

    y = c(
       0.85,
       0.85,
      -1.35,
       1.25,
      -0.35,
      -0.35,
       0.25
    )
  )


  p <- ggplot() +

    geom_polygon(
      data = circle_data,
      aes(
        x = x,
        y = y,
        group = disease,
        fill = disease
      ),
      alpha = 0.28,
      color = NA
    ) +

    geom_path(
      data = circle_data,
      aes(
        x = x,
        y = y,
        group = disease,
        color = disease
      ),
      linewidth = 1.4
    ) +

    geom_text(
      data = counts,
      aes(
        x = x,
        y = y,
        label = n
      ),
      size = 7
    ) +

    annotate(
      "text",
      x = -1.75,
      y = 2.15,
      label = "AD",
      color = unname(
        disease_colors["AD"]
      ),
      size = 9,
      fontface = "bold"
    ) +

    annotate(
      "text",
      x = 1.75,
      y = 2.15,
      label = "HD",
      color = unname(
        disease_colors["HD"]
      ),
      size = 9,
      fontface = "bold"
    ) +

    annotate(
      "text",
      x = 0,
      y = -2.45,
      label = "PD",
      color = unname(
        disease_colors["PD"]
      ),
      size = 9,
      fontface = "bold"
    ) +

    scale_fill_manual(
      values = disease_colors
    ) +

    scale_color_manual(
      values = disease_colors
    ) +

    coord_fixed(
      xlim = c(
        -2.5,
        2.5
      ),
      ylim = c(
        -2.7,
        2.5
      ),
      clip = "off"
    ) +

    labs(
      title = title,
      subtitle = subtitle
    ) +

    theme_void() +

    theme(

      legend.position = "none",

      plot.title = element_text(
        size = 22,
        face = "bold"
      ),

      plot.subtitle = element_text(
        size = 14
      ),

      plot.margin = margin(
        20,
        20,
        20,
        20
      )
    )


  return(p)
}


# ==============================================================================
# 7. CARREGAR RESULTADOS
# ==============================================================================

results_list <- lapply(
  datasets,
  read_deseq_result
)


names(
  results_list
) <- datasets


# ==============================================================================
# 8. UNIVERSO TESTÁVEL POR DATASET
#
# Gene testável:
# possui p-value e log2FoldChange válidos.
# ==============================================================================

testable_ids <- lapply(

  results_list,

  function(dt) {

    unique(

      dt[
        !is.na(pvalue) &
          !is.na(log2FoldChange),
        gene_id
      ]

    )
  }
)


# ==============================================================================
# 9. UNIVERSO COMUM TESTÁVEL
# ==============================================================================

common_universe <- Reduce(
  intersect,
  testable_ids
)


cat("\n")
cat("============================================================\n")
cat("INTEGRAÇÃO TRANSDIAGNÓSTICA\n")
cat("============================================================\n")


for (dataset_name in datasets) {

  cat(
    dataset_name,
    " genes testáveis: ",
    length(
      testable_ids[[dataset_name]]
    ),
    "\n",
    sep = ""
  )
}


cat(
  "Universo comum testável: ",
  length(common_universe),
  "\n",
  sep = ""
)


# ==============================================================================
# 10. PREPARAR MATRIZES PARA INTEGRAÇÃO
# ==============================================================================

prepare_for_merge <- function(
  dt,
  prefix
) {

  x <- dt[

    gene_id %in% common_universe,

    .(

      gene_id,
      symbol,

      log2FoldChange,
      stat,
      pvalue,
      padj,

      DEG,
      direction,

      DEG_sensitivity,
      direction_sensitivity
    )
  ]


  columns_to_rename <- setdiff(
    names(x),
    "gene_id"
  )


  setnames(

    x,

    old = columns_to_rename,

    new = paste0(
      columns_to_rename,
      "_",
      prefix
    )
  )


  return(x)
}


ad_merge <- prepare_for_merge(
  results_list[["AD"]],
  "AD"
)


hd_merge <- prepare_for_merge(
  results_list[["HD"]],
  "HD"
)


pd_merge <- prepare_for_merge(
  results_list[["PD"]],
  "PD"
)


# ==============================================================================
# 11. MATRIZ INTEGRADA
# ==============================================================================

integrated <- Reduce(

  function(x, y) {

    merge(
      x,
      y,
      by = "gene_id",
      all = FALSE,
      sort = FALSE
    )
  },

  list(
    ad_merge,
    hd_merge,
    pd_merge
  )
)


if (
  nrow(integrated) !=
    length(common_universe)
) {

  stop(
    "Matriz integrada possui ",
    nrow(integrated),
    " genes, mas o universo comum possui ",
    length(common_universe),
    "."
  )
}


# ==============================================================================
# 12. SÍMBOLO CONSENSO
# ==============================================================================

integrated[
  ,
  symbol := first_nonmissing(
    symbol_AD,
    symbol_HD,
    symbol_PD
  )
]


# ==============================================================================
# 13. PADRÃO DO CRITÉRIO PRINCIPAL
# ==============================================================================

integrated[
  ,
  n_DEG_primary :=
    as.integer(DEG_AD) +
    as.integer(DEG_HD) +
    as.integer(DEG_PD)
]


integrated[
  ,
  primary_pattern := fcase(

    DEG_AD &
      DEG_HD &
      DEG_PD,

    "AD+HD+PD",


    DEG_AD &
      DEG_HD,

    "AD+HD",


    DEG_AD &
      DEG_PD,

    "AD+PD",


    DEG_HD &
      DEG_PD,

    "HD+PD",


    DEG_AD,

    "AD only",


    DEG_HD,

    "HD only",


    DEG_PD,

    "PD only",


    default = "None"
  )
]


# ==============================================================================
# 14. PADRÃO DA ANÁLISE DE SENSIBILIDADE
#
# FDR < 0.05 e |log2FC| > 1
# ==============================================================================

integrated[
  ,
  n_DEG_sensitivity :=
    as.integer(
      DEG_sensitivity_AD
    ) +
    as.integer(
      DEG_sensitivity_HD
    ) +
    as.integer(
      DEG_sensitivity_PD
    )
]


integrated[
  ,
  sensitivity_pattern := fcase(

    DEG_sensitivity_AD &
      DEG_sensitivity_HD &
      DEG_sensitivity_PD,

    "AD+HD+PD",


    DEG_sensitivity_AD &
      DEG_sensitivity_HD,

    "AD+HD",


    DEG_sensitivity_AD &
      DEG_sensitivity_PD,

    "AD+PD",


    DEG_sensitivity_HD &
      DEG_sensitivity_PD,

    "HD+PD",


    DEG_sensitivity_AD,

    "AD only",


    DEG_sensitivity_HD,

    "HD only",


    DEG_sensitivity_PD,

    "PD only",


    default = "None"
  )
]


# ==============================================================================
# 15. DIREÇÃO DO EFEITO
# ==============================================================================

integrated[
  ,
  effect_AD := effect_direction(
    log2FoldChange_AD
  )
]


integrated[
  ,
  effect_HD := effect_direction(
    log2FoldChange_HD
  )
]


integrated[
  ,
  effect_PD := effect_direction(
    log2FoldChange_PD
  )
]


integrated[
  ,
  directional_pattern := paste(
    effect_AD,
    effect_HD,
    effect_PD,
    sep = "/"
  )
]


# ==============================================================================
# 16. CONCORDÂNCIA DIRECIONAL NAS TRÊS DOENÇAS
# ==============================================================================

integrated[
  ,
  concordant_all_three := (

    effect_AD == effect_HD &

      effect_HD == effect_PD &

      effect_AD != "Zero" &

      !is.na(effect_AD)
  )
]


# ==============================================================================
# 17. RESUMOS DOS PADRÕES
# ==============================================================================

primary_pattern_summary <- make_pattern_summary(
  integrated$primary_pattern,
  "primary_pattern"
)


sensitivity_pattern_summary <- make_pattern_summary(
  integrated$sensitivity_pattern,
  "sensitivity_pattern"
)


# ==============================================================================
# 18. DEGs COMPARTILHADOS HD + PD
# ==============================================================================

shared_hd_pd <- integrated[
  DEG_HD == TRUE &
    DEG_PD == TRUE
]


shared_hd_pd[
  ,
  hd_pd_same_direction := (

    effect_HD == effect_PD &

      effect_HD != "Zero"
  )
]


setorder(
  shared_hd_pd,
  padj_HD,
  padj_PD
)


shared_hd_pd_same_direction <- shared_hd_pd[
  hd_pd_same_direction == TRUE
]


shared_hd_pd_opposite_direction <- shared_hd_pd[
  hd_pd_same_direction == FALSE
]


# ==============================================================================
# 19. GENES COM MESMO SINAL NAS TRÊS DOENÇAS
# ==============================================================================

concordant_three <- integrated[
  concordant_all_three == TRUE
]


concordant_three[
  ,
  total_abs_log2FC :=
    abs(log2FoldChange_AD) +
    abs(log2FoldChange_HD) +
    abs(log2FoldChange_PD)
]


setorder(
  concordant_three,
  -n_DEG_primary,
  -total_abs_log2FC
)


# ==============================================================================
# 20. CORRELAÇÃO ENTRE log2FC
# ==============================================================================

correlation_summary <- rbindlist(

  list(

    data.table(

      pair = "AD-HD",

      spearman_rho = cor(
        integrated$log2FoldChange_AD,
        integrated$log2FoldChange_HD,
        method = "spearman",
        use = "complete.obs"
      )
    ),


    data.table(

      pair = "AD-PD",

      spearman_rho = cor(
        integrated$log2FoldChange_AD,
        integrated$log2FoldChange_PD,
        method = "spearman",
        use = "complete.obs"
      )
    ),


    data.table(

      pair = "HD-PD",

      spearman_rho = cor(
        integrated$log2FoldChange_HD,
        integrated$log2FoldChange_PD,
        method = "spearman",
        use = "complete.obs"
      )
    )

  )
)


# ==============================================================================
# 21. INTEGRAÇÃO GSEA
# ==============================================================================

integrate_gsea <- function(
  file_path,
  database_name
) {

  dt <- fread(
    file_path
  )


  required <- c(

    "dataset",
    "ID",
    "Description",

    "NES",
    "p.adjust",

    "direction",
    "significant"
  )


  missing <- setdiff(
    required,
    names(dt)
  )


  if (length(missing) > 0) {

    stop(
      database_name,
      " GSEA: colunas ausentes: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }


  dt[
    ,
    significant := as_logical_safe(
      significant
    )
  ]


  sig <- dt[
    significant == TRUE &
      dataset %in% datasets
  ]


  if (nrow(sig) == 0) {

    return(
      list(
        long = data.table(),
        wide = data.table()
      )
    )
  }


  shared_ids <- sig[
    ,
    .(

      n_sig_diseases = uniqueN(
        dataset
      ),

      datasets_shared = paste(
        sort(
          unique(dataset)
        ),
        collapse = "+"
      ),

      same_direction = (
        uniqueN(direction) == 1
      ),

      min_FDR = min(
        p.adjust,
        na.rm = TRUE
      ),

      max_abs_NES = max(
        abs(NES),
        na.rm = TRUE
      )
    ),
    by = .(
      ID,
      Description
    )
  ][
    n_sig_diseases >= 2
  ]


  if (nrow(shared_ids) == 0) {

    return(
      list(
        long = data.table(),
        wide = data.table()
      )
    )
  }


  shared_long <- merge(

    sig,

    shared_ids[
      ,
      .(

        ID,
        Description,

        n_sig_diseases,
        datasets_shared,

        same_direction,

        min_FDR,
        max_abs_NES
      )
    ],

    by = c(
      "ID",
      "Description"
    ),

    all = FALSE
  )


  shared_wide <- dcast(

    shared_long,

    ID +
      Description +
      n_sig_diseases +
      datasets_shared +
      same_direction +
      min_FDR +
      max_abs_NES ~

      dataset,

    value.var = c(
      "NES",
      "p.adjust",
      "direction"
    )
  )


  setorder(
    shared_wide,
    -n_sig_diseases,
    min_FDR,
    -max_abs_NES
  )


  return(
    list(
      long = shared_long,
      wide = shared_wide
    )
  )
}


gsea_go <- integrate_gsea(
  gsea_files[["GO_BP"]],
  "GO_BP"
)


gsea_kegg <- integrate_gsea(
  gsea_files[["KEGG"]],
  "KEGG"
)


# ==============================================================================
# 22. INTEGRAÇÃO ORA HD + PD
# ==============================================================================

integrate_ora_hd_pd <- function(
  file_path,
  database_name
) {

  dt <- fread(
    file_path
  )


  if (nrow(dt) == 0) {

    return(
      data.table()
    )
  }


  required <- c(

    "dataset",
    "gene_set",

    "ID",
    "Description",

    "p.adjust",

    "significant"
  )


  missing <- setdiff(
    required,
    names(dt)
  )


  if (length(missing) > 0) {

    stop(
      database_name,
      " ORA: colunas ausentes: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }


  dt[
    ,
    significant := as_logical_safe(
      significant
    )
  ]


  sig <- dt[
    significant == TRUE &
      dataset %in% c(
        "HD",
        "PD"
      )
  ]


  if (nrow(sig) == 0) {

    return(
      data.table()
    )
  }


  shared <- sig[
    ,
    .(

      n_datasets = uniqueN(
        dataset
      ),

      datasets_shared = paste(
        sort(
          unique(dataset)
        ),
        collapse = "+"
      ),

      best_FDR = min(
        p.adjust,
        na.rm = TRUE
      )
    ),
    by = .(
      gene_set,
      ID,
      Description
    )
  ][
    n_datasets == 2
  ]


  setorder(
    shared,
    gene_set,
    best_FDR
  )


  return(shared)
}


ora_go_shared <- integrate_ora_hd_pd(
  ora_files[["GO_BP"]],
  "GO_BP"
)


ora_kegg_shared <- integrate_ora_hd_pd(
  ora_files[["KEGG"]],
  "KEGG"
)


# ==============================================================================
# 23. SALVAR TABELAS
# ==============================================================================

fwrite(
  integrated,
  file.path(
    output_dir,
    "gene_integration.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  primary_pattern_summary,
  file.path(
    output_dir,
    "primary_pattern_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  sensitivity_pattern_summary,
  file.path(
    output_dir,
    "sensitivity_pattern_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  shared_hd_pd,
  file.path(
    output_dir,
    "shared_primary_HD_PD.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  concordant_three,
  file.path(
    output_dir,
    "directionally_concordant_AD_HD_PD.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  correlation_summary,
  file.path(
    output_dir,
    "log2FC_spearman_correlations.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  gsea_go$wide,
  file.path(
    output_dir,
    "shared_GSEA_GO_BP.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  gsea_kegg$wide,
  file.path(
    output_dir,
    "shared_GSEA_KEGG.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  ora_go_shared,
  file.path(
    output_dir,
    "shared_ORA_HD_PD_GO_BP.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  ora_kegg_shared,
  file.path(
    output_dir,
    "shared_ORA_HD_PD_KEGG.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 24. GRÁFICO DOS PADRÕES DE DEGs
# ==============================================================================

primary_plot_dt <- primary_pattern_summary[
  primary_pattern != "None"
]


primary_plot_dt[
  ,
  primary_pattern := factor(
    primary_pattern,
    levels = pattern_levels[
      pattern_levels != "None"
    ]
  )
]


p_primary <- ggplot(

  primary_plot_dt,

  aes(
    x = primary_pattern,
    y = n_genes,
    fill = primary_pattern
  )

) +

  geom_col(
    width = 0.75
  ) +

  geom_text(
    aes(
      label = n_genes
    ),
    vjust = -0.35,
    size = 4
  ) +

  scale_fill_manual(
    values = pattern_colors,
    drop = FALSE
  ) +

  labs(

    title = "Distribuição dos DEGs entre AD, HD e PD",

    subtitle = paste0(
      "Critério principal; universo comum testável = ",
      format(
        length(common_universe),
        big.mark = "."
      ),
      " genes"
    ),

    x = NULL,

    y = "Número de genes"
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(

    legend.position = "none",

    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),

    plot.title = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    plot_dir,
    "primary_deg_patterns.png"
  ),
  p_primary,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ==============================================================================
# 25. VENN - CRITÉRIO PRINCIPAL
#
# FDR < 0.05 e |log2FC| > 0.58
# ==============================================================================

p_venn_primary <- create_venn_plot(

  pattern_dt = primary_pattern_summary,

  pattern_column = "primary_pattern",

  title = "Interseção dos DEGs entre AD, HD e PD",

  subtitle = paste0(

    "FDR < 0,05 e |log2FC| > 0,58; universo comum testável = ",

    format(
      length(common_universe),
      big.mark = "."
    ),

    " genes"
  )
)


ggsave(

  file.path(
    plot_dir,
    "venn_primary_AD_HD_PD.png"
  ),

  p_venn_primary,

  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)


# ==============================================================================
# 26. VENN - ANÁLISE DE SENSIBILIDADE
#
# FDR < 0.05 e |log2FC| > 1
# ==============================================================================

p_venn_sensitivity <- create_venn_plot(

  pattern_dt = sensitivity_pattern_summary,

  pattern_column = "sensitivity_pattern",

  title = "Interseção dos DEGs — análise de sensibilidade",

  subtitle = paste0(

    "FDR < 0,05 e |log2FC| > 1; universo comum testável = ",

    format(
      length(common_universe),
      big.mark = "."
    ),

    " genes"
  )
)


ggsave(

  file.path(
    plot_dir,
    "venn_sensitivity_AD_HD_PD.png"
  ),

  p_venn_sensitivity,

  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)


# ==============================================================================
# 27. SCATTER HD vs PD
# ==============================================================================

scatter_dt <- copy(
  integrated
)


scatter_dt[
  ,
  scatter_class := fcase(

    DEG_HD &
      DEG_PD,

    "HD+PD",


    DEG_HD,

    "HD",


    DEG_PD,

    "PD",


    default = "Other"
  )
]


scatter_colors <- c(

  "Other" = "#D0D0D0",

  "HD" = "#277DA1",

  "PD" = "#E75480",

  "HD+PD" = "#DF9203"
)


p_hd_pd <- ggplot(

  scatter_dt,

  aes(
    x = log2FoldChange_HD,
    y = log2FoldChange_PD,
    color = scatter_class
  )

) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +

  geom_point(
    alpha = 0.55,
    size = 1.3
  ) +

  scale_color_manual(
    values = scatter_colors
  ) +

  labs(

    title = "Concordância de efeito - HD vs PD",

    subtitle = paste0(
      "Universo comum testável = ",
      format(
        length(common_universe),
        big.mark = "."
      ),
      " genes"
    ),

    x = "log2FoldChange - HD",

    y = "log2FoldChange - PD",

    color = "Classificação"
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


ggsave(

  file.path(
    plot_dir,
    "log2FC_HD_vs_PD.png"
  ),

  p_hd_pd,

  width = 8,
  height = 6.5,
  dpi = 300,
  bg = "white"
)


# ==============================================================================
# 28. HEATMAP DE NES DOS TERMOS GSEA COMPARTILHADOS
# ==============================================================================

plot_shared_gsea <- function(
  gsea_result,
  database_name,
  top_n = 25
) {

  if (nrow(gsea_result$long) == 0) {

    message(
      database_name,
      ": nenhum termo GSEA significativo ",
      "em pelo menos duas doenças."
    )


    return(
      invisible(NULL)
    )
  }


  ranking <- unique(

    gsea_result$long[
      ,
      .(

        ID,
        Description,

        n_sig_diseases,

        min_FDR,

        max_abs_NES
      )
    ]
  )


  setorder(
    ranking,
    -n_sig_diseases,
    min_FDR,
    -max_abs_NES
  )


  ranking <- head(
    ranking,
    top_n
  )


  ranking[
    ,
    plot_label := Description
  ]


  duplicated_labels <- (

    duplicated(
      ranking$plot_label
    ) |

      duplicated(
        ranking$plot_label,
        fromLast = TRUE
      )
  )


  ranking[
    duplicated_labels,
    plot_label := paste0(
      Description,
      " [",
      ID,
      "]"
    )
  ]


  plot_dt <- merge(

    gsea_result$long,

    ranking[
      ,
      .(
        ID,
        plot_label
      )
    ],

    by = "ID",

    all = FALSE
  )


  plot_dt[
    ,
    plot_label := factor(
      plot_label,
      levels = rev(
        ranking$plot_label
      )
    )
  ]


  p <- ggplot(

    plot_dt,

    aes(
      x = dataset,
      y = plot_label,
      fill = NES
    )

  ) +

    geom_tile(
      color = "white",
      linewidth = 0.5
    ) +

    geom_text(
      aes(
        label = sprintf(
          "%.2f",
          NES
        )
      ),
      size = 3
    ) +

    scale_fill_gradient2(

      low = "#2166AC",

      mid = "white",

      high = "#B2182B",

      midpoint = 0,

      na.value = "grey95"
    ) +

    labs(

      title = paste0(
        "GSEA compartilhado - ",
        database_name
      ),

      subtitle = paste0(
        "Termos significativos (BH-FDR < ",
        alpha,
        ") em pelo menos duas doenças"
      ),

      x = NULL,

      y = NULL,

      fill = "NES"
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
      )
    )


  ggsave(

    file.path(
      plot_dir,
      paste0(
        "shared_GSEA_",
        database_name,
        ".png"
      )
    ),

    p,

    width = 9,

    height = max(
      6,
      0.32 *
        nrow(ranking) +
        2
    ),

    dpi = 300,

    bg = "white"
  )


  invisible(NULL)
}


plot_shared_gsea(
  gsea_go,
  "GO_BP"
)


plot_shared_gsea(
  gsea_kegg,
  "KEGG"
)


# ==============================================================================
# 29. FUNÇÕES DO RESUMO GSEA
# ==============================================================================

count_gsea_all3 <- function(
  gsea_result
) {

  if (nrow(gsea_result$wide) == 0) {

    return(0L)
  }


  sum(
    gsea_result$wide$n_sig_diseases == 3
  )
}


count_gsea_all3_same_direction <- function(
  gsea_result
) {

  if (nrow(gsea_result$wide) == 0) {

    return(0L)
  }


  sum(

    gsea_result$wide$n_sig_diseases == 3 &

      gsea_result$wide$same_direction == TRUE
  )
}


# ==============================================================================
# 30. RESUMO FINAL DA INTEGRAÇÃO
# ==============================================================================

integration_summary <- data.table(

  metric = c(

    "common_testable_universe",

    "primary_AD",
    "primary_HD",
    "primary_PD",

    "primary_AD_only",
    "primary_HD_only",
    "primary_PD_only",

    "primary_AD_HD",
    "primary_AD_PD",
    "primary_HD_PD",

    "primary_AD_HD_PD",

    "shared_HD_PD_same_direction",
    "shared_HD_PD_opposite_direction",

    "directionally_concordant_all_three",

    "sensitivity_AD",
    "sensitivity_HD",
    "sensitivity_PD",

    "sensitivity_AD_only",
    "sensitivity_HD_only",
    "sensitivity_PD_only",

    "sensitivity_AD_HD",
    "sensitivity_AD_PD",
    "sensitivity_HD_PD",

    "sensitivity_AD_HD_PD",

    "shared_GSEA_GO_BP_at_least_2",
    "shared_GSEA_GO_BP_all_3",
    "shared_GSEA_GO_BP_all_3_same_direction",

    "shared_GSEA_KEGG_at_least_2",
    "shared_GSEA_KEGG_all_3",
    "shared_GSEA_KEGG_all_3_same_direction",

    "shared_ORA_HD_PD_GO_BP",
    "shared_ORA_HD_PD_KEGG"
  ),


  value = c(

    length(
      common_universe
    ),


    sum(
      integrated$DEG_AD,
      na.rm = TRUE
    ),

    sum(
      integrated$DEG_HD,
      na.rm = TRUE
    ),

    sum(
      integrated$DEG_PD,
      na.rm = TRUE
    ),


    sum(
      integrated$primary_pattern ==
        "AD only"
    ),

    sum(
      integrated$primary_pattern ==
        "HD only"
    ),

    sum(
      integrated$primary_pattern ==
        "PD only"
    ),


    sum(
      integrated$primary_pattern ==
        "AD+HD"
    ),

    sum(
      integrated$primary_pattern ==
        "AD+PD"
    ),

    sum(
      integrated$primary_pattern ==
        "HD+PD"
    ),


    sum(
      integrated$primary_pattern ==
        "AD+HD+PD"
    ),


    nrow(
      shared_hd_pd_same_direction
    ),

    nrow(
      shared_hd_pd_opposite_direction
    ),


    nrow(
      concordant_three
    ),


    sum(
      integrated$DEG_sensitivity_AD,
      na.rm = TRUE
    ),

    sum(
      integrated$DEG_sensitivity_HD,
      na.rm = TRUE
    ),

    sum(
      integrated$DEG_sensitivity_PD,
      na.rm = TRUE
    ),


    sum(
      integrated$sensitivity_pattern ==
        "AD only"
    ),

    sum(
      integrated$sensitivity_pattern ==
        "HD only"
    ),

    sum(
      integrated$sensitivity_pattern ==
        "PD only"
    ),


    sum(
      integrated$sensitivity_pattern ==
        "AD+HD"
    ),

    sum(
      integrated$sensitivity_pattern ==
        "AD+PD"
    ),

    sum(
      integrated$sensitivity_pattern ==
        "HD+PD"
    ),


    sum(
      integrated$sensitivity_pattern ==
        "AD+HD+PD"
    ),


    nrow(
      gsea_go$wide
    ),

    count_gsea_all3(
      gsea_go
    ),

    count_gsea_all3_same_direction(
      gsea_go
    ),


    nrow(
      gsea_kegg$wide
    ),

    count_gsea_all3(
      gsea_kegg
    ),

    count_gsea_all3_same_direction(
      gsea_kegg
    ),


    nrow(
      ora_go_shared
    ),

    nrow(
      ora_kegg_shared
    )
  )
)


fwrite(

  integration_summary,

  file.path(
    output_dir,
    "integration_summary.tsv"
  ),

  sep = "\t",

  quote = FALSE,

  na = "NA"
)


# ==============================================================================
# 31. RESUMO NO TERMINAL
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("RESUMO FINAL - INTEGRAÇÃO TRANSDIAGNÓSTICA\n")
cat("============================================================\n")


print(
  integration_summary
)


cat("\n")
cat("Correlações de Spearman dos log2FC:\n")


print(
  correlation_summary
)


cat("\n")
cat("Figuras geradas:\n")

cat(
  " - primary_deg_patterns.png\n"
)

cat(
  " - venn_primary_AD_HD_PD.png\n"
)

cat(
  " - venn_sensitivity_AD_HD_PD.png\n"
)

cat(
  " - log2FC_HD_vs_PD.png\n"
)

cat(
  " - shared_GSEA_GO_BP.png\n"
)

cat(
  " - shared_GSEA_KEGG.png\n"
)


cat("\n")
cat("============================================================\n")
cat("07_cross_disease_integration.R CONCLUÍDO\n")
cat("============================================================\n")