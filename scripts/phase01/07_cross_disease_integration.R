#!/usr/bin/env Rscript

# ==============================================================================
# 07_cross_disease_integration.R
#
# INTEGRAÇÃO TRANSDIAGNÓSTICA — AD, HD e PD
#
# OBJETIVOS
# ------------------------------------------------------------------------------
# 1. Manter uma análise transdiagnóstica harmonizada baseada no universo comum
#    testável entre AD, HD e PD.
#
# 2. Manter os diagramas de Venn do universo comum testável:
#       - critério principal:
#           FDR < 0.05 e |log2FC| > 0.58
#       - análise de sensibilidade:
#           FDR < 0.05 e |log2FC| > 1
#
# 3. Acrescentar um Venn "article-like", comparável conceitualmente ao artigo:
#       p-value < 0.05 e |log2FC| > 1
#
#    IMPORTANTE:
#    esse Venn usa as LISTAS COMPLETAS de genes que satisfazem esse critério
#    em cada dataset, sem restringir ao universo comum AD∩HD∩PD.
#
# 4. Corrigir a análise de concordância HD vs PD:
#       - usar o universo pareado HD∩PD;
#       - não deixar AD determinar quais genes entram no scatter;
#       - calcular Spearman nesse universo pareado;
#       - mostrar rho e n no gráfico;
#       - adicionar linha y = x para referência de magnitude;
#       - distinguir DEGs de HD, PD, compartilhados e não-DEGs.
#
# 5. Calcular correlações de Spearman em universos PAREADOS:
#       AD∩HD, AD∩PD e HD∩PD.
#
# 6. Manter integração funcional por GSEA e ORA.
# 7. Auditar por GeneID os dez candidatos de Alharbi et al. (2025)
#    e a identidade exata da interseção nominal AD-HD-PD.
# 8. Detalhar os termos de ORA compartilhados em HD e PD, por All/Up/Down.
# 9. Realizar ORA exploratória com os genes da interseção nominal tripla,
#    usando o universo testável comum como fundo; não é análise confirmatória.
# 10. Garantir títulos legíveis nos três diagramas de Venn.
#
# CHAVE GÊNICA
# ------------------------------------------------------------------------------
# NCBI GeneID é a chave primária.
# Símbolos gênicos são utilizados somente como anotação.
#
# ENTRADAS
# ------------------------------------------------------------------------------
# results/deseq2/AD_deseq2_results.tsv
# results/deseq2/HD_deseq2_results.tsv
# results/deseq2/PD_deseq2_results.tsv
# results/deseq2/deseq2_summary.tsv
#
# results/enrichment/GSEA_GO_BP.tsv
# results/enrichment/GSEA_KEGG.tsv
# results/enrichment/ORA_GO_BP.tsv
# results/enrichment/ORA_KEGG.tsv
#
# SAÍDAS
# ------------------------------------------------------------------------------
# results/integration/
#   gene_integration.tsv
#   hd_pd_pairwise_integration.tsv
#   primary_pattern_summary.tsv
#   sensitivity_pattern_summary.tsv
#   venn_regions.tsv
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
#   venn_common_primary_AD_HD_PD.png
#   venn_common_sensitivity_AD_HD_PD.png
#   venn_article_like_nominal_AD_HD_PD.png
#   log2FC_HD_vs_PD_pairwise.png
#   shared_GSEA_GO_BP.png
#   shared_GSEA_KEGG.png
#
# OBSERVAÇÃO SOBRE CACHE
# ------------------------------------------------------------------------------
# Este script remove explicitamente os PNGs que ele próprio vai recriar antes
# de chamar ggsave(). O relatório HTML deve ser renderizado com cache = FALSE.
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

primary_lfc <- 0.58
sensitivity_lfc <- 1

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

scatter_colors <- c(
  "Not_DEG" = "#CFCFCF",
  "HD_only" = "#277DA1",
  "PD_only" = "#E75480",
  "HD+PD" = "#DF9203"
)

scatter_alpha <- c(
  "Not_DEG" = 0.28,
  "HD_only" = 0.62,
  "PD_only" = 0.62,
  "HD+PD" = 0.85
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

deseq_summary_file <- file.path(
  "results/deseq2",
  "deseq2_summary.tsv"
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
  deseq_summary_file,
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

  # Critério nominal comparável ao empregado no artigo:
  # p < 0.05 e |log2FC| > 1
  dt[
    ,
    DEG_article_like := (
      !is.na(pvalue) &
        pvalue < alpha &
        !is.na(log2FoldChange) &
        abs(log2FoldChange) > sensitivity_lfc
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
# 6.5. Resumo completo dos padrões
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
# 6.6. Calcular as sete regiões exclusivas de um Venn de três conjuntos
# ------------------------------------------------------------------------------

calculate_venn_regions <- function(sets) {

  if (!all(
    c(
      "AD",
      "HD",
      "PD"
    ) %in% names(sets)
  )) {
    stop(
      "calculate_venn_regions(): os conjuntos devem se chamar AD, HD e PD."
    )
  }

  AD <- unique(
    as.character(
      sets[["AD"]]
    )
  )

  HD <- unique(
    as.character(
      sets[["HD"]]
    )
  )

  PD <- unique(
    as.character(
      sets[["PD"]]
    )
  )

  triple <- Reduce(
    intersect,
    list(
      AD,
      HD,
      PD
    )
  )

  AD_HD <- setdiff(
    intersect(
      AD,
      HD
    ),
    triple
  )

  AD_PD <- setdiff(
    intersect(
      AD,
      PD
    ),
    triple
  )

  HD_PD <- setdiff(
    intersect(
      HD,
      PD
    ),
    triple
  )

  AD_only <- setdiff(
    AD,
    union(
      HD,
      PD
    )
  )

  HD_only <- setdiff(
    HD,
    union(
      AD,
      PD
    )
  )

  PD_only <- setdiff(
    PD,
    union(
      AD,
      HD
    )
  )

  regions <- data.table(
    pattern = c(
      "AD only",
      "HD only",
      "PD only",
      "AD+HD",
      "AD+PD",
      "HD+PD",
      "AD+HD+PD"
    ),
    n_genes = c(
      length(AD_only),
      length(HD_only),
      length(PD_only),
      length(AD_HD),
      length(AD_PD),
      length(HD_PD),
      length(triple)
    )
  )

  totals <- data.table(
    dataset = c(
      "AD",
      "HD",
      "PD"
    ),
    n_genes = c(
      length(AD),
      length(HD),
      length(PD)
    )
  )

  reconstructed_AD <- (
    length(AD_only) +
      length(AD_HD) +
      length(AD_PD) +
      length(triple)
  )

  reconstructed_HD <- (
    length(HD_only) +
      length(AD_HD) +
      length(HD_PD) +
      length(triple)
  )

  reconstructed_PD <- (
    length(PD_only) +
      length(AD_PD) +
      length(HD_PD) +
      length(triple)
  )

  if (reconstructed_AD != length(AD)) {
    stop(
      "Erro interno no Venn: regiões de AD não reproduzem o total de AD."
    )
  }

  if (reconstructed_HD != length(HD)) {
    stop(
      "Erro interno no Venn: regiões de HD não reproduzem o total de HD."
    )
  }

  if (reconstructed_PD != length(PD)) {
    stop(
      "Erro interno no Venn: regiões de PD não reproduzem o total de PD."
    )
  }

  return(
    list(
      regions = regions,
      totals = totals
    )
  )
}


# ------------------------------------------------------------------------------
# 6.7. Criar diagrama de Venn customizado
#
# IMPORTANTE:
# os círculos têm o mesmo tamanho.
# A área NÃO é proporcional ao número de genes.
# ------------------------------------------------------------------------------

create_venn_plot <- function(
  venn_result,
  title,
  subtitle
) {

  dt <- copy(
    venn_result$regions
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

  # Evitar corte dos titulos longos no PNG e no relatorio HTML.
  titulo_curto <- paste(strwrap(title, width = 38), collapse = "\n")
  subtitulo_curto <- paste(strwrap(subtitle, width = 75), collapse = "\n")

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
      title = titulo_curto,
      subtitle = subtitulo_curto
    ) +

    theme_void() +

    theme(
      legend.position = "none",

      plot.title = element_text(
        size = 17,
        face = "bold",
        hjust = 0.5,
        lineheight = 1.12,
        margin = margin(b = 7)
      ),

      plot.subtitle = element_text(
        size = 11,
        hjust = 0.5,
        lineheight = 1.15,
        margin = margin(b = 10)
      ),

      plot.title.position = "plot",
      plot.margin = margin(35, 65, 30, 65)
    )

  return(p)
}


# ------------------------------------------------------------------------------
# 6.8. Preparar resultado de um dataset para merge
# ------------------------------------------------------------------------------

prepare_for_merge <- function(
  dt,
  prefix,
  universe
) {

  x <- dt[
    gene_id %in% universe,
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
      direction_sensitivity,
      DEG_article_like
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


# ------------------------------------------------------------------------------
# 6.9. Correlação pareada
# ------------------------------------------------------------------------------

compute_pairwise_correlation <- function(
  dataset_a,
  dataset_b,
  results_list,
  testable_ids
) {

  pair_universe <- intersect(
    testable_ids[[dataset_a]],
    testable_ids[[dataset_b]]
  )

  a <- results_list[[dataset_a]][
    gene_id %in% pair_universe,
    .(
      gene_id,
      log2FoldChange
    )
  ]

  b <- results_list[[dataset_b]][
    gene_id %in% pair_universe,
    .(
      gene_id,
      log2FoldChange
    )
  ]

  setnames(
    a,
    "log2FoldChange",
    "lfc_a"
  )

  setnames(
    b,
    "log2FoldChange",
    "lfc_b"
  )

  pair_dt <- merge(
    a,
    b,
    by = "gene_id",
    all = FALSE,
    sort = FALSE
  )

  rho <- cor(
    pair_dt$lfc_a,
    pair_dt$lfc_b,
    method = "spearman",
    use = "complete.obs"
  )

  return(
    data.table(
      pair = paste0(
        dataset_a,
        "-",
        dataset_b
      ),
      universe = paste0(
        dataset_a,
        "∩",
        dataset_b
      ),
      n_genes = nrow(
        pair_dt
      ),
      spearman_rho = rho
    )
  )
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
# 8. AUDITAR TOTAIS DO 05_DESEQ2.R
# ==============================================================================

deseq_summary <- fread(
  deseq_summary_file
)

required_summary_columns <- c(
  "dataset",
  "n_DEG_primary",
  "n_DEG_sensitivity"
)

missing_summary_columns <- setdiff(
  required_summary_columns,
  names(
    deseq_summary
  )
)

if (length(missing_summary_columns) > 0) {
  stop(
    "deseq2_summary.tsv não contém: ",
    paste(
      missing_summary_columns,
      collapse = ", "
    )
  )
}

for (dataset_name in datasets) {

  dt <- results_list[[dataset_name]]

  observed_primary <- sum(
    dt$DEG,
    na.rm = TRUE
  )

  observed_sensitivity <- sum(
    dt$DEG_sensitivity,
    na.rm = TRUE
  )

  expected_primary <- deseq_summary[
    dataset == dataset_name,
    n_DEG_primary
  ]

  expected_sensitivity <- deseq_summary[
    dataset == dataset_name,
    n_DEG_sensitivity
  ]

  if (length(expected_primary) != 1) {
    stop(
      dataset_name,
      ": número inesperado de linhas em deseq2_summary.tsv."
    )
  }

  if (observed_primary != expected_primary) {
    stop(
      dataset_name,
      ": DEG no arquivo completo = ",
      observed_primary,
      ", mas deseq2_summary.tsv informa ",
      expected_primary,
      "."
    )
  }

  if (observed_sensitivity != expected_sensitivity) {
    stop(
      dataset_name,
      ": DEG_sensitivity no arquivo completo = ",
      observed_sensitivity,
      ", mas deseq2_summary.tsv informa ",
      expected_sensitivity,
      "."
    )
  }
}


# ==============================================================================
# 9. UNIVERSOS TESTÁVEIS
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

common_universe <- Reduce(
  intersect,
  testable_ids
)

hd_pd_universe <- intersect(
  testable_ids[["HD"]],
  testable_ids[["PD"]]
)

cat("\n")
cat("============================================================\n")
cat("UNIVERSOS TESTÁVEIS\n")
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
  "AD ∩ HD ∩ PD: ",
  length(common_universe),
  " genes\n",
  sep = ""
)

cat(
  "HD ∩ PD: ",
  length(hd_pd_universe),
  " genes\n",
  sep = ""
)


# ==============================================================================
# 10. MATRIZ INTEGRADA — UNIVERSO COMUM AD ∩ HD ∩ PD
# ==============================================================================

ad_common <- prepare_for_merge(
  results_list[["AD"]],
  "AD",
  common_universe
)

hd_common <- prepare_for_merge(
  results_list[["HD"]],
  "HD",
  common_universe
)

pd_common <- prepare_for_merge(
  results_list[["PD"]],
  "PD",
  common_universe
)

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
    ad_common,
    hd_common,
    pd_common
  )
)

if (
  nrow(integrated) !=
    length(common_universe)
) {
  stop(
    "Matriz comum possui ",
    nrow(integrated),
    " genes, mas common_universe possui ",
    length(common_universe),
    "."
  )
}

integrated[
  ,
  symbol := first_nonmissing(
    symbol_AD,
    symbol_HD,
    symbol_PD
  )
]


# ==============================================================================
# 11. PADRÕES DE DEGs — UNIVERSO COMUM
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
# 12. DIREÇÃO DOS EFEITOS — UNIVERSO COMUM
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
# 13. RESUMOS DOS PADRÕES — UNIVERSO COMUM
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
# 14. VENNS — UNIVERSO COMUM TESTÁVEL
# ==============================================================================

common_primary_sets <- list(
  AD = integrated[
    DEG_AD == TRUE,
    gene_id
  ],
  HD = integrated[
    DEG_HD == TRUE,
    gene_id
  ],
  PD = integrated[
    DEG_PD == TRUE,
    gene_id
  ]
)

common_sensitivity_sets <- list(
  AD = integrated[
    DEG_sensitivity_AD == TRUE,
    gene_id
  ],
  HD = integrated[
    DEG_sensitivity_HD == TRUE,
    gene_id
  ],
  PD = integrated[
    DEG_sensitivity_PD == TRUE,
    gene_id
  ]
)

venn_common_primary <- calculate_venn_regions(
  common_primary_sets
)

venn_common_sensitivity <- calculate_venn_regions(
  common_sensitivity_sets
)


# ==============================================================================
# 15. VENN — COMPARAÇÃO NOMINAL "ARTICLE-LIKE"
#
# p < 0.05 e |log2FC| > 1
#
# Usa as listas COMPLETAS de cada dataset.
# ==============================================================================

article_like_sets <- lapply(
  results_list,
  function(dt) {

    unique(
      dt[
        DEG_article_like == TRUE,
        gene_id
      ]
    )
  }
)

venn_article_like <- calculate_venn_regions(
  article_like_sets
)



# ==============================================================================
# 15b. GENES DA INTERSECAO NOMINAL E COMPARACAO COM ALHARBI ET AL. (2025)
#
# Tabela 3 do artigo: DOI 10.3390/cimb47120976.
# Os efeitos do artigo sao apresentados como "logFC"; nao se presume que a
# escala seja diretamente comparavel ao log2FoldChange do nosso DESeq2.
# A correspondencia gene a gene e feita EXCLUSIVAMENTE por NCBI GeneID.
# O simbolo serve para localizar os IDs candidatos na anotacao e indicar
# possiveis ambiguidades / discrepancias entre as anotacoes.
# ==============================================================================

artigo_2025 <- data.table(
  gene_artigo = c("H19", "CCL2", "CSF3", "IL17REL", "MMP9",
                  "PDLIM1", "MMRN1", "SLPI", "S100A8", "S100A9"),
  logFC_artigo_AD = c(-2.335, -2.810, -1.013, -1.230, -1.508,
                       1.135, -1.262, -1.159, -1.674, -1.351),
  logFC_artigo_HD = c( 1.887, -1.743,  2.319, -1.612,  3.235,
                       1.617, -2.954, -2.718,  1.430,  1.032),
  logFC_artigo_PD = c(-1.2741803, -1.4988129, -1.28162, 1.5461572,
                      -1.1023493, -1.2294081, -1.2365235, -1.9889283,
                      -1.49514, -1.0651068),
  p_artigo_AD = c(1.20e-4,8.75e-4,4.87e-2,3.36e-2,2.53e-2,
                  2.43e-2,3.16e-2,3.47e-2,4.20e-3,2.19e-2),
  p_artigo_HD = c(1.93e-7,9.45e-5,3.79e-3,1.30e-6,3.01e-11,
                  5.53e-3,4.68e-4,1.00e-12,2.07e-4,1.79e-3),
  p_artigo_PD = c(1.08e-6,8.06e-5,4.45e-2,1.49e-6,2.68e-6,
                  3.55e-4,5.85e-4,4.52e-6,1.02e-3,2.00e-2)
)

# Em caso de divergencia entre anotacoes, nao selecionar arbitrariamente IDs.
# O OrgDb usado no enriquecimento da Fase 01 e a fonte preferencial.
resolver_geneid <- function(simbolo) {
  ids_db <- character(0)
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    ids_db <- tryCatch({
      z <- suppressMessages(AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = simbolo, column = "ENTREZID", keytype = "SYMBOL",
        multiVals = "list"
      ))
      unique(as.character(unlist(z, use.names = FALSE)))
    }, error = function(e) character(0))
    ids_db <- ids_db[!is.na(ids_db) & nzchar(ids_db)]
  }
  ids_locais <- unique(unlist(lapply(results_list, function(dt) {
    dt[!is.na(symbol) & toupper(trimws(symbol)) == toupper(simbolo), gene_id]
  }), use.names = FALSE))
  ids_locais <- as.character(ids_locais[!is.na(ids_locais)])

  if (length(ids_db) == 1L) {
    status <- if (length(ids_locais) > 0L &&
                  !ids_db[[1]] %in% ids_locais) {
      "divergencia_orgdb_vs_simbolos_locais"
    } else "mapeamento_orgdb_univoco"
    return(list(gene_id = ids_db[[1]], status = status,
                ids_locais = paste(sort(ids_locais), collapse = ";")))
  }
  if (length(ids_db) > 1L) {
    return(list(gene_id = NA_character_, status = "orgdb_ambiguo_revisar",
                ids_locais = paste(sort(ids_locais), collapse = ";")))
  }
  if (length(ids_locais) == 1L) {
    return(list(gene_id = ids_locais[[1]], status = "simbolo_local_univoco_sem_orgdb",
                ids_locais = ids_locais[[1]]))
  }
  status <- if (length(ids_locais) > 1L) "simbolo_local_ambiguo_revisar" else
    "simbolo_nao_localizado"
  list(gene_id = NA_character_, status = status,
       ids_locais = paste(sort(ids_locais), collapse = ";"))
}

resolvidos <- lapply(artigo_2025$gene_artigo, resolver_geneid)
artigo_2025[, gene_id := vapply(resolvidos, `[[`, character(1), "gene_id")]
artigo_2025[, status_mapeamento := vapply(resolvidos, `[[`, character(1), "status")]
artigo_2025[, gene_ids_achados_nos_arquivos :=
              vapply(resolvidos, `[[`, character(1), "ids_locais")]

# Inclui os dez genes do artigo MESMO quando ausentes/não testaveis nas nossas
# matrizes; NA significa indisponibilidade, nunca ausencia de efeito.
audit_geneids <- function(ids) {
  d <- data.table(gene_id = as.character(ids))
  for (doenca in datasets) {
    x <- results_list[[doenca]][, .(
      gene_id, symbol, log2FoldChange, pvalue, padj,
      DEG_article_like, DEG, DEG_sensitivity
    )]
    setnames(x, setdiff(names(x), "gene_id"),
             paste0(setdiff(names(x), "gene_id"), "_", doenca))
    d <- merge(d, x, by = "gene_id", all.x = TRUE, sort = FALSE)
    effect_col <- paste0("log2FoldChange_", doenca)
    p_col <- paste0("pvalue_", doenca)
    d[, (paste0("testavel_", doenca)) :=
        is.finite(get(effect_col)) & is.finite(get(p_col))]
    d[, (paste0("direcao_", doenca)) := effect_direction(get(effect_col))]
    d[, (paste0("status_nominal_", doenca)) := fcase(
      is.na(get(paste0("symbol_", doenca))) & is.na(get(effect_col)),
      "ausente_ou_sem_estatistica",
      !get(paste0("testavel_", doenca)), "nao_testavel",
      get(paste0("DEG_article_like_", doenca)) == TRUE, "DEG_nominal",
      default = "nao_DEG_nominal"
    )]
  }
  d[, simbolo_nosso := first_nonmissing(symbol_AD, symbol_HD, symbol_PD)]
  d
}

ids_triplo_nominal <- sort(Reduce(intersect, article_like_sets))
if (length(ids_triplo_nominal) != 8L) {
  warning("Interseção nominal da Fase 01: ", length(ids_triplo_nominal),
          " GeneIDs (o relatório anterior tinha 8); verificar versões dos dados.")
}

genes_triplo_nominal <- audit_geneids(ids_triplo_nominal)
genes_triplo_nominal[, no_artigo_2025 :=
                      gene_id %in% artigo_2025$gene_id[!is.na(artigo_2025$gene_id)]]
genes_triplo_nominal[, mesma_direcao_nas_tres :=
                      direcao_AD == direcao_HD & direcao_HD == direcao_PD &
                      direcao_AD %in% c("Up", "Down")]
setorder(genes_triplo_nominal, gene_id)

artigo_2025[, ordem_original := seq_len(.N)]
artigo_auditado <- merge(
  artigo_2025,
  audit_geneids(unique(artigo_2025$gene_id[!is.na(artigo_2025$gene_id)])),
  by = "gene_id", all.x = TRUE, sort = FALSE
)
setorder(artigo_auditado, ordem_original)
artigo_auditado[, ordem_original := NULL]

for (doenca in datasets) {
  observado <- artigo_auditado[[paste0("log2FoldChange_", doenca)]]
  artigo <- artigo_auditado[[paste0("logFC_artigo_", doenca)]]
  artigo_auditado[, (paste0("mesma_direcao_artigo_", doenca)) :=
                   ifelse(is.finite(observado) & is.finite(artigo),
                          sign(observado) == sign(artigo), NA)]
}
artigo_auditado[, reproducao_nominal_tripla :=
                   !is.na(gene_id) & gene_id %in% ids_triplo_nominal]
artigo_auditado[, concordancia_direcional_tripla :=
   mesma_direcao_artigo_AD == TRUE & mesma_direcao_artigo_HD == TRUE &
   mesma_direcao_artigo_PD == TRUE]

comparacao_artigo_resumo <- data.table(
  indicador = c("artigo_nominais_AD", "artigo_nominais_HD", "artigo_nominais_PD",
     "nossos_nominais_AD", "nossos_nominais_HD", "nossos_nominais_PD",
     "artigo_intersecao_tripla", "nossa_intersecao_tripla",
     "genes_do_artigo_reproduzidos_nominal", "genes_nossos_ausentes_artigo",
     "genes_artigo_sem_mapeamento_univoco",
     "genes_artigo_reproduzidos_mesma_direcao_tripla"),
  valor = c(262L, 1581L, 537L,
      length(article_like_sets$AD), length(article_like_sets$HD),
      length(article_like_sets$PD), 10L, length(ids_triplo_nominal),
      sum(artigo_auditado$reproducao_nominal_tripla, na.rm = TRUE),
      sum(!genes_triplo_nominal$no_artigo_2025, na.rm = TRUE),
      sum(is.na(artigo_auditado$gene_id)),
      sum(artigo_auditado$reproducao_nominal_tripla &
          artigo_auditado$concordancia_direcional_tripla, na.rm = TRUE))
)

fwrite(genes_triplo_nominal,
       file.path(output_dir, "genes_intersecao_nominal_tripla.tsv"), sep="\t", na="NA")
fwrite(artigo_auditado,
       file.path(output_dir, "dez_genes_Alharbi_2025_comparacao.tsv"), sep="\t", na="NA")
fwrite(comparacao_artigo_resumo,
       file.path(output_dir, "comparacao_Alharbi_2025_resumo.tsv"), sep="\t", na="NA")

cat("\nInterseção nominal tripla — nomes e GeneIDs:\n")
print(genes_triplo_nominal[, .(gene_id, simbolo_nosso, no_artigo_2025,
          direcao_AD, direcao_HD, direcao_PD)])
cat("\nCandidatos de Alharbi 2025 — reprodução nominal:\n")
print(artigo_auditado[, .(gene_artigo, gene_id, status_mapeamento,
               reproducao_nominal_tripla, concordancia_direcional_tripla)])

# ==============================================================================
# 16. AUDITORIA DOS VENNS
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("AUDITORIA DOS DIAGRAMAS DE VENN\n")
cat("============================================================\n")

cat("\nVenn — universo comum / critério principal:\n")
print(
  venn_common_primary$totals
)
print(
  venn_common_primary$regions
)

cat("\nVenn — universo comum / sensibilidade:\n")
print(
  venn_common_sensitivity$totals
)
print(
  venn_common_sensitivity$regions
)

cat("\nVenn — listas completas / critério nominal article-like:\n")
print(
  venn_article_like$totals
)
print(
  venn_article_like$regions
)

cat("\n")


# ==============================================================================
# 17. MATRIZ PAREADA HD ∩ PD
# ==============================================================================

hd_pair <- prepare_for_merge(
  results_list[["HD"]],
  "HD",
  hd_pd_universe
)

pd_pair <- prepare_for_merge(
  results_list[["PD"]],
  "PD",
  hd_pd_universe
)

hd_pd_pairwise <- merge(
  hd_pair,
  pd_pair,
  by = "gene_id",
  all = FALSE,
  sort = FALSE
)

if (
  nrow(hd_pd_pairwise) !=
    length(hd_pd_universe)
) {
  stop(
    "Matriz HD∩PD possui ",
    nrow(hd_pd_pairwise),
    " genes, mas hd_pd_universe possui ",
    length(hd_pd_universe),
    "."
  )
}

hd_pd_pairwise[
  ,
  symbol := first_nonmissing(
    symbol_HD,
    symbol_PD
  )
]

hd_pd_pairwise[
  ,
  effect_HD := effect_direction(
    log2FoldChange_HD
  )
]

hd_pd_pairwise[
  ,
  effect_PD := effect_direction(
    log2FoldChange_PD
  )
]

hd_pd_pairwise[
  ,
  primary_class := fcase(

    DEG_HD &
      DEG_PD,
    "HD+PD",

    DEG_HD,
    "HD_only",

    DEG_PD,
    "PD_only",

    default = "Not_DEG"
  )
]

hd_pd_pairwise[
  ,
  same_direction := (
    effect_HD == effect_PD &
      effect_HD != "Zero" &
      !is.na(effect_HD)
  )
]


# ==============================================================================
# 18. DEGs COMPARTILHADOS HD + PD — UNIVERSO PAREADO
# ==============================================================================

shared_hd_pd <- hd_pd_pairwise[
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
# 20. CORRELAÇÕES DE SPEARMAN — UNIVERSOS PAREADOS
# ==============================================================================

correlation_summary <- rbindlist(
  list(
    compute_pairwise_correlation(
      "AD",
      "HD",
      results_list,
      testable_ids
    ),
    compute_pairwise_correlation(
      "AD",
      "PD",
      results_list,
      testable_ids
    ),
    compute_pairwise_correlation(
      "HD",
      "PD",
      results_list,
      testable_ids
    )
  )
)

rho_hd_pd <- correlation_summary[
  pair == "HD-PD",
  spearman_rho
]

n_hd_pd_pair <- correlation_summary[
  pair == "HD-PD",
  n_genes
]


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

integrate_ora_hd_pd <- function(file_path, database_name) {
  dt <- fread(file_path)
  required <- c("dataset", "gene_set", "ID", "Description", "p.adjust", "significant")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0L) {
    stop(database_name, " ORA: colunas ausentes: ", paste(missing, collapse = ", "))
  }
  if (!"geneID" %in% names(dt)) dt[, geneID := NA_character_]
  if (!"Count" %in% names(dt)) dt[, Count := NA_integer_]
  dt[, significant := as_logical_safe(significant)]
  sig <- dt[significant == TRUE & dataset %in% c("HD", "PD") &
            !is.na(ID) & is.finite(p.adjust)]
  # Se um mesmo ID constar mais de uma vez, preservar o menor FDR dentro
  # da mesma categoria e do mesmo dataset, sem contar duplicatas.
  setorder(sig, p.adjust)
  sig <- unique(sig, by = c("dataset", "gene_set", "ID"))

  empty <- data.table(gene_set = character(), ID = character(),
     Description = character(), n_datasets = integer(),
     datasets_shared = character(), best_FDR = numeric(),
     FDR_HD = numeric(), FDR_PD = numeric(),
     geneID_HD = character(), geneID_PD = character())
  if (nrow(sig) == 0L) return(empty)

  h <- sig[dataset == "HD", .(gene_set, ID, Description_HD=Description,
                                FDR_HD=p.adjust, geneID_HD=geneID, Count_HD=Count)]
  p <- sig[dataset == "PD", .(gene_set, ID, Description_PD=Description,
                                FDR_PD=p.adjust, geneID_PD=geneID, Count_PD=Count)]
  if (nrow(h) == 0L || nrow(p) == 0L) return(empty)
  shared <- merge(h, p, by = c("gene_set", "ID"), all = FALSE)
  if (nrow(shared) == 0L) return(empty)
  shared[, Description := fifelse(!is.na(Description_HD) & nzchar(Description_HD),
                                 Description_HD, Description_PD)]
  shared[, `:=`(n_datasets=2L, datasets_shared="HD+PD",
                best_FDR=pmin(FDR_HD, FDR_PD))]
  setcolorder(shared, c("gene_set", "ID", "Description", "n_datasets",
           "datasets_shared", "best_FDR", "FDR_HD", "FDR_PD",
           "Count_HD", "Count_PD", "geneID_HD", "geneID_PD",
           "Description_HD", "Description_PD"))
  setorder(shared, gene_set, best_FDR)
  shared
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
# 22b. ORA COMPLEMENTAR DOS GENES DA INTERSECAO NOMINAL TRIPLA
#
# Comparacao EXPLORATORIA: os genes nao foram selecionados por FDR.
# Universo de fundo: genes testaveis nos tres datasets (common_universe).
# Evitar interpretar eventual ausencia de termos como ausencia de biologia.
# KEGG pode depender de acesso remoto e falhar; a tabela de status distingue
# falha de execucao de ausencia de termos significativos.
# ==============================================================================

bases_ora <- c("GO_BP", "KEGG")
nominal_ora_status <- data.table(
  base = bases_ora,
  situacao = "nao_executada",
  mensagem = "",
  n_genes_entrada = length(ids_triplo_nominal),
  n_universo = length(common_universe),
  n_termos_testados = NA_integer_,
  n_termos_FDR005 = NA_integer_
)

# A tabela sempre tem cabecalho, mesmo quando nao ha resultados.
modelo_ora_nominal <- data.table(
  database=character(), ID=character(), Description=character(),
  GeneRatio=character(), BgRatio=character(), pvalue=numeric(),
  p.adjust=numeric(), qvalue=numeric(), Count=integer(),
  geneID=character(), significant=logical()
)

executar_ora_nominal <- function(base) {
  if (!requireNamespace("clusterProfiler", quietly=TRUE) ||
      (base == "GO_BP" && !requireNamespace("org.Hs.eg.db", quietly=TRUE))) {
    return(list(situacao="pacote_indisponivel",
                mensagem="Requer clusterProfiler; GO requer org.Hs.eg.db",
                tabela=copy(modelo_ora_nominal)))
  }
  if (length(ids_triplo_nominal) == 0L) {
    return(list(situacao="sem_genes_entrada", mensagem="Interseção vazia",
                tabela=copy(modelo_ora_nominal)))
  }
  z <- tryCatch({
    if (base == "GO_BP") {
      clusterProfiler::enrichGO(
        gene = ids_triplo_nominal, universe = common_universe,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
        ont = "BP", pAdjustMethod = "BH", pvalueCutoff=1,
        qvalueCutoff=1, minGSSize=10, maxGSSize=500,
        readable=FALSE)
    } else {
      clusterProfiler::enrichKEGG(
        gene = ids_triplo_nominal, universe = common_universe,
        organism = "hsa", keyType="ncbi-geneid", pAdjustMethod="BH",
        pvalueCutoff=1, qvalueCutoff=1, minGSSize=10,
        maxGSSize=500, use_internal_data=FALSE)
    }
  }, error=function(e) e)
  if (inherits(z, "error")) {
    return(list(situacao="erro_execucao", mensagem=conditionMessage(z),
                tabela=copy(modelo_ora_nominal)))
  }
  if (is.null(z)) {
    return(list(situacao="executada_sem_termos_testaveis",
                mensagem="Nenhum termo retornado pelos filtros da base",
                tabela=copy(modelo_ora_nominal)))
  }
  tab <- as.data.table(as.data.frame(z))
  if (nrow(tab) == 0L) {
    return(list(situacao="executada_sem_termos_testaveis",
                mensagem="Nenhum termo retornado pelos filtros da base",
                tabela=copy(modelo_ora_nominal)))
  }
  tab[, database := base]
  tab[, significant := is.finite(p.adjust) & p.adjust < alpha]
  # Conservar todos os termos testados, nao apenas FDR-significativos.
  for (col in setdiff(names(modelo_ora_nominal), names(tab))) {
    tab[, (col) := NA]
  }
  setcolorder(tab, c(names(modelo_ora_nominal),
                      setdiff(names(tab), names(modelo_ora_nominal))))
  setorder(tab, p.adjust, pvalue)
  list(situacao="executada", mensagem="ORA nominal comparativa",
       tabela=tab)
}

ora_nominal_results <- lapply(bases_ora, executar_ora_nominal)
names(ora_nominal_results) <- bases_ora
for (base in bases_ora) {
  z <- ora_nominal_results[[base]]
  i <- match(base, nominal_ora_status$base)
  nominal_ora_status[i, `:=`(
    situacao=z$situacao, mensagem=z$mensagem,
    n_termos_testados=if (startsWith(z$situacao,"executada")) nrow(z$tabela) else NA_integer_,
    n_termos_FDR005=if (startsWith(z$situacao,"executada"))
       sum(z$tabela$significant, na.rm=TRUE) else NA_integer_
  )]
  fwrite(z$tabela,
         file.path(output_dir, paste0("ORA_exploratoria_8_genes_nominais_",base,".tsv")),
         sep="\t", na="NA")
}
fwrite(nominal_ora_status,
       file.path(output_dir, "ORA_exploratoria_nominal_status.tsv"),
       sep="\t", na="NA")

# Quantificar o que foi efetivamente compartilhado na ORA principal HD-PD,
# distinguindo termos All / Up / Down para nao somar a mesma via varias vezes.
ora_hd_pd_resumo <- rbindlist(lapply(bases_ora, function(base) {
  x <- if (base == "GO_BP") ora_go_shared else ora_kegg_shared
  data.table(database=base, gene_set=c("All","Up","Down"),
    n_termos_compartilhados=vapply(c("All","Up","Down"),
      function(g) sum(x$gene_set==g,na.rm=TRUE),integer(1)))
}))
fwrite(ora_hd_pd_resumo,
       file.path(output_dir, "ORA_compartilhada_HD_PD_resumo_por_conjunto.tsv"),
       sep="\t", na="NA")

# ==============================================================================
# 23. TABELA ÚNICA DAS REGIÕES DOS VENNS
# ==============================================================================

venn_regions <- rbindlist(
  list(
    copy(
      venn_common_primary$regions
    )[
      ,
      `:=`(
        venn = "common_primary",
        definition = "Common testable universe; FDR<0.05; |log2FC|>0.58"
      )
    ],

    copy(
      venn_common_sensitivity$regions
    )[
      ,
      `:=`(
        venn = "common_sensitivity",
        definition = "Common testable universe; FDR<0.05; |log2FC|>1"
      )
    ],

    copy(
      venn_article_like$regions
    )[
      ,
      `:=`(
        venn = "article_like_nominal",
        definition = "Full DEG lists; p<0.05; |log2FC|>1"
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

setcolorder(
  venn_regions,
  c(
    "venn",
    "definition",
    "pattern",
    "n_genes"
  )
)


# ==============================================================================
# 24. SALVAR TABELAS
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
  hd_pd_pairwise,
  file.path(
    output_dir,
    "hd_pd_pairwise_integration.tsv"
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
  venn_regions,
  file.path(
    output_dir,
    "venn_regions.tsv"
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
# 25. LIMPAR FIGURAS QUE SERÃO RECRIADAS
# ==============================================================================

plot_files <- c(
  "primary_deg_patterns.png",
  "venn_common_primary_AD_HD_PD.png",
  "venn_common_sensitivity_AD_HD_PD.png",
  "venn_article_like_nominal_AD_HD_PD.png",
  "log2FC_HD_vs_PD_pairwise.png",
  "shared_GSEA_GO_BP.png",
  "shared_GSEA_KEGG.png"
)

plot_paths <- file.path(
  plot_dir,
  plot_files
)

existing_plots <- plot_paths[
  file.exists(
    plot_paths
  )
]

if (length(existing_plots) > 0) {
  unlink(
    existing_plots,
    force = TRUE
  )
}


# ==============================================================================
# 26. GRÁFICO DOS PADRÕES DE DEGs — UNIVERSO COMUM
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
      "Critério principal; universo comum testável AD∩HD∩PD = ",
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
# 27. VENN — UNIVERSO COMUM / CRITÉRIO PRINCIPAL
# ==============================================================================

p_venn_common_primary <- create_venn_plot(
  venn_result = venn_common_primary,

  title = "DEGs no universo comum — critério principal",

  subtitle = paste0(
    "FDR < 0,05 e |log2FC| > 0,58; AD∩HD∩PD = ",
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
    "venn_common_primary_AD_HD_PD.png"
  ),
  p_venn_common_primary,
  width = 12.5,
  height = 9,
  dpi = 300,
  bg = "white"
)


# ==============================================================================
# 28. VENN — UNIVERSO COMUM / SENSIBILIDADE
# ==============================================================================

p_venn_common_sensitivity <- create_venn_plot(
  venn_result = venn_common_sensitivity,

  title = "DEGs no universo comum — critério estrito",

  subtitle = paste0(
    "FDR < 0,05 e |log2FC| > 1; AD∩HD∩PD = ",
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
    "venn_common_sensitivity_AD_HD_PD.png"
  ),
  p_venn_common_sensitivity,
  width = 12.5,
  height = 9,
  dpi = 300,
  bg = "white"
)


# ==============================================================================
# 29. VENN — COMPARAÇÃO NOMINAL "ARTICLE-LIKE"
# ==============================================================================

article_totals <- setNames(
  venn_article_like$totals$n_genes,
  venn_article_like$totals$dataset
)

p_venn_article_like <- create_venn_plot(
  venn_result = venn_article_like,

  title = "Interseção nominal — comparação com o artigo",

  subtitle = paste0(
    "p < 0,05 e |log2FC| > 1; listas completas: AD = ",
    article_totals[["AD"]],
    ", HD = ",
    format(
      article_totals[["HD"]],
      big.mark = "."
    ),
    ", PD = ",
    format(
      article_totals[["PD"]],
      big.mark = "."
    )
  )
)

ggsave(
  file.path(
    plot_dir,
    "venn_article_like_nominal_AD_HD_PD.png"
  ),
  p_venn_article_like,
  width = 12.5,
  height = 9,
  dpi = 300,
  bg = "white"
)


# ==============================================================================
# 30. SCATTER — CONCORDÂNCIA DOS EFEITOS HD vs PD
#
# Agora usa o universo pareado HD∩PD.
# ==============================================================================

scatter_dt <- copy(
  hd_pd_pairwise
)

scatter_dt[
  ,
  primary_class := factor(
    primary_class,
    levels = c(
      "Not_DEG",
      "HD_only",
      "PD_only",
      "HD+PD"
    )
  )
]

p_hd_pd <- ggplot(
  scatter_dt,
  aes(
    x = log2FoldChange_HD,
    y = log2FoldChange_PD
  )
) +

  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "solid",
    linewidth = 0.45,
    color = "#7A7A7A"
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
    aes(
      color = primary_class,
      alpha = primary_class
    ),
    size = 1.35
  ) +

  scale_color_manual(
    values = scatter_colors,
    labels = c(
      "Not_DEG" = "Not_DEG",
      "HD_only" = "HD",
      "PD_only" = "PD",
      "HD+PD" = "HD+PD"
    )
  ) +

  scale_alpha_manual(
    values = scatter_alpha,
    guide = "none"
  ) +

  labs(
    title = "Concordância de efeito — HD vs PD",

    subtitle = paste0(
      "Universo pareado HD∩PD = ",
      format(
        n_hd_pd_pair,
        big.mark = "."
      ),
      " genes; Spearman ρ = ",
      sprintf(
        "%.3f",
        rho_hd_pd
      ),
      "; DEGs compartilhados = ",
      nrow(
        shared_hd_pd
      )
    ),

    x = "log2FoldChange — HD",
    y = "log2FoldChange — PD",
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
    "log2FC_HD_vs_PD_pairwise.png"
  ),
  p_hd_pd,
  width = 8.5,
  height = 6.8,
  dpi = 300,
  bg = "white"
)


# ==============================================================================
# 31. HEATMAP DE NES DOS TERMOS GSEA COMPARTILHADOS
# ==============================================================================

plot_shared_gsea <- function(
  gsea_result,
  database_name,
  top_n = 25
) {

  if (nrow(gsea_result$long) == 0) {
    message(
      database_name,
      ": nenhum termo GSEA significativo em pelo menos duas doenças."
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
        "GSEA compartilhado — ",
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
# 32. FUNÇÕES DO RESUMO GSEA
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
# 33. RESUMO FINAL DA INTEGRAÇÃO
# ==============================================================================

get_venn_region <- function(
  venn_result,
  pattern_name
) {

  value <- venn_result$regions[
    pattern == pattern_name,
    n_genes
  ]

  if (length(value) == 0) {
    return(0L)
  }

  as.integer(
    value[1]
  )
}

integration_summary <- data.table(
  metric = c(
    "common_testable_universe",
    "hd_pd_pairwise_universe",

    "common_primary_AD",
    "common_primary_HD",
    "common_primary_PD",

    "common_primary_AD_only",
    "common_primary_HD_only",
    "common_primary_PD_only",
    "common_primary_AD_HD",
    "common_primary_AD_PD",
    "common_primary_HD_PD",
    "common_primary_AD_HD_PD",

    "common_sensitivity_AD",
    "common_sensitivity_HD",
    "common_sensitivity_PD",
    "common_sensitivity_HD_PD",

    "article_like_AD",
    "article_like_HD",
    "article_like_PD",
    "article_like_AD_HD",
    "article_like_AD_PD",
    "article_like_HD_PD",
    "article_like_AD_HD_PD",

    "pairwise_HD_PD_shared_primary",
    "pairwise_HD_PD_shared_same_direction",
    "pairwise_HD_PD_shared_opposite_direction",

    "directionally_concordant_all_three",

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

    length(
      hd_pd_universe
    ),

    venn_common_primary$totals[
      dataset == "AD",
      n_genes
    ],

    venn_common_primary$totals[
      dataset == "HD",
      n_genes
    ],

    venn_common_primary$totals[
      dataset == "PD",
      n_genes
    ],

    get_venn_region(
      venn_common_primary,
      "AD only"
    ),

    get_venn_region(
      venn_common_primary,
      "HD only"
    ),

    get_venn_region(
      venn_common_primary,
      "PD only"
    ),

    get_venn_region(
      venn_common_primary,
      "AD+HD"
    ),

    get_venn_region(
      venn_common_primary,
      "AD+PD"
    ),

    get_venn_region(
      venn_common_primary,
      "HD+PD"
    ),

    get_venn_region(
      venn_common_primary,
      "AD+HD+PD"
    ),

    venn_common_sensitivity$totals[
      dataset == "AD",
      n_genes
    ],

    venn_common_sensitivity$totals[
      dataset == "HD",
      n_genes
    ],

    venn_common_sensitivity$totals[
      dataset == "PD",
      n_genes
    ],

    get_venn_region(
      venn_common_sensitivity,
      "HD+PD"
    ),

    venn_article_like$totals[
      dataset == "AD",
      n_genes
    ],

    venn_article_like$totals[
      dataset == "HD",
      n_genes
    ],

    venn_article_like$totals[
      dataset == "PD",
      n_genes
    ],

    get_venn_region(
      venn_article_like,
      "AD+HD"
    ),

    get_venn_region(
      venn_article_like,
      "AD+PD"
    ),

    get_venn_region(
      venn_article_like,
      "HD+PD"
    ),

    get_venn_region(
      venn_article_like,
      "AD+HD+PD"
    ),

    nrow(
      shared_hd_pd
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

# Novas métricas auditáveis, sem alterar as definicoes originais.
integration_summary <- rbind(
  integration_summary,
  data.table(
    metric = c(
      "article_like_n_genes_Alharbi_2025_reproduzidos",
      "article_like_n_genes_nossos_exclusivos_vs_Alharbi",
      "ORA_nominal_GO_BP_termos_testados",
      "ORA_nominal_GO_BP_termos_FDR005",
      "ORA_nominal_KEGG_termos_testados",
      "ORA_nominal_KEGG_termos_FDR005"
    ),
    value = c(
      sum(artigo_auditado$reproducao_nominal_tripla, na.rm=TRUE),
      sum(!genes_triplo_nominal$no_artigo_2025, na.rm=TRUE),
      nominal_ora_status[base=="GO_BP",n_termos_testados],
      nominal_ora_status[base=="GO_BP",n_termos_FDR005],
      nominal_ora_status[base=="KEGG",n_termos_testados],
      nominal_ora_status[base=="KEGG",n_termos_FDR005]
    )
  ),
  use.names = TRUE
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
# 34. RESUMO NO TERMINAL
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("RESUMO FINAL — INTEGRAÇÃO TRANSDIAGNÓSTICA\n")
cat("============================================================\n")

print(
  integration_summary
)

cat("\n")
cat("Correlações de Spearman — universos pareados:\n")

print(
  correlation_summary
)

cat("\n")
cat("HD vs PD — universo pareado:\n")

cat(
  "Genes: ",
  n_hd_pd_pair,
  "\n",
  sep = ""
)

cat(
  "Spearman rho: ",
  sprintf(
    "%.6f",
    rho_hd_pd
  ),
  "\n",
  sep = ""
)

cat(
  "DEGs compartilhados: ",
  nrow(
    shared_hd_pd
  ),
  "\n",
  sep = ""
)

cat(
  "Mesma direção: ",
  nrow(
    shared_hd_pd_same_direction
  ),
  "\n",
  sep = ""
)

cat(
  "Direção oposta: ",
  nrow(
    shared_hd_pd_opposite_direction
  ),
  "\n",
  sep = ""
)

cat("\n")
cat("Figuras geradas:\n")
cat(" - primary_deg_patterns.png\n")
cat(" - venn_common_primary_AD_HD_PD.png\n")
cat(" - venn_common_sensitivity_AD_HD_PD.png\n")
cat(" - venn_article_like_nominal_AD_HD_PD.png\n")
cat(" - log2FC_HD_vs_PD_pairwise.png\n")
cat(" - shared_GSEA_GO_BP.png\n")
cat(" - shared_GSEA_KEGG.png\n")

cat("\n")
cat("============================================================\n")
cat("07_cross_disease_integration.R CONCLUÍDO\n")
cat("============================================================\n")

cat("\nAUDITORIA COMPLEMENTAR FASE 01 CONCLUIDA\n")
cat("Interseção nominal tripla: ", length(ids_triplo_nominal), " genes\n", sep="")
cat("Genes de Alharbi 2025 reproduzidos: ",
    sum(artigo_auditado$reproducao_nominal_tripla, na.rm=TRUE), " / 10\n", sep="")
cat("ORA compartilhada HD-PD (GO BP / KEGG, All / Up / Down):\n")
print(ora_hd_pd_resumo)
cat("ORA nominal tripla — situacao por base:\n")
print(nominal_ora_status)
