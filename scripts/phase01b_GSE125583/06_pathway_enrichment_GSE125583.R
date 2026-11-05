#!/usr/bin/env Rscript

# ==============================================================================
# 06_pathway_enrichment_GSE125583.R
#
# FASE 01b — ENRIQUECIMENTO FUNCIONAL DO GSE125583
#
# OBJETIVO
# ------------------------------------------------------------------------------
# Investigar processos biológicos e vias associados à assinatura de expressão
# diferencial do GSE125583 usando duas abordagens complementares:
#
# 1) ORA (Over-Representation Analysis)
#    - usa os DEGs do critério principal;
#    - executada para All, Up e Down separadamente;
#    - bancos: GO Biological Process e KEGG.
#
# 2) GSEA (Gene Set Enrichment Analysis)
#    - usa TODOS os genes testados;
#    - ranking pela estatística de Wald ("stat") do DESeq2;
#    - bancos: GO Biological Process e KEGG.
#
# BACKGROUND DA ORA
# ------------------------------------------------------------------------------
# O universo é formado pelos genes realmente testados no DESeq2 após o
# pré-filtro da etapa 04. Genes que não foram testados não entram no background.
#
# CRITÉRIO DE INTERPRETAÇÃO
# ------------------------------------------------------------------------------
# BH-FDR (p.adjust) < 0,05.
#
# INTERPRETAÇÃO DA DIREÇÃO NA GSEA
# ------------------------------------------------------------------------------
# NES > 0 -> conjunto associado à extremidade positiva do ranking
#            (genes com maior expressão em AD).
#
# NES < 0 -> conjunto associado à extremidade negativa do ranking
#            (genes com menor expressão em AD).
#
# IMPORTANTE:
# nomes de vias KEGG que mencionam infecções/doenças não devem ser interpretados
# literalmente como evidência de que essas condições estão presentes. Essas vias
# podem compartilhar módulos imunes, inflamatórios ou de sinalização.
#
# ENTRADA
# ------------------------------------------------------------------------------
# results/candidate_GSE125583/deseq2/GSE125583_deseq2_results.tsv
#
# SAÍDAS
# ------------------------------------------------------------------------------
# results/candidate_GSE125583/enrichment/
#   ora_GO_BP_all.tsv
#   ora_GO_BP_up.tsv
#   ora_GO_BP_down.tsv
#   ora_KEGG_all.tsv
#   ora_KEGG_up.tsv
#   ora_KEGG_down.tsv
#   gsea_GO_BP.tsv
#   gsea_KEGG.tsv
#   enrichment_summary_GSE125583.tsv
#   06_pathway_enrichment_GSE125583.txt
#
# results/candidate_GSE125583/plots/enrichment/
#   ora_GO_BP_all.png
#   ora_GO_BP_up.png
#   ora_GO_BP_down.png
#   ora_KEGG_all.png
#   ora_KEGG_up.png
#   ora_KEGG_down.png
#   gsea_GO_BP.png
#   gsea_KEGG.png
#
# O QUE ESTE SCRIPT NÃO FAZ
# ------------------------------------------------------------------------------
# - não altera resultados do DESeq2;
# - não muda thresholds;
# - não usa o modelo ajustado como análise funcional principal;
# - não integra ainda AD com HD e PD.
# ==============================================================================


# ==============================================================================
# 1. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table",
  "clusterProfiler",
  "org.Hs.eg.db",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Pacotes necessários não instalados: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

deseq_results_file <- file.path(
  "results",
  "candidate_GSE125583",
  "deseq2",
  "GSE125583_deseq2_results.tsv"
)

output_dir <- file.path(
  "results",
  "candidate_GSE125583",
  "enrichment"
)

plot_dir <- file.path(
  "results",
  "candidate_GSE125583",
  "plots",
  "enrichment"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

audit_file <- file.path(
  output_dir,
  "06_pathway_enrichment_GSE125583.txt"
)

summary_file <- file.path(
  output_dir,
  "enrichment_summary_GSE125583.tsv"
)


# ==============================================================================
# 3. CONSTANTES ANALÍTICAS
# ==============================================================================

FDR_CUTOFF <- 0.05
N_TERMS_PLOT <- 15L

# Tamanho permitido dos conjuntos de genes na GSEA.
GSEA_MIN_GS_SIZE <- 10L
GSEA_MAX_GS_SIZE <- 500L


# ==============================================================================
# 4. FUNÇÕES AUXILIARES
# ==============================================================================

validate_file <- function(path) {

  if (!file.exists(path)) {
    stop("Arquivo não encontrado: ", path)
  }

  info <- file.info(path)

  if (is.na(info$size) || info$size == 0) {
    stop("Arquivo vazio ou tamanho indisponível: ", path)
  }

  invisible(TRUE)
}


format_int <- function(x) {

  format(
    x,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


# ------------------------------------------------------------------------------
# Salvar objeto enrichResult/gseaResult em TSV
# ------------------------------------------------------------------------------

save_enrichment_table <- function(object, output_file) {

  if (is.null(object)) {

    dt <- data.table()

    fwrite(
      dt,
      output_file,
      sep = "\t",
      quote = FALSE
    )

    return(dt)
  }

  dt <- as.data.table(
    as.data.frame(object)
  )

  fwrite(
    dt,
    output_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  dt
}


count_significant_terms <- function(dt) {

  if (
    nrow(dt) == 0 ||
    !"p.adjust" %in% names(dt)
  ) {
    return(0L)
  }

  sum(
    !is.na(dt$p.adjust) &
      dt$p.adjust < FDR_CUTOFF
  )
}


enrichment_theme <- function() {

  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      plot.subtitle = element_text(
        size = 10
      ),
      axis.title = element_text(
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 9
      ),
      legend.title = element_text(
        face = "bold"
      )
    )
}


# ------------------------------------------------------------------------------
# ORA — GO BP
# ------------------------------------------------------------------------------

run_enrich_go <- function(genes, universe) {

  if (length(genes) == 0) {
    return(NULL)
  }

  tryCatch(
    enrichGO(
      gene = genes,
      universe = universe,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = FALSE
    ),
    error = function(e) {

      warning(
        "Falha em enrichGO: ",
        conditionMessage(e)
      )

      NULL
    }
  )
}


# ------------------------------------------------------------------------------
# ORA — KEGG
# ------------------------------------------------------------------------------

run_enrich_kegg <- function(genes, universe) {

  if (length(genes) == 0) {
    return(NULL)
  }

  tryCatch(
    enrichKEGG(
      gene = genes,
      universe = universe,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    ),
    error = function(e) {

      warning(
        "Falha em enrichKEGG: ",
        conditionMessage(e)
      )

      NULL
    }
  )
}


# ------------------------------------------------------------------------------
# GSEA — GO BP
# ------------------------------------------------------------------------------

run_gsea_go <- function(gene_list) {

  tryCatch(
    gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      exponent = 1,
      minGSSize = GSEA_MIN_GS_SIZE,
      maxGSSize = GSEA_MAX_GS_SIZE,
      eps = 0,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = TRUE,
      by = "fgsea"
    ),
    error = function(e) {

      warning(
        "Falha em gseGO: ",
        conditionMessage(e)
      )

      NULL
    }
  )
}


# ------------------------------------------------------------------------------
# GSEA — KEGG
# ------------------------------------------------------------------------------

run_gsea_kegg <- function(gene_list) {

  tryCatch(
    gseKEGG(
      geneList = gene_list,
      organism = "hsa",
      keyType = "ncbi-geneid",
      exponent = 1,
      minGSSize = GSEA_MIN_GS_SIZE,
      maxGSSize = GSEA_MAX_GS_SIZE,
      eps = 0,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = TRUE,
      by = "fgsea"
    ),
    error = function(e) {

      warning(
        "Falha em gseKEGG: ",
        conditionMessage(e)
      )

      NULL
    }
  )
}


# ------------------------------------------------------------------------------
# Plot ORA
# ------------------------------------------------------------------------------

plot_ora <- function(
  ora_dt,
  title,
  subtitle,
  output_file
) {

  if (
    nrow(ora_dt) == 0 ||
    !"p.adjust" %in% names(ora_dt)
  ) {
    return(FALSE)
  }

  sig <- copy(
    ora_dt[
      !is.na(p.adjust) &
        p.adjust < FDR_CUTOFF
    ]
  )

  if (nrow(sig) == 0) {
    return(FALSE)
  }

  # GeneRatio vem como "número_de_hits/tamanho_da_lista".
  sig[
    ,
    GeneRatio_numeric := vapply(
      strsplit(
        as.character(GeneRatio),
        "/",
        fixed = TRUE
      ),
      function(x) {
        as.numeric(x[1]) /
          as.numeric(x[2])
      },
      numeric(1)
    )
  ]

  sig[
    ,
    minus_log10_fdr := -log10(
      pmax(
        p.adjust,
        .Machine$double.xmin
      )
    )
  ]

  setorder(
    sig,
    p.adjust,
    -Count
  )

  sig <- head(
    sig,
    N_TERMS_PLOT
  )

  sig[
    ,
    Description_plot := factor(
      Description,
      levels = rev(Description)
    )
  ]

  p <- ggplot(
    sig,
    aes(
      x = GeneRatio_numeric,
      y = Description_plot,
      size = Count,
      color = minus_log10_fdr
    )
  ) +
    geom_point(
      alpha = 0.85
    ) +
    scale_color_viridis_c(
      option = "viridis"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Gene ratio",
      y = NULL,
      size = "Genes",
      color = "-log10(FDR)"
    ) +
    enrichment_theme()

  ggsave(
    filename = output_file,
    plot = p,
    width = 10,
    height = 7,
    dpi = 300
  )

  TRUE
}


# ------------------------------------------------------------------------------
# Plot GSEA
# ------------------------------------------------------------------------------

plot_gsea <- function(
  gsea_dt,
  title,
  output_file
) {

  if (
    nrow(gsea_dt) == 0 ||
    !"p.adjust" %in% names(gsea_dt) ||
    !"NES" %in% names(gsea_dt)
  ) {
    return(FALSE)
  }

  sig <- copy(
    gsea_dt[
      !is.na(p.adjust) &
        p.adjust < FDR_CUTOFF &
        !is.na(NES)
    ]
  )

  if (nrow(sig) == 0) {
    return(FALSE)
  }

  positive <- sig[
    NES > 0
  ][
    order(
      p.adjust,
      -NES
    )
  ]

  negative <- sig[
    NES < 0
  ][
    order(
      p.adjust,
      NES
    )
  ]

  # Mostrar termos dos dois sentidos sempre que possível.
  n_each <- ceiling(
    N_TERMS_PLOT / 2
  )

  selected <- rbind(
    head(positive, n_each),
    head(negative, n_each),
    fill = TRUE
  )

  if (nrow(selected) > N_TERMS_PLOT) {

    selected <- selected[
      order(
        p.adjust,
        -abs(NES)
      )
    ][
      seq_len(N_TERMS_PLOT)
    ]
  }

  selected[
    ,
    Direction := fifelse(
      NES > 0,
      "AD-up",
      "AD-down"
    )
  ]

  selected[
    ,
    minus_log10_fdr := -log10(
      pmax(
        p.adjust,
        .Machine$double.xmin
      )
    )
  ]

  setorder(
    selected,
    NES
  )

  selected[
    ,
    Description_plot := factor(
      Description,
      levels = Description
    )
  ]

  p <- ggplot(
    selected,
    aes(
      x = NES,
      y = Description_plot,
      size = minus_log10_fdr,
      color = Direction
    )
  ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.4
    ) +
    geom_point(
      alpha = 0.85
    ) +
    scale_color_manual(
      values = c(
        "AD-up" = "#B2182B",
        "AD-down" = "#2166AC"
      )
    ) +
    labs(
      title = title,
      subtitle = "BH-FDR < 0,05; ranking pela estatística de Wald",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      size = "-log10(FDR)",
      color = "Direção"
    ) +
    enrichment_theme()

  ggsave(
    filename = output_file,
    plot = p,
    width = 11,
    height = 7.5,
    dpi = 300
  )

  TRUE
}


# ==============================================================================
# 5. LER E VALIDAR RESULTADOS DO DESeq2
# ==============================================================================

validate_file(
  deseq_results_file
)

results_dt <- fread(
  deseq_results_file,
  na.strings = c(
    "NA",
    ""
  ),
  showProgress = FALSE
)

required_columns <- c(
  "GeneID",
  "Symbol",
  "stat",
  "pvalue",
  "padj",
  "log2FoldChange",
  "is_primary",
  "primary_direction"
)

missing_columns <- setdiff(
  required_columns,
  names(results_dt)
)

if (length(missing_columns) > 0) {
  stop(
    "Arquivo do DESeq2 sem colunas obrigatórias: ",
    paste(missing_columns, collapse = ", ")
  )
}

results_dt[
  ,
  GeneID := as.character(GeneID)
]

if (
  any(
    is.na(results_dt$GeneID) |
      trimws(results_dt$GeneID) == ""
  )
) {
  stop(
    "Existem GeneIDs ausentes nos resultados do DESeq2."
  )
}

if (anyDuplicated(results_dt$GeneID)) {
  stop(
    "Existem GeneIDs duplicados nos resultados do DESeq2."
  )
}


# ==============================================================================
# 6. DEFINIR UNIVERSO E LISTAS DA ORA
# ==============================================================================

# Background = todos os genes testados no DESeq2.
universe_genes <- unique(
  results_dt$GeneID
)

primary_all <- unique(
  results_dt[
    is_primary == TRUE,
    GeneID
  ]
)

primary_up <- unique(
  results_dt[
    primary_direction == "Up",
    GeneID
  ]
)

primary_down <- unique(
  results_dt[
    primary_direction == "Down",
    GeneID
  ]
)

if (!all(primary_all %in% universe_genes)) {
  stop(
    "Existem DEGs principais fora do universo testado."
  )
}

if (!all(primary_up %in% primary_all)) {
  stop(
    "Existem genes Up fora dos DEGs principais."
  )
}

if (!all(primary_down %in% primary_all)) {
  stop(
    "Existem genes Down fora dos DEGs principais."
  )
}

if (length(intersect(primary_up, primary_down)) > 0) {
  stop(
    "Há genes simultaneamente classificados como Up e Down."
  )
}


# ==============================================================================
# 7. ORA — GO BIOLOGICAL PROCESS
# ==============================================================================

ora_go_all <- run_enrich_go(
  primary_all,
  universe_genes
)

ora_go_up <- run_enrich_go(
  primary_up,
  universe_genes
)

ora_go_down <- run_enrich_go(
  primary_down,
  universe_genes
)


# ==============================================================================
# 8. ORA — KEGG
# ==============================================================================

ora_kegg_all <- run_enrich_kegg(
  primary_all,
  universe_genes
)

ora_kegg_up <- run_enrich_kegg(
  primary_up,
  universe_genes
)

ora_kegg_down <- run_enrich_kegg(
  primary_down,
  universe_genes
)


# ==============================================================================
# 9. SALVAR RESULTADOS DA ORA
# ==============================================================================

ora_go_all_dt <- save_enrichment_table(
  ora_go_all,
  file.path(
    output_dir,
    "ora_GO_BP_all.tsv"
  )
)

ora_go_up_dt <- save_enrichment_table(
  ora_go_up,
  file.path(
    output_dir,
    "ora_GO_BP_up.tsv"
  )
)

ora_go_down_dt <- save_enrichment_table(
  ora_go_down,
  file.path(
    output_dir,
    "ora_GO_BP_down.tsv"
  )
)

ora_kegg_all_dt <- save_enrichment_table(
  ora_kegg_all,
  file.path(
    output_dir,
    "ora_KEGG_all.tsv"
  )
)

ora_kegg_up_dt <- save_enrichment_table(
  ora_kegg_up,
  file.path(
    output_dir,
    "ora_KEGG_up.tsv"
  )
)

ora_kegg_down_dt <- save_enrichment_table(
  ora_kegg_down,
  file.path(
    output_dir,
    "ora_KEGG_down.tsv"
  )
)


# ==============================================================================
# 10. CONSTRUIR RANKING DA GSEA
# ==============================================================================

# A estatística de Wald preserva magnitude e direção do efeito e é adequada
# para ordenar os genes na GSEA.
gsea_source <- results_dt[
  !is.na(stat) &
    is.finite(stat)
]

if (nrow(gsea_source) == 0) {
  stop(
    "Nenhum gene com estatística de Wald válida para GSEA."
  )
}

if (anyDuplicated(gsea_source$GeneID)) {
  stop(
    "GeneIDs duplicados no ranking da GSEA."
  )
}

gene_list <- gsea_source$stat

names(gene_list) <- gsea_source$GeneID

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)


# ==============================================================================
# 11. GSEA — GO BP E KEGG
# ==============================================================================

gsea_go <- run_gsea_go(
  gene_list
)

gsea_kegg <- run_gsea_kegg(
  gene_list
)

gsea_go_dt <- save_enrichment_table(
  gsea_go,
  file.path(
    output_dir,
    "gsea_GO_BP.tsv"
  )
)

gsea_kegg_dt <- save_enrichment_table(
  gsea_kegg,
  file.path(
    output_dir,
    "gsea_KEGG.tsv"
  )
)


# ==============================================================================
# 12. GERAR FIGURAS DA ORA
# ==============================================================================

plot_ora_go_all <- plot_ora(
  ora_go_all_dt,
  "GO BP ORA — AD (All DEGs)",
  "DEGs principais; background = genes testados no DESeq2",
  file.path(
    plot_dir,
    "ora_GO_BP_all.png"
  )
)

plot_ora_go_up <- plot_ora(
  ora_go_up_dt,
  "GO BP ORA — AD (Up)",
  "DEGs principais com maior expressão em AD",
  file.path(
    plot_dir,
    "ora_GO_BP_up.png"
  )
)

plot_ora_go_down <- plot_ora(
  ora_go_down_dt,
  "GO BP ORA — AD (Down)",
  "DEGs principais com menor expressão em AD",
  file.path(
    plot_dir,
    "ora_GO_BP_down.png"
  )
)

plot_ora_kegg_all <- plot_ora(
  ora_kegg_all_dt,
  "KEGG ORA — AD (All DEGs)",
  "DEGs principais; background = genes testados no DESeq2",
  file.path(
    plot_dir,
    "ora_KEGG_all.png"
  )
)

plot_ora_kegg_up <- plot_ora(
  ora_kegg_up_dt,
  "KEGG ORA — AD (Up)",
  "DEGs principais com maior expressão em AD",
  file.path(
    plot_dir,
    "ora_KEGG_up.png"
  )
)

plot_ora_kegg_down <- plot_ora(
  ora_kegg_down_dt,
  "KEGG ORA — AD (Down)",
  "DEGs principais com menor expressão em AD",
  file.path(
    plot_dir,
    "ora_KEGG_down.png"
  )
)


# ==============================================================================
# 13. GERAR FIGURAS DA GSEA
# ==============================================================================

plot_gsea_go <- plot_gsea(
  gsea_go_dt,
  "GSEA — GO Biological Process — AD",
  file.path(
    plot_dir,
    "gsea_GO_BP.png"
  )
)

plot_gsea_kegg <- plot_gsea(
  gsea_kegg_dt,
  "GSEA — KEGG — AD",
  file.path(
    plot_dir,
    "gsea_KEGG.png"
  )
)


# ==============================================================================
# 14. RESUMO QUANTITATIVO
# ==============================================================================

summary_dt <- data.table(
  analysis = c(
    "ORA_GO_BP_all",
    "ORA_GO_BP_up",
    "ORA_GO_BP_down",
    "ORA_KEGG_all",
    "ORA_KEGG_up",
    "ORA_KEGG_down",
    "GSEA_GO_BP",
    "GSEA_KEGG"
  ),

  n_input_genes = c(
    length(primary_all),
    length(primary_up),
    length(primary_down),
    length(primary_all),
    length(primary_up),
    length(primary_down),
    length(gene_list),
    length(gene_list)
  ),

  n_terms_total = c(
    nrow(ora_go_all_dt),
    nrow(ora_go_up_dt),
    nrow(ora_go_down_dt),
    nrow(ora_kegg_all_dt),
    nrow(ora_kegg_up_dt),
    nrow(ora_kegg_down_dt),
    nrow(gsea_go_dt),
    nrow(gsea_kegg_dt)
  ),

  n_terms_FDR_lt_005 = c(
    count_significant_terms(ora_go_all_dt),
    count_significant_terms(ora_go_up_dt),
    count_significant_terms(ora_go_down_dt),
    count_significant_terms(ora_kegg_all_dt),
    count_significant_terms(ora_kegg_up_dt),
    count_significant_terms(ora_kegg_down_dt),
    count_significant_terms(gsea_go_dt),
    count_significant_terms(gsea_kegg_dt)
  ),

  plot_created = c(
    plot_ora_go_all,
    plot_ora_go_up,
    plot_ora_go_down,
    plot_ora_kegg_all,
    plot_ora_kegg_up,
    plot_ora_kegg_down,
    plot_gsea_go,
    plot_gsea_kegg
  )
)

fwrite(
  summary_dt,
  summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 15. CONTAR DIREÇÃO DOS TERMOS DA GSEA
# ==============================================================================

count_gsea_direction <- function(dt, direction = c("positive", "negative")) {

  direction <- match.arg(direction)

  if (
    nrow(dt) == 0 ||
    !all(
      c(
        "p.adjust",
        "NES"
      ) %in% names(dt)
    )
  ) {
    return(0L)
  }

  if (direction == "positive") {

    return(
      sum(
        dt$p.adjust < FDR_CUTOFF &
          dt$NES > 0,
        na.rm = TRUE
      )
    )
  }

  sum(
    dt$p.adjust < FDR_CUTOFF &
      dt$NES < 0,
    na.rm = TRUE
  )
}

gsea_go_positive <- count_gsea_direction(
  gsea_go_dt,
  "positive"
)

gsea_go_negative <- count_gsea_direction(
  gsea_go_dt,
  "negative"
)

gsea_kegg_positive <- count_gsea_direction(
  gsea_kegg_dt,
  "positive"
)

gsea_kegg_negative <- count_gsea_direction(
  gsea_kegg_dt,
  "negative"
)


# ==============================================================================
# 16. RELATÓRIO TEXTUAL
# ==============================================================================

sink(
  audit_file
)

cat("============================================================\n")
cat("06 — ENRIQUECIMENTO FUNCIONAL GSE125583\n")
cat("============================================================\n\n")


cat("------------------------------------------------------------\n")
cat("1. ENTRADA\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes testados / background ORA: ",
  format_int(length(universe_genes)),
  "\n",
  sep = ""
)

cat(
  "DEGs principais — All: ",
  format_int(length(primary_all)),
  "\n",
  sep = ""
)

cat(
  "DEGs principais — Up: ",
  format_int(length(primary_up)),
  "\n",
  sep = ""
)

cat(
  "DEGs principais — Down: ",
  format_int(length(primary_down)),
  "\n",
  sep = ""
)

cat(
  "Genes no ranking GSEA: ",
  format_int(length(gene_list)),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("2. ORA — GO BIOLOGICAL PROCESS\n")
cat("------------------------------------------------------------\n")

cat(
  "All — termos totais: ",
  format_int(nrow(ora_go_all_dt)),
  "\n",
  sep = ""
)

cat(
  "All — FDR < 0,05: ",
  format_int(count_significant_terms(ora_go_all_dt)),
  "\n",
  sep = ""
)

cat(
  "Up — termos totais: ",
  format_int(nrow(ora_go_up_dt)),
  "\n",
  sep = ""
)

cat(
  "Up — FDR < 0,05: ",
  format_int(count_significant_terms(ora_go_up_dt)),
  "\n",
  sep = ""
)

cat(
  "Down — termos totais: ",
  format_int(nrow(ora_go_down_dt)),
  "\n",
  sep = ""
)

cat(
  "Down — FDR < 0,05: ",
  format_int(count_significant_terms(ora_go_down_dt)),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("3. ORA — KEGG\n")
cat("------------------------------------------------------------\n")

cat(
  "All — termos totais: ",
  format_int(nrow(ora_kegg_all_dt)),
  "\n",
  sep = ""
)

cat(
  "All — FDR < 0,05: ",
  format_int(count_significant_terms(ora_kegg_all_dt)),
  "\n",
  sep = ""
)

cat(
  "Up — termos totais: ",
  format_int(nrow(ora_kegg_up_dt)),
  "\n",
  sep = ""
)

cat(
  "Up — FDR < 0,05: ",
  format_int(count_significant_terms(ora_kegg_up_dt)),
  "\n",
  sep = ""
)

cat(
  "Down — termos totais: ",
  format_int(nrow(ora_kegg_down_dt)),
  "\n",
  sep = ""
)

cat(
  "Down — FDR < 0,05: ",
  format_int(count_significant_terms(ora_kegg_down_dt)),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("4. GSEA — GO BIOLOGICAL PROCESS\n")
cat("------------------------------------------------------------\n")

cat(
  "Termos totais: ",
  format_int(nrow(gsea_go_dt)),
  "\n",
  sep = ""
)

cat(
  "FDR < 0,05: ",
  format_int(count_significant_terms(gsea_go_dt)),
  "\n",
  sep = ""
)

cat(
  "NES > 0 e FDR < 0,05: ",
  format_int(gsea_go_positive),
  "\n",
  sep = ""
)

cat(
  "NES < 0 e FDR < 0,05: ",
  format_int(gsea_go_negative),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("5. GSEA — KEGG\n")
cat("------------------------------------------------------------\n")

cat(
  "Termos totais: ",
  format_int(nrow(gsea_kegg_dt)),
  "\n",
  sep = ""
)

cat(
  "FDR < 0,05: ",
  format_int(count_significant_terms(gsea_kegg_dt)),
  "\n",
  sep = ""
)

cat(
  "NES > 0 e FDR < 0,05: ",
  format_int(gsea_kegg_positive),
  "\n",
  sep = ""
)

cat(
  "NES < 0 e FDR < 0,05: ",
  format_int(gsea_kegg_negative),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("6. INTERPRETAÇÃO METODOLÓGICA\n")
cat("------------------------------------------------------------\n")

cat(
  "- ORA usa os DEGs principais e os genes testados no DESeq2 como background.\n"
)

cat(
  "- ORA foi executada para All, Up e Down separadamente.\n"
)

cat(
  "- GSEA usa todos os genes com estatística de Wald válida.\n"
)

cat(
  "- NES positivo indica enriquecimento associado à extremidade AD-up; NES negativo indica enriquecimento associado à extremidade AD-down.\n"
)

cat(
  "- A significância funcional é interpretada com BH-FDR < 0,05.\n"
)

cat(
  "- ORA e GSEA são análises complementares e não precisam produzir conjuntos idênticos de termos.\n"
)

cat(
  "- Nomes de vias KEGG relacionados a doenças ou infecções não devem ser interpretados literalmente como diagnóstico clínico.\n"
)

cat(
  "- A integração funcional com HD e PD será realizada em etapa posterior.\n"
)

cat("\n")
cat("============================================================\n")
cat("06_pathway_enrichment_GSE125583.R CONCLUÍDO\n")
cat("============================================================\n")

sink()


# ==============================================================================
# 17. RESUMO NO TERMINAL
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("06 — ENRIQUECIMENTO GSE125583 CONCLUÍDO\n")
cat("============================================================\n")

cat(
  "Background ORA: ",
  format_int(length(universe_genes)),
  " genes\n",
  sep = ""
)

cat(
  "DEGs principais: ",
  format_int(length(primary_all)),
  " (",
  format_int(length(primary_up)),
  " Up / ",
  format_int(length(primary_down)),
  " Down)\n",
  sep = ""
)

cat(
  "ORA GO BP significativos — All/Up/Down: ",
  count_significant_terms(ora_go_all_dt),
  "/",
  count_significant_terms(ora_go_up_dt),
  "/",
  count_significant_terms(ora_go_down_dt),
  "\n",
  sep = ""
)

cat(
  "ORA KEGG significativos — All/Up/Down: ",
  count_significant_terms(ora_kegg_all_dt),
  "/",
  count_significant_terms(ora_kegg_up_dt),
  "/",
  count_significant_terms(ora_kegg_down_dt),
  "\n",
  sep = ""
)

cat(
  "GSEA GO BP significativos: ",
  count_significant_terms(gsea_go_dt),
  "\n",
  sep = ""
)

cat(
  "GSEA KEGG significativos: ",
  count_significant_terms(gsea_kegg_dt),
  "\n",
  sep = ""
)

cat(
  "Resumo: ",
  summary_file,
  "\n",
  sep = ""
)

cat(
  "Relatório: ",
  audit_file,
  "\n",
  sep = ""
)

cat("============================================================\n")
