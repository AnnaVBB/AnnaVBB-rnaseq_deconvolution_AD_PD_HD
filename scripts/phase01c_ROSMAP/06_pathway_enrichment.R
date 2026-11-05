# ==============================================================================
# 06_pathway_enrichment_ROSMAP.R
# ==============================================================================
# Fase 01c — ROSMAP
# Enriquecimento funcional dos modelos principal e ajustado
#
# Entradas:
#   results/candidate_ROSMAP/deseq2/ROSMAP_deseq2_results.tsv
#   results/candidate_ROSMAP/deseq2/ROSMAP_deseq2_results_adjusted.tsv
#
# Estratégia:
#   1) ORA: GO Biological Process e KEGG
#      - alvo: DEGs do critério principal
#      - universo: TODOS os genes testados pelo DESeq2
#      - análise separada: all / up / down
#
#   2) GSEA: GO Biological Process e KEGG
#      - usa TODOS os genes com Wald statistic disponível
#      - ranking = stat do DESeq2
#      - NES > 0: enriquecimento na direção AD
#      - NES < 0: enriquecimento na direção Control / reduzido em AD
#
#   3) Comparação principal × ajustado
#      - termos significativos
#      - interseção
#      - Jaccard
#      - concordância da direção do NES
#
# Critério funcional:
#   BH-FDR / p.adjust < 0,05
#
# IMPORTANTE:
# O objetivo desta etapa não é escolher o modelo com "mais vias", e sim testar
# se a convergência funcional é mais estável que a classificação gene-a-gene.
#
# Como algumas listas ROSMAP são muito pequenas, ORA com <5 genes é pulada.
# A GSEA é especialmente importante porque não depende de cutoff binário.
# ==============================================================================


# ==============================================================================
# 0. PACOTES
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
# 1. PARÂMETROS
# ==============================================================================

FUNCTIONAL_ALPHA <- 0.05
MIN_ORA_GENES <- 5L

GSEA_MIN_GS_SIZE <- 10L
GSEA_MAX_GS_SIZE <- 500L

TOP_PLOT_TERMS <- 15L


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

candidate_dir <- file.path("results", "candidate_ROSMAP")
deseq2_dir <- file.path(candidate_dir, "deseq2")

enrichment_dir <- file.path(candidate_dir, "enrichment")
ora_dir <- file.path(enrichment_dir, "ORA")
gsea_dir <- file.path(enrichment_dir, "GSEA")
comparison_dir <- file.path(enrichment_dir, "model_comparison")
plot_dir <- file.path(enrichment_dir, "plots")

for (d in c(enrichment_dir, ora_dir, gsea_dir, comparison_dir, plot_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

main_results_input <- file.path(
  deseq2_dir,
  "ROSMAP_deseq2_results.tsv"
)

adjusted_results_input <- file.path(
  deseq2_dir,
  "ROSMAP_deseq2_results_adjusted.tsv"
)

summary_output <- file.path(
  enrichment_dir,
  "ROSMAP_enrichment_summary.tsv"
)

ora_status_output <- file.path(
  ora_dir,
  "ROSMAP_ORA_analysis_status.tsv"
)

gsea_status_output <- file.path(
  gsea_dir,
  "ROSMAP_GSEA_analysis_status.tsv"
)

comparison_summary_output <- file.path(
  comparison_dir,
  "ROSMAP_functional_stability_summary.tsv"
)

comparison_terms_output <- file.path(
  comparison_dir,
  "ROSMAP_functional_stability_terms.tsv"
)

gsea_stability_plot_output <- file.path(
  plot_dir,
  "GSEA_NES_main_vs_adjusted_ROSMAP.png"
)

report_output <- file.path(
  enrichment_dir,
  "06_pathway_enrichment_ROSMAP.txt"
)

session_info_output <- file.path(
  enrichment_dir,
  "ROSMAP_enrichment_sessionInfo.txt"
)


# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================

assert_file <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo não encontrado:\n", path)
  }
  invisible(TRUE)
}

assert_columns <- function(dt, required, object_name) {
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop(
      object_name,
      ": colunas obrigatórias ausentes: ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}

fmt_int <- function(x) {
  format(
    round(as.numeric(x), 0),
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

fmt_num <- function(x, digits = 4) {
  format(
    round(as.numeric(x), digits),
    nsmall = digits,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

safe_jaccard <- function(x, y) {
  x <- unique(x[!is.na(x)])
  y <- unique(y[!is.na(y)])
  u <- union(x, y)
  if (length(u) == 0) return(NA_real_)
  length(intersect(x, y)) / length(u)
}

cp_to_dt <- function(x) {
  if (is.null(x)) return(data.table())
  out <- tryCatch(as.data.frame(x), error = function(e) NULL)
  if (is.null(out) || nrow(out) == 0) return(data.table())
  as.data.table(out)
}

write_empty_or_table <- function(dt, path, columns = character()) {
  if (nrow(dt) == 0 && length(names(dt)) == 0 && length(columns) > 0) {
    dt <- as.data.table(
      setNames(
        replicate(length(columns), character(), simplify = FALSE),
        columns
      )
    )
  }

  fwrite(
    dt,
    path,
    sep = "\t",
    na = "NA"
  )
}

save_empty_plot <- function(title, subtitle, output_file) {
  p <- ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "Nenhum termo significativo\n(BH-FDR < 0,05)",
      size = 5
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(
      title = title,
      subtitle = subtitle
    ) +
    theme_void()

  ggsave(
    output_file,
    p,
    width = 8,
    height = 5,
    dpi = 300
  )
}

parse_ratio <- function(x) {
  vapply(
    x,
    function(z) {
      z <- strsplit(as.character(z), "/", fixed = TRUE)[[1]]
      if (length(z) != 2) return(NA_real_)
      as.numeric(z[1]) / as.numeric(z[2])
    },
    numeric(1)
  )
}

plot_ora <- function(dt, title, output_file) {

  sig <- dt[
    !is.na(p.adjust) &
      p.adjust < FUNCTIONAL_ALPHA
  ]

  if (nrow(sig) == 0) {
    save_empty_plot(title, "ORA", output_file)
    return(invisible(NULL))
  }

  sig[, GeneRatio_numeric := parse_ratio(GeneRatio)]
  sig <- sig[
    order(p.adjust, -GeneRatio_numeric)
  ][
    seq_len(min(.N, TOP_PLOT_TERMS))
  ]

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
      color = -log10(p.adjust)
    )
  ) +
    geom_point(alpha = 0.85) +
    labs(
      title = title,
      x = "Gene ratio",
      y = NULL,
      size = "Genes",
      color = "-log10(FDR)"
    ) +
    theme_bw(base_size = 11)

  ggsave(
    output_file,
    p,
    width = 10,
    height = 7,
    dpi = 300
  )
}

plot_gsea <- function(dt, title, output_file) {

  sig <- dt[
    !is.na(p.adjust) &
      p.adjust < FUNCTIONAL_ALPHA &
      is.finite(NES)
  ]

  if (nrow(sig) == 0) {
    save_empty_plot(title, "GSEA", output_file)
    return(invisible(NULL))
  }

  sig[, abs_NES := abs(NES)]

  sig <- sig[
    order(p.adjust, -abs_NES)
  ][
    seq_len(min(.N, TOP_PLOT_TERMS))
  ]

  sig[, direction := fifelse(NES > 0, "AD", "Control")]

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
      x = NES,
      y = Description_plot,
      fill = direction
    )
  ) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    labs(
      title = title,
      subtitle = "NES > 0: direção AD | NES < 0: direção Control",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      fill = "Direção"
    ) +
    theme_bw(base_size = 11)

  ggsave(
    output_file,
    p,
    width = 10,
    height = 7,
    dpi = 300
  )
}


# ==============================================================================
# 4. LER E VALIDAR RESULTADOS DESeq2
# ==============================================================================

assert_file(main_results_input)
assert_file(adjusted_results_input)

main <- fread(
  main_results_input,
  na.strings = c("NA", "")
)

adjusted <- fread(
  adjusted_results_input,
  na.strings = c("NA", "")
)

required_result_cols <- c(
  "GeneID",
  "log2FoldChange",
  "stat",
  "pvalue",
  "padj",
  "DEG",
  "direction"
)

assert_columns(
  main,
  required_result_cols,
  "ROSMAP_deseq2_results.tsv"
)

assert_columns(
  adjusted,
  required_result_cols,
  "ROSMAP_deseq2_results_adjusted.tsv"
)

main[, GeneID := as.character(GeneID)]
adjusted[, GeneID := as.character(GeneID)]

if (
  anyDuplicated(main$GeneID) ||
  anyDuplicated(adjusted$GeneID)
) {
  stop("GeneID duplicado nos resultados do DESeq2.")
}

if (!setequal(main$GeneID, adjusted$GeneID)) {
  stop(
    "Modelo principal e ajustado não contêm o mesmo universo de genes."
  )
}

adjusted <- adjusted[
  match(main$GeneID, GeneID)
]

if (anyNA(adjusted$GeneID)) {
  stop("Falha ao alinhar resultados principal e ajustado.")
}

universe_geneids <- unique(main$GeneID)
n_universe <- length(universe_geneids)


# ==============================================================================
# 5. CONJUNTOS DE DEGs PARA ORA
# ==============================================================================

get_deg_sets <- function(res_dt) {
  list(
    all = unique(
      res_dt[
        DEG == TRUE,
        GeneID
      ]
    ),
    up = unique(
      res_dt[
        DEG == TRUE &
          direction == "Up",
        GeneID
      ]
    ),
    down = unique(
      res_dt[
        DEG == TRUE &
          direction == "Down",
        GeneID
      ]
    )
  )
}

main_sets <- get_deg_sets(main)
adjusted_sets <- get_deg_sets(adjusted)


# ==============================================================================
# 6. ORA — FUNÇÕES
# ==============================================================================

run_ora_go <- function(
  genes,
  universe,
  model_name,
  direction_name
) {

  output_file <- file.path(
    ora_dir,
    paste0(
      "ORA_GO_BP_",
      model_name,
      "_",
      direction_name,
      ".tsv"
    )
  )

  plot_file <- file.path(
    plot_dir,
    paste0(
      "ORA_GO_BP_",
      model_name,
      "_",
      direction_name,
      ".png"
    )
  )

  if (length(genes) < MIN_ORA_GENES) {

    write_empty_or_table(
      data.table(),
      output_file,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "pvalue",
        "p.adjust",
        "qvalue",
        "geneID",
        "Count"
      )
    )

    save_empty_plot(
      paste0(
        "ROSMAP — GO:BP ORA | ",
        model_name,
        " | ",
        direction_name
      ),
      paste0(
        "Análise não executada: ",
        length(genes),
        " genes (< ",
        MIN_ORA_GENES,
        ")"
      ),
      plot_file
    )

    return(
      list(
        result = data.table(),
        status = data.table(
          analysis_type = "ORA",
          database = "GO_BP",
          model = model_name,
          direction = direction_name,
          n_input_genes = length(genes),
          status = "skipped_too_few_genes",
          message = paste0("n < ", MIN_ORA_GENES)
        )
      )
    )
  }

  obj <- tryCatch(
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
    error = function(e) e
  )

  if (inherits(obj, "error")) {

    write_empty_or_table(
      data.table(),
      output_file,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "pvalue",
        "p.adjust",
        "qvalue",
        "geneID",
        "Count"
      )
    )

    save_empty_plot(
      paste0(
        "ROSMAP — GO:BP ORA | ",
        model_name,
        " | ",
        direction_name
      ),
      paste0("Erro: ", conditionMessage(obj)),
      plot_file
    )

    return(
      list(
        result = data.table(),
        status = data.table(
          analysis_type = "ORA",
          database = "GO_BP",
          model = model_name,
          direction = direction_name,
          n_input_genes = length(genes),
          status = "error",
          message = conditionMessage(obj)
        )
      )
    )
  }

  dt <- cp_to_dt(obj)

  write_empty_or_table(
    dt,
    output_file
  )

  plot_ora(
    dt,
    paste0(
      "ROSMAP — GO:BP ORA | ",
      model_name,
      " | ",
      direction_name
    ),
    plot_file
  )

  list(
    result = dt,
    status = data.table(
      analysis_type = "ORA",
      database = "GO_BP",
      model = model_name,
      direction = direction_name,
      n_input_genes = length(genes),
      status = "completed",
      message = NA_character_
    )
  )
}


run_ora_kegg <- function(
  genes,
  universe,
  model_name,
  direction_name
) {

  output_file <- file.path(
    ora_dir,
    paste0(
      "ORA_KEGG_",
      model_name,
      "_",
      direction_name,
      ".tsv"
    )
  )

  plot_file <- file.path(
    plot_dir,
    paste0(
      "ORA_KEGG_",
      model_name,
      "_",
      direction_name,
      ".png"
    )
  )

  if (length(genes) < MIN_ORA_GENES) {

    write_empty_or_table(
      data.table(),
      output_file,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "pvalue",
        "p.adjust",
        "qvalue",
        "geneID",
        "Count"
      )
    )

    save_empty_plot(
      paste0(
        "ROSMAP — KEGG ORA | ",
        model_name,
        " | ",
        direction_name
      ),
      paste0(
        "Análise não executada: ",
        length(genes),
        " genes (< ",
        MIN_ORA_GENES,
        ")"
      ),
      plot_file
    )

    return(
      list(
        result = data.table(),
        status = data.table(
          analysis_type = "ORA",
          database = "KEGG",
          model = model_name,
          direction = direction_name,
          n_input_genes = length(genes),
          status = "skipped_too_few_genes",
          message = paste0("n < ", MIN_ORA_GENES)
        )
      )
    )
  }

  obj <- tryCatch(
    enrichKEGG(
      gene = genes,
      universe = universe,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = GSEA_MIN_GS_SIZE,
      maxGSSize = GSEA_MAX_GS_SIZE
    ),
    error = function(e) e
  )

  if (inherits(obj, "error")) {

    write_empty_or_table(
      data.table(),
      output_file,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "pvalue",
        "p.adjust",
        "qvalue",
        "geneID",
        "Count"
      )
    )

    save_empty_plot(
      paste0(
        "ROSMAP — KEGG ORA | ",
        model_name,
        " | ",
        direction_name
      ),
      paste0("Erro: ", conditionMessage(obj)),
      plot_file
    )

    return(
      list(
        result = data.table(),
        status = data.table(
          analysis_type = "ORA",
          database = "KEGG",
          model = model_name,
          direction = direction_name,
          n_input_genes = length(genes),
          status = "error",
          message = conditionMessage(obj)
        )
      )
    )
  }

  dt <- cp_to_dt(obj)

  write_empty_or_table(
    dt,
    output_file
  )

  plot_ora(
    dt,
    paste0(
      "ROSMAP — KEGG ORA | ",
      model_name,
      " | ",
      direction_name
    ),
    plot_file
  )

  list(
    result = dt,
    status = data.table(
      analysis_type = "ORA",
      database = "KEGG",
      model = model_name,
      direction = direction_name,
      n_input_genes = length(genes),
      status = "completed",
      message = NA_character_
    )
  )
}


# ==============================================================================
# 7. EXECUTAR ORA
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "06 — ENRIQUECIMENTO FUNCIONAL ROSMAP\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Universo funcional: ",
  fmt_int(n_universe),
  " genes.\n\n",
  sep = ""
)

ora_results <- list()
ora_status <- list()

models_sets <- list(
  main = main_sets,
  adjusted = adjusted_sets
)

for (model_name in names(models_sets)) {

  set_list <- models_sets[[model_name]]

  for (direction_name in names(set_list)) {

    genes <- set_list[[direction_name]]

    cat(
      "ORA ",
      model_name,
      " / ",
      direction_name,
      ": ",
      length(genes),
      " genes.\n",
      sep = ""
    )

    go_run <- run_ora_go(
      genes,
      universe_geneids,
      model_name,
      direction_name
    )

    kegg_run <- run_ora_kegg(
      genes,
      universe_geneids,
      model_name,
      direction_name
    )

    ora_results[[paste(
      model_name,
      direction_name,
      "GO_BP",
      sep = "__"
    )]] <- go_run$result

    ora_results[[paste(
      model_name,
      direction_name,
      "KEGG",
      sep = "__"
    )]] <- kegg_run$result

    ora_status[[length(ora_status) + 1L]] <- go_run$status
    ora_status[[length(ora_status) + 1L]] <- kegg_run$status
  }
}

ora_status_dt <- rbindlist(
  ora_status,
  fill = TRUE
)

fwrite(
  ora_status_dt,
  ora_status_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 8. GSEA — RANKINGS
# ==============================================================================

make_ranked_list <- function(res_dt) {

  ranked <- res_dt[
    !is.na(stat) &
      is.finite(stat),
    .(
      GeneID,
      stat
    )
  ]

  if (anyDuplicated(ranked$GeneID)) {
    stop("GeneID duplicado ao preparar ranking GSEA.")
  }

  ranked <- ranked[
    order(-stat)
  ]

  gene_list <- ranked$stat
  names(gene_list) <- ranked$GeneID

  gene_list
}

main_ranked <- make_ranked_list(main)
adjusted_ranked <- make_ranked_list(adjusted)


# ==============================================================================
# 9. GSEA — FUNÇÕES
# ==============================================================================

run_gsea_go <- function(
  gene_list,
  model_name
) {

  output_file <- file.path(
    gsea_dir,
    paste0(
      "GSEA_GO_BP_",
      model_name,
      ".tsv"
    )
  )

  plot_file <- file.path(
    plot_dir,
    paste0(
      "GSEA_GO_BP_",
      model_name,
      ".png"
    )
  )

  obj <- tryCatch(
    gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = GSEA_MIN_GS_SIZE,
      maxGSSize = GSEA_MAX_GS_SIZE,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE,
      eps = 0,
      by = "fgsea"
    ),
    error = function(e) e
  )

  if (inherits(obj, "error")) {

    write_empty_or_table(
      data.table(),
      output_file,
      c(
        "ID",
        "Description",
        "setSize",
        "enrichmentScore",
        "NES",
        "pvalue",
        "p.adjust",
        "qvalue",
        "rank",
        "leading_edge",
        "core_enrichment"
      )
    )

    save_empty_plot(
      paste0(
        "ROSMAP — GSEA GO:BP | ",
        model_name
      ),
      paste0("Erro: ", conditionMessage(obj)),
      plot_file
    )

    return(
      list(
        result = data.table(),
        status = data.table(
          analysis_type = "GSEA",
          database = "GO_BP",
          model = model_name,
          direction = "ranked_all_genes",
          n_input_genes = length(gene_list),
          status = "error",
          message = conditionMessage(obj)
        )
      )
    )
  }

  dt <- cp_to_dt(obj)

  write_empty_or_table(
    dt,
    output_file
  )

  plot_gsea(
    dt,
    paste0(
      "ROSMAP — GSEA GO:BP | ",
      model_name
    ),
    plot_file
  )

  list(
    result = dt,
    status = data.table(
      analysis_type = "GSEA",
      database = "GO_BP",
      model = model_name,
      direction = "ranked_all_genes",
      n_input_genes = length(gene_list),
      status = "completed",
      message = NA_character_
    )
  )
}


run_gsea_kegg <- function(
  gene_list,
  model_name
) {

  output_file <- file.path(
    gsea_dir,
    paste0(
      "GSEA_KEGG_",
      model_name,
      ".tsv"
    )
  )

  plot_file <- file.path(
    plot_dir,
    paste0(
      "GSEA_KEGG_",
      model_name,
      ".png"
    )
  )

  obj <- tryCatch(
    gseKEGG(
      geneList = gene_list,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = GSEA_MIN_GS_SIZE,
      maxGSSize = GSEA_MAX_GS_SIZE,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE,
      eps = 0,
      by = "fgsea"
    ),
    error = function(e) e
  )

  if (inherits(obj, "error")) {

    write_empty_or_table(
      data.table(),
      output_file,
      c(
        "ID",
        "Description",
        "setSize",
        "enrichmentScore",
        "NES",
        "pvalue",
        "p.adjust",
        "qvalue",
        "rank",
        "leading_edge",
        "core_enrichment"
      )
    )

    save_empty_plot(
      paste0(
        "ROSMAP — GSEA KEGG | ",
        model_name
      ),
      paste0("Erro: ", conditionMessage(obj)),
      plot_file
    )

    return(
      list(
        result = data.table(),
        status = data.table(
          analysis_type = "GSEA",
          database = "KEGG",
          model = model_name,
          direction = "ranked_all_genes",
          n_input_genes = length(gene_list),
          status = "error",
          message = conditionMessage(obj)
        )
      )
    )
  }

  dt <- cp_to_dt(obj)

  write_empty_or_table(
    dt,
    output_file
  )

  plot_gsea(
    dt,
    paste0(
      "ROSMAP — GSEA KEGG | ",
      model_name
    ),
    plot_file
  )

  list(
    result = dt,
    status = data.table(
      analysis_type = "GSEA",
      database = "KEGG",
      model = model_name,
      direction = "ranked_all_genes",
      n_input_genes = length(gene_list),
      status = "completed",
      message = NA_character_
    )
  )
}


# ==============================================================================
# 10. EXECUTAR GSEA
# ==============================================================================

cat("\nExecutando GSEA...\n")

gsea_results <- list()
gsea_status <- list()

rankings <- list(
  main = main_ranked,
  adjusted = adjusted_ranked
)

for (model_name in names(rankings)) {

  gene_list <- rankings[[model_name]]

  cat(
    "GSEA ",
    model_name,
    ": ",
    fmt_int(length(gene_list)),
    " genes ranqueados.\n",
    sep = ""
  )

  go_run <- run_gsea_go(
    gene_list,
    model_name
  )

  kegg_run <- run_gsea_kegg(
    gene_list,
    model_name
  )

  gsea_results[[paste(
    model_name,
    "GO_BP",
    sep = "__"
  )]] <- go_run$result

  gsea_results[[paste(
    model_name,
    "KEGG",
    sep = "__"
  )]] <- kegg_run$result

  gsea_status[[length(gsea_status) + 1L]] <- go_run$status
  gsea_status[[length(gsea_status) + 1L]] <- kegg_run$status
}

gsea_status_dt <- rbindlist(
  gsea_status,
  fill = TRUE
)

fwrite(
  gsea_status_dt,
  gsea_status_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 11. RESUMO DAS ANÁLISES
# ==============================================================================

count_sig <- function(dt) {

  if (
    nrow(dt) == 0 ||
      !"p.adjust" %in% names(dt)
  ) {
    return(0L)
  }

  sum(
    !is.na(dt$p.adjust) &
      dt$p.adjust < FUNCTIONAL_ALPHA
  )
}

count_gsea_direction <- function(dt, positive = TRUE) {

  if (
    nrow(dt) == 0 ||
      !all(c("p.adjust", "NES") %in% names(dt))
  ) {
    return(0L)
  }

  if (positive) {
    return(
      sum(
        !is.na(dt$p.adjust) &
          dt$p.adjust < FUNCTIONAL_ALPHA &
          dt$NES > 0,
        na.rm = TRUE
      )
    )
  }

  sum(
    !is.na(dt$p.adjust) &
      dt$p.adjust < FUNCTIONAL_ALPHA &
      dt$NES < 0,
    na.rm = TRUE
  )
}

summary_rows <- list()

for (nm in names(ora_results)) {

  dt <- ora_results[[nm]]
  pieces <- strsplit(nm, "__", fixed = TRUE)[[1]]

  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    analysis_type = "ORA",
    model = pieces[1],
    direction = pieces[2],
    database = pieces[3],
    significant_terms_FDR_0.05 = count_sig(dt),
    significant_positive_NES = NA_integer_,
    significant_negative_NES = NA_integer_
  )
}

for (nm in names(gsea_results)) {

  dt <- gsea_results[[nm]]
  pieces <- strsplit(nm, "__", fixed = TRUE)[[1]]

  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    analysis_type = "GSEA",
    model = pieces[1],
    direction = "ranked_all_genes",
    database = pieces[2],
    significant_terms_FDR_0.05 = count_sig(dt),
    significant_positive_NES = count_gsea_direction(dt, TRUE),
    significant_negative_NES = count_gsea_direction(dt, FALSE)
  )
}

summary_dt <- rbindlist(
  summary_rows,
  fill = TRUE
)

fwrite(
  summary_dt,
  summary_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 12. COMPARAÇÃO FUNCIONAL PRINCIPAL × AJUSTADO
# ==============================================================================

get_sig_ids <- function(dt) {

  if (
    nrow(dt) == 0 ||
      !all(c("ID", "p.adjust") %in% names(dt))
  ) {
    return(character())
  }

  unique(
    as.character(
      dt[
        !is.na(p.adjust) &
          p.adjust < FUNCTIONAL_ALPHA,
        ID
      ]
    )
  )
}


# ------------------------------------------------------------------------------
# Comparação termo a termo da GSEA.
# Inclui termos significativos em pelo menos um dos modelos.
# ------------------------------------------------------------------------------

compare_gsea_database <- function(database_name) {

  main_dt <- gsea_results[[paste(
    "main",
    database_name,
    sep = "__"
  )]]

  adjusted_dt <- gsea_results[[paste(
    "adjusted",
    database_name,
    sep = "__"
  )]]

  if (nrow(main_dt) == 0 && nrow(adjusted_dt) == 0) {
    return(data.table())
  }

  main_sub <- if (nrow(main_dt) > 0) {
    main_dt[
      ,
      .(
        ID = as.character(ID),
        Description_main = Description,
        NES_main = NES,
        padj_main = p.adjust
      )
    ]
  } else {
    data.table(
      ID = character(),
      Description_main = character(),
      NES_main = numeric(),
      padj_main = numeric()
    )
  }

  adjusted_sub <- if (nrow(adjusted_dt) > 0) {
    adjusted_dt[
      ,
      .(
        ID = as.character(ID),
        Description_adjusted = Description,
        NES_adjusted = NES,
        padj_adjusted = p.adjust
      )
    ]
  } else {
    data.table(
      ID = character(),
      Description_adjusted = character(),
      NES_adjusted = numeric(),
      padj_adjusted = numeric()
    )
  }

  cmp <- merge(
    main_sub,
    adjusted_sub,
    by = "ID",
    all = TRUE,
    sort = FALSE
  )

  cmp[
    ,
    Description := fifelse(
      !is.na(Description_main),
      Description_main,
      Description_adjusted
    )
  ]

  cmp[
    ,
    significant_main :=
      !is.na(padj_main) &
      padj_main < FUNCTIONAL_ALPHA
  ]

  cmp[
    ,
    significant_adjusted :=
      !is.na(padj_adjusted) &
      padj_adjusted < FUNCTIONAL_ALPHA
  ]

  cmp[
    ,
    significant_either :=
      significant_main |
      significant_adjusted
  ]

  cmp <- cmp[
    significant_either == TRUE
  ]

  cmp[
    ,
    shared_significant :=
      significant_main &
      significant_adjusted
  ]

  cmp[
    ,
    NES_same_direction :=
      shared_significant &
      sign(NES_main) == sign(NES_adjusted)
  ]

  cmp[
    ,
    database := database_name
  ]

  cmp
}

gsea_cmp_go <- compare_gsea_database("GO_BP")
gsea_cmp_kegg <- compare_gsea_database("KEGG")

gsea_comparison_terms <- rbindlist(
  list(
    gsea_cmp_go,
    gsea_cmp_kegg
  ),
  fill = TRUE
)

write_empty_or_table(
  gsea_comparison_terms,
  comparison_terms_output,
  c(
    "ID",
    "Description",
    "database",
    "NES_main",
    "padj_main",
    "NES_adjusted",
    "padj_adjusted",
    "significant_main",
    "significant_adjusted",
    "shared_significant",
    "NES_same_direction"
  )
)


# ------------------------------------------------------------------------------
# Resumo de estabilidade:
#   ORA: conjunto ALL de DEGs
#   GSEA: todos os genes ranqueados
# ------------------------------------------------------------------------------

comparison_rows <- list()

for (db in c("GO_BP", "KEGG")) {

  main_ora <- ora_results[[paste(
    "main",
    "all",
    db,
    sep = "__"
  )]]

  adjusted_ora <- ora_results[[paste(
    "adjusted",
    "all",
    db,
    sep = "__"
  )]]

  main_ids <- get_sig_ids(main_ora)
  adjusted_ids <- get_sig_ids(adjusted_ora)

  comparison_rows[[length(comparison_rows) + 1L]] <- data.table(
    analysis_type = "ORA_all_DEGs",
    database = db,
    significant_main = length(main_ids),
    significant_adjusted = length(adjusted_ids),
    shared_significant = length(intersect(main_ids, adjusted_ids)),
    union_significant = length(union(main_ids, adjusted_ids)),
    Jaccard = safe_jaccard(main_ids, adjusted_ids),
    shared_same_direction = NA_integer_,
    shared_opposite_direction = NA_integer_
  )

  main_gsea <- gsea_results[[paste(
    "main",
    db,
    sep = "__"
  )]]

  adjusted_gsea <- gsea_results[[paste(
    "adjusted",
    db,
    sep = "__"
  )]]

  main_ids <- get_sig_ids(main_gsea)
  adjusted_ids <- get_sig_ids(adjusted_gsea)
  shared_ids <- intersect(main_ids, adjusted_ids)

  cmp_dt <- if (db == "GO_BP") {
    gsea_cmp_go
  } else {
    gsea_cmp_kegg
  }

  n_same <- if (nrow(cmp_dt) == 0) {
    0L
  } else {
    sum(cmp_dt$NES_same_direction, na.rm = TRUE)
  }

  n_shared <- length(shared_ids)

  comparison_rows[[length(comparison_rows) + 1L]] <- data.table(
    analysis_type = "GSEA_ranked_all_genes",
    database = db,
    significant_main = length(main_ids),
    significant_adjusted = length(adjusted_ids),
    shared_significant = n_shared,
    union_significant = length(union(main_ids, adjusted_ids)),
    Jaccard = safe_jaccard(main_ids, adjusted_ids),
    shared_same_direction = n_same,
    shared_opposite_direction = max(0L, n_shared - n_same)
  )
}

comparison_summary <- rbindlist(
  comparison_rows,
  fill = TRUE
)

fwrite(
  comparison_summary,
  comparison_summary_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 13. FIGURA DE ESTABILIDADE DA GSEA
# ==============================================================================

if (nrow(gsea_comparison_terms) == 0) {

  save_empty_plot(
    "ROSMAP — Estabilidade funcional da GSEA",
    "Nenhum termo significativo em GO:BP ou KEGG",
    gsea_stability_plot_output
  )

} else {

  plot_dt <- gsea_comparison_terms[
    is.finite(NES_main) &
      is.finite(NES_adjusted)
  ]

  if (nrow(plot_dt) == 0) {

    save_empty_plot(
      "ROSMAP — Estabilidade funcional da GSEA",
      "Não há termos com NES disponível nos dois modelos",
      gsea_stability_plot_output
    )

  } else {

    plot_dt[
      ,
      significance_class := fifelse(
        significant_main &
          significant_adjusted,
        "Significativo em ambos",
        fifelse(
          significant_main,
          "Somente principal",
          "Somente ajustado"
        )
      )
    ]

    p <- ggplot(
      plot_dt,
      aes(
        x = NES_main,
        y = NES_adjusted,
        color = significance_class,
        shape = database
      )
    ) +
      geom_hline(
        yintercept = 0,
        linewidth = 0.3
      ) +
      geom_vline(
        xintercept = 0,
        linewidth = 0.3
      ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = 2
      ) +
      geom_point(
        alpha = 0.75,
        size = 2.4
      ) +
      labs(
        title = "ROSMAP — Estabilidade funcional da GSEA",
        subtitle = "Termos significativos em pelo menos um modelo",
        x = "NES — modelo principal",
        y = "NES — modelo ajustado",
        color = "Significância",
        shape = "Base"
      ) +
      theme_bw(
        base_size = 12
      )

    ggsave(
      gsea_stability_plot_output,
      p,
      width = 8,
      height = 7,
      dpi = 300
    )
  }
}


# ==============================================================================
# 14. SESSION INFO
# ==============================================================================

sink(session_info_output)
print(sessionInfo())
sink()


# ==============================================================================
# 15. RELATÓRIO TEXTUAL
# ==============================================================================

get_summary_value <- function(
  analysis_type_value,
  model_value,
  direction_value,
  database_value,
  field
) {

  row <- summary_dt[
    analysis_type == analysis_type_value &
      model == model_value &
      direction == direction_value &
      database == database_value
  ]

  if (nrow(row) == 0) {
    return(NA)
  }

  row[[field]][1]
}


sink(report_output)

cat(
  "============================================================\n"
)
cat(
  "FASE 01c — 06_pathway_enrichment_ROSMAP.R\n"
)
cat(
  "ENRIQUECIMENTO FUNCIONAL — MODELO PRINCIPAL E AJUSTADO\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat("1. ENTRADA\n")
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Genes testados / universo funcional: ",
  fmt_int(n_universe),
  "\n",
  sep = ""
)

cat(
  "DEGs principais — modelo principal: ",
  fmt_int(length(main_sets$all)),
  "\n",
  sep = ""
)

cat(
  "  Up em AD: ",
  fmt_int(length(main_sets$up)),
  "\n",
  sep = ""
)

cat(
  "  Down em AD: ",
  fmt_int(length(main_sets$down)),
  "\n\n",
  sep = ""
)

cat(
  "DEGs principais — modelo ajustado: ",
  fmt_int(length(adjusted_sets$all)),
  "\n",
  sep = ""
)

cat(
  "  Up em AD: ",
  fmt_int(length(adjusted_sets$up)),
  "\n",
  sep = ""
)

cat(
  "  Down em AD: ",
  fmt_int(length(adjusted_sets$down)),
  "\n\n",
  sep = ""
)


cat(
  "------------------------------------------------------------\n"
)
cat("2. METODOLOGIA\n")
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Bases:\n"
)
cat(
  "  - GO Biological Process (GO:BP)\n"
)
cat(
  "  - KEGG\n\n"
)

cat(
  "Critério funcional:\n"
)
cat(
  "  BH-FDR / p.adjust < ",
  FUNCTIONAL_ALPHA,
  "\n\n",
  sep = ""
)

cat(
  "ORA:\n"
)
cat(
  "  alvo = DEGs do critério principal\n"
)
cat(
  "  universo = todos os ",
  fmt_int(n_universe),
  " genes testados\n",
  sep = ""
)
cat(
  "  mínimo para executar ORA = ",
  MIN_ORA_GENES,
  " genes\n\n",
  sep = ""
)

cat(
  "GSEA:\n"
)
cat(
  "  ranking = Wald statistic do DESeq2\n"
)
cat(
  "  usa todos os genes com statistic disponível\n"
)
cat(
  "  minGSSize = ",
  GSEA_MIN_GS_SIZE,
  "\n",
  sep = ""
)
cat(
  "  maxGSSize = ",
  GSEA_MAX_GS_SIZE,
  "\n\n",
  sep = ""
)

cat(
  "Interpretação GSEA:\n"
)
cat(
  "  NES > 0 → enriquecimento na direção AD\n"
)
cat(
  "  NES < 0 → enriquecimento na direção Control / reduzido em AD\n\n"
)


cat(
  "------------------------------------------------------------\n"
)
cat("3. ORA — TERMOS SIGNIFICATIVOS\n")
cat(
  "------------------------------------------------------------\n"
)

for (model_name in c("main", "adjusted")) {

  cat(
    "Modelo ",
    model_name,
    ":\n",
    sep = ""
  )

  for (direction_name in c("all", "up", "down")) {

    go_n <- get_summary_value(
      "ORA",
      model_name,
      direction_name,
      "GO_BP",
      "significant_terms_FDR_0.05"
    )

    kegg_n <- get_summary_value(
      "ORA",
      model_name,
      direction_name,
      "KEGG",
      "significant_terms_FDR_0.05"
    )

    cat(
      "  ",
      direction_name,
      ": GO:BP = ",
      ifelse(is.na(go_n), "NA", fmt_int(go_n)),
      "; KEGG = ",
      ifelse(is.na(kegg_n), "NA", fmt_int(kegg_n)),
      "\n",
      sep = ""
    )
  }

  cat("\n")
}

cat(
  "OBSERVAÇÃO: listas ORA com < ",
  MIN_ORA_GENES,
  " genes são puladas e registradas no arquivo de status.\n\n",
  sep = ""
)


cat(
  "------------------------------------------------------------\n"
)
cat("4. GSEA — TERMOS SIGNIFICATIVOS\n")
cat(
  "------------------------------------------------------------\n"
)

for (model_name in c("main", "adjusted")) {

  cat(
    "Modelo ",
    model_name,
    ":\n",
    sep = ""
  )

  for (db in c("GO_BP", "KEGG")) {

    total_n <- get_summary_value(
      "GSEA",
      model_name,
      "ranked_all_genes",
      db,
      "significant_terms_FDR_0.05"
    )

    pos_n <- get_summary_value(
      "GSEA",
      model_name,
      "ranked_all_genes",
      db,
      "significant_positive_NES"
    )

    neg_n <- get_summary_value(
      "GSEA",
      model_name,
      "ranked_all_genes",
      db,
      "significant_negative_NES"
    )

    cat(
      "  ",
      db,
      ": total = ",
      ifelse(is.na(total_n), "NA", fmt_int(total_n)),
      "; NES+ = ",
      ifelse(is.na(pos_n), "NA", fmt_int(pos_n)),
      "; NES- = ",
      ifelse(is.na(neg_n), "NA", fmt_int(neg_n)),
      "\n",
      sep = ""
    )
  }

  cat("\n")
}


cat(
  "------------------------------------------------------------\n"
)
cat("5. ESTABILIDADE FUNCIONAL PRINCIPAL × AJUSTADO\n")
cat(
  "------------------------------------------------------------\n"
)

for (i in seq_len(nrow(comparison_summary))) {

  row <- comparison_summary[i]

  cat(
    row$analysis_type,
    " | ",
    row$database,
    "\n",
    sep = ""
  )

  cat(
    "  significativo principal: ",
    fmt_int(row$significant_main),
    "\n",
    sep = ""
  )

  cat(
    "  significativo ajustado: ",
    fmt_int(row$significant_adjusted),
    "\n",
    sep = ""
  )

  cat(
    "  compartilhados: ",
    fmt_int(row$shared_significant),
    "\n",
    sep = ""
  )

  cat(
    "  união: ",
    fmt_int(row$union_significant),
    "\n",
    sep = ""
  )

  cat(
    "  Jaccard: ",
    ifelse(
      is.na(row$Jaccard),
      "NA",
      fmt_num(row$Jaccard, 4)
    ),
    "\n",
    sep = ""
  )

  if (row$analysis_type == "GSEA_ranked_all_genes") {

    cat(
      "  compartilhados com mesma direção do NES: ",
      fmt_int(row$shared_same_direction),
      "\n",
      sep = ""
    )

    cat(
      "  compartilhados com direção oposta do NES: ",
      fmt_int(row$shared_opposite_direction),
      "\n",
      sep = ""
    )
  }

  cat("\n")
}


cat(
  "------------------------------------------------------------\n"
)
cat("6. ARQUIVOS PRINCIPAIS\n")
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Resumo geral:\n  ",
  summary_output,
  "\n\n",
  sep = ""
)

cat(
  "Status ORA:\n  ",
  ora_status_output,
  "\n\n",
  sep = ""
)

cat(
  "Status GSEA:\n  ",
  gsea_status_output,
  "\n\n",
  sep = ""
)

cat(
  "Estabilidade funcional — resumo:\n  ",
  comparison_summary_output,
  "\n\n",
  sep = ""
)

cat(
  "Estabilidade funcional — termos:\n  ",
  comparison_terms_output,
  "\n\n",
  sep = ""
)

cat(
  "Figura NES principal × ajustado:\n  ",
  gsea_stability_plot_output,
  "\n\n",
  sep = ""
)


cat(
  "------------------------------------------------------------\n"
)
cat("7. INTERPRETAÇÃO PLANEJADA\n")
cat(
  "------------------------------------------------------------\n"
)

cat(
  "- A ORA avalia sobre-representação entre os genes que cruzaram o cutoff.\n"
)

cat(
  "- A GSEA avalia deslocamentos coordenados usando todos os genes ranqueados.\n"
)

cat(
  "- Como as listas de DEGs do ROSMAP são pequenas e pouco sobrepostas entre\n"
)
cat(
  "  os modelos, a GSEA é especialmente informativa para avaliar estabilidade\n"
)
cat(
  "  funcional sem depender exclusivamente do cutoff binário de DEG.\n"
)

cat(
  "- A estabilidade funcional é quantificada por sobreposição de termos,\n"
)
cat(
  "  Jaccard e concordância da direção do NES.\n"
)

cat(
  "- O enriquecimento não deve ser usado isoladamente para escolher a coorte\n"
)
cat(
  "  principal de AD.\n"
)

cat(
  "- A próxima etapa deverá comparar ROSMAP e GSE125583 em nível de genes,\n"
)
cat(
  "  efeitos e vias, conforme os critérios definidos antes da inspeção dos\n"
)
cat(
  "  resultados do ROSMAP.\n\n"
)

cat(
  "============================================================\n"
)
cat(
  "06_pathway_enrichment_ROSMAP.R CONCLUÍDO\n"
)
cat(
  "============================================================\n"
)

sink()


# ==============================================================================
# 16. RESUMO NO TERMINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "ENRIQUECIMENTO ROSMAP CONCLUÍDO\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Universo funcional: ",
  fmt_int(n_universe),
  " genes\n",
  sep = ""
)

cat(
  "DEGs principal: ",
  fmt_int(length(main_sets$all)),
  "\n",
  sep = ""
)

cat(
  "DEGs ajustado: ",
  fmt_int(length(adjusted_sets$all)),
  "\n\n",
  sep = ""
)

cat(
  "Resumo:\n  ",
  summary_output,
  "\n\n",
  sep = ""
)

cat(
  "Estabilidade funcional:\n  ",
  comparison_summary_output,
  "\n\n",
  sep = ""
)

cat(
  "Relatório:\n  ",
  report_output,
  "\n\n",
  sep = ""
)

cat(
  "Próxima etapa após revisar os resultados:\n"
)
cat(
  "  comparação ROSMAP × GSE125583\n"
)
