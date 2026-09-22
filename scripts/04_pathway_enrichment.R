#!/usr/bin/env Rscript

# ============================================================
# 04_pathway_enrichment.R
#[3]
# Enriquecimento funcional dos DEGs principais:
#   - GO Biological Process
#   - KEGG
#
# Estratégia:
#   - entrada: deg_all.csv produzido por 03_exploratory_deg.R
#   - DEG principal: adj.P.Val < 0.05 & |logFC| > 0.58
#   - análises separadas:
#       * todos os DEGs
#       * genes Up
#       * genes Down
#   - background: todos os genes testados no próprio dataset
#
# Uso:
#
# Rscript scripts/04_pathway_enrichment.R \
#   GSE64810 \
#   results/01_core_original/GSE64810/deg/deg_all.csv \
#   results/01_core_original/GSE64810/enrichment
#
# ============================================================


# ------------------------------------------------------------
# 1. Pacotes
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})


# ------------------------------------------------------------
# 2. Argumentos
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "\nUso:\n",
      "Rscript scripts/04_pathway_enrichment.R ",
      "<dataset> <deg_all.csv> <output_dir>\n"
    )
  )
}

dataset <- args[1]
deg_file <- args[2]
output_dir <- args[3]


# ------------------------------------------------------------
# 3. Configuração do dataset
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
  stop("Dataset não reconhecido: ", dataset)
}

disease <- dataset_config[[dataset]]$disease
disease_label <- dataset_config[[dataset]]$disease_label


# ------------------------------------------------------------
# 4. Diretórios
# ------------------------------------------------------------

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

plot_dir <- file.path(
  output_dir,
  "plots"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 5. Ler DEG
# ------------------------------------------------------------

message("Lendo: ", deg_file)

deg <- read_csv(
  deg_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 6. Validação
# ------------------------------------------------------------

required_columns <- c(
  "gene_id",
  "logFC",
  "P.Value",
  "adj.P.Val",
  "significant_primary",
  "direction_primary"
)

missing_columns <- setdiff(
  required_columns,
  colnames(deg)
)

if (length(missing_columns) > 0) {
  stop(
    "Colunas ausentes no arquivo DEG: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 7. Conferir genes principais
# ------------------------------------------------------------

deg_primary <- deg %>%
  filter(
    significant_primary == TRUE
  )

deg_up <- deg_primary %>%
  filter(
    direction_primary == "Up"
  )

deg_down <- deg_primary %>%
  filter(
    direction_primary == "Down"
  )

message("")
message("============================================")
message("ENRIQUECIMENTO — ", dataset)
message("============================================")
message("Doença: ", disease_label)
message("Genes testados: ", nrow(deg))
message("DEGs principais: ", nrow(deg_primary))
message("  Up: ", nrow(deg_up))
message("  Down: ", nrow(deg_down))
message("============================================")
message("")


# ------------------------------------------------------------
# 8. Limpar Ensembl
# ------------------------------------------------------------

clean_ensembl <- function(x) {

  x <- as.character(x)

  x <- sub(
    "\\..*$",
    "",
    x
  )

  x[
    is.na(x) |
    x == ""
  ] <- NA_character_

  x
}


deg <- deg %>%
  mutate(
    ensembl_clean = clean_ensembl(gene_id)
  )

deg_primary <- deg %>%
  filter(
    significant_primary == TRUE
  )

deg_up <- deg_primary %>%
  filter(
    direction_primary == "Up"
  )

deg_down <- deg_primary %>%
  filter(
    direction_primary == "Down"
  )


# ------------------------------------------------------------
# 9. Conversão Ensembl -> Entrez
# ------------------------------------------------------------

convert_to_entrez <- function(ensembl_ids) {

  ensembl_ids <- unique(
    na.omit(
      ensembl_ids
    )
  )

  if (length(ensembl_ids) == 0) {

    return(
      tibble(
        ENSEMBL = character(),
        ENTREZID = character()
      )
    )
  }

  converted <- suppressMessages(
    bitr(
      ensembl_ids,
      fromType = "ENSEMBL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )

  converted %>%
    distinct(
      ENSEMBL,
      ENTREZID
    )
}


# ------------------------------------------------------------
# 10. Background
# ------------------------------------------------------------

background_map <- convert_to_entrez(
  deg$ensembl_clean
)

background_entrez <- unique(
  background_map$ENTREZID
)

message(
  "Genes no background com ENTREZ: ",
  length(background_entrez)
)


# ------------------------------------------------------------
# 11. Função para preparar conjunto de genes
# ------------------------------------------------------------

prepare_gene_set <- function(data) {

  mapping <- convert_to_entrez(
    data$ensembl_clean
  )

  unique(
    mapping$ENTREZID
  )
}


genes_all <- prepare_gene_set(
  deg_primary
)

genes_up <- prepare_gene_set(
  deg_up
)

genes_down <- prepare_gene_set(
  deg_down
)

message("")
message("DEGs convertidos para ENTREZ:")
message("  All:  ", length(genes_all))
message("  Up:   ", length(genes_up))
message("  Down: ", length(genes_down))
message("")


# ------------------------------------------------------------
# 12. Função: tabela vazia
# ------------------------------------------------------------

empty_enrichment_table <- function() {

  tibble(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer()
  )
}


# ------------------------------------------------------------
# 13. Função: GO Biological Process
# ------------------------------------------------------------

run_go_bp <- function(
  genes,
  universe
) {

  if (length(genes) < 2) {

    message(
      "GO BP não executado: menos de 2 genes."
    )

    return(
      empty_enrichment_table()
    )
  }

  result <- suppressMessages(
    enrichGO(
      gene = genes,
      universe = universe,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
  )

  result_df <- as.data.frame(
    result
  )

  if (nrow(result_df) == 0) {
    return(
      empty_enrichment_table()
    )
  }

  as_tibble(
    result_df
  )
}


# ------------------------------------------------------------
# 14. Função: KEGG
# ------------------------------------------------------------

run_kegg <- function(
  genes,
  universe
) {

  if (length(genes) < 2) {

    message(
      "KEGG não executado: menos de 2 genes."
    )

    return(
      empty_enrichment_table()
    )
  }

  result <- suppressMessages(
    enrichKEGG(
      gene = genes,
      universe = universe,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  )

  result_df <- as.data.frame(
    result
  )

  if (nrow(result_df) == 0) {
    return(
      empty_enrichment_table()
    )
  }

  as_tibble(
    result_df
  )
}


# ------------------------------------------------------------
# 15. Executar enriquecimentos
# ------------------------------------------------------------

message("Executando GO BP...")

go_all <- run_go_bp(
  genes_all,
  background_entrez
)

go_up <- run_go_bp(
  genes_up,
  background_entrez
)

go_down <- run_go_bp(
  genes_down,
  background_entrez
)


message("Executando KEGG...")

kegg_all <- run_kegg(
  genes_all,
  background_entrez
)

kegg_up <- run_kegg(
  genes_up,
  background_entrez
)

kegg_down <- run_kegg(
  genes_down,
  background_entrez
)


# ------------------------------------------------------------
# 16. Salvar tabelas
# ------------------------------------------------------------

write_csv(
  go_all,
  file.path(
    output_dir,
    "go_bp_all.csv"
  )
)

write_csv(
  go_up,
  file.path(
    output_dir,
    "go_bp_up.csv"
  )
)

write_csv(
  go_down,
  file.path(
    output_dir,
    "go_bp_down.csv"
  )
)

write_csv(
  kegg_all,
  file.path(
    output_dir,
    "kegg_all.csv"
  )
)

write_csv(
  kegg_up,
  file.path(
    output_dir,
    "kegg_up.csv"
  )
)

write_csv(
  kegg_down,
  file.path(
    output_dir,
    "kegg_down.csv"
  )
)


# ------------------------------------------------------------
# 17. Função para plotar resultados
# ------------------------------------------------------------
save_dotplot <- function(
  enrichment_df,
  title,
  filename
) {

  if (nrow(enrichment_df) == 0) {

    message(
      "Plot não produzido: ",
      title,
      " — nenhum termo enriquecido."
    )

    return(
      invisible(NULL)
    )
  }

  plot_data <- enrichment_df %>%
    arrange(
      p.adjust
    ) %>%
    slice_head(
      n = 15
    ) %>%
    mutate(
      Description = factor(
        Description,
        levels = rev(
          unique(Description)
        )
      ),
      minus_log10_fdr = -log10(
        pmax(
          p.adjust,
          .Machine$double.xmin
        )
      )
    )

  p <- ggplot(
    plot_data,
    aes(
      x = minus_log10_fdr,
      y = Description,
      size = Count,
      color = p.adjust
    )
  ) +
    geom_point(
      alpha = 0.85
    ) +
    scale_color_gradient(
      low = "red",
      high = "blue"
    ) +
    labs(
      title = title,
      x = expression(-log[10]("FDR")),
      y = NULL,
      size = "Genes",
      color = "FDR"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(
        size = 8
      )
    )

  ggsave(
    filename = filename,
    plot = p,
    width = 10,
    height = 7,
    dpi = 300
  )
}


# ------------------------------------------------------------
# 18. Plots GO
# ------------------------------------------------------------

save_dotplot(
  go_all,
  paste0(
    dataset,
    " — GO BP — todos os DEGs"
  ),
  file.path(
    plot_dir,
    "go_bp_all.png"
  )
)

save_dotplot(
  go_up,
  paste0(
    dataset,
    " — GO BP — Up"
  ),
  file.path(
    plot_dir,
    "go_bp_up.png"
  )
)

save_dotplot(
  go_down,
  paste0(
    dataset,
    " — GO BP — Down"
  ),
  file.path(
    plot_dir,
    "go_bp_down.png"
  )
)


# ------------------------------------------------------------
# 19. Plots KEGG
# ------------------------------------------------------------

save_dotplot(
  kegg_all,
  paste0(
    dataset,
    " — KEGG — todos os DEGs"
  ),
  file.path(
    plot_dir,
    "kegg_all.png"
  )
)

save_dotplot(
  kegg_up,
  paste0(
    dataset,
    " — KEGG — Up"
  ),
  file.path(
    plot_dir,
    "kegg_up.png"
  )
)

save_dotplot(
  kegg_down,
  paste0(
    dataset,
    " — KEGG — Down"
  ),
  file.path(
    plot_dir,
    "kegg_down.png"
  )
)


# ------------------------------------------------------------
# 20. Resumo
# ------------------------------------------------------------

summary_table <- tibble(

  dataset = dataset,

  disease = disease,

  genes_tested = nrow(deg),

  background_ensembl = length(
    unique(
      na.omit(
        deg$ensembl_clean
      )
    )
  ),

  background_entrez = length(
    background_entrez
  ),

  deg_total = nrow(
    deg_primary
  ),

  deg_up = nrow(
    deg_up
  ),

  deg_down = nrow(
    deg_down
  ),

  deg_entrez_total = length(
    genes_all
  ),

  deg_entrez_up = length(
    genes_up
  ),

  deg_entrez_down = length(
    genes_down
  ),

  mapping_rate_total = ifelse(
  nrow(deg_primary) > 0,
  length(genes_all) / nrow(deg_primary) * 100,
  NA_real_
  ),

  mapping_rate_up = ifelse(
    nrow(deg_up) > 0,
    length(genes_up) / nrow(deg_up) * 100,
    NA_real_
  ),

  mapping_rate_down = ifelse(
    nrow(deg_down) > 0,
    length(genes_down) / nrow(deg_down) * 100,
    NA_real_
  ),

  go_all_terms = nrow(
    go_all
  ),

  go_up_terms = nrow(
    go_up
  ),

  go_down_terms = nrow(
    go_down
  ),

  kegg_all_terms = nrow(
    kegg_all
  ),

  kegg_up_terms = nrow(
    kegg_up
  ),

  kegg_down_terms = nrow(
    kegg_down
  )
)

write_csv(
  summary_table,
  file.path(
    output_dir,
    "enrichment_summary.csv"
  )
)


# ------------------------------------------------------------
# 21. Resumo no terminal
# ------------------------------------------------------------

message("")
message("============================================")
message("ENRIQUECIMENTO FINALIZADO — ", dataset)
message("============================================")

message(
  "Background ENTREZ: ",
  length(background_entrez)
)

message(
  "DEGs ENTREZ: ",
  length(genes_all)
)

message("")
message("GO BP:")
message("  All:  ", nrow(go_all))
message("  Up:   ", nrow(go_up))
message("  Down: ", nrow(go_down))

message("")
message("KEGG:")
message("  All:  ", nrow(kegg_all))
message("  Up:   ", nrow(kegg_up))
message("  Down: ", nrow(kegg_down))

message("")
message(
  "Resultados: ",
  output_dir
)

message("============================================")
