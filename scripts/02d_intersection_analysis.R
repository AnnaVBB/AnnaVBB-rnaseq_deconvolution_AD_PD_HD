library(dplyr)
library(readr)
library(ggplot2)
library(ggVennDiagram)

results_dir <- "results/01_core_original"

datasets <- c(
  "GSE53697",
  "GSE64810",
  "GSE68719"
)

disease_names <- c(
  GSE53697 = "AD",
  GSE64810 = "HD",
  GSE68719 = "PD"
)

deg_lists <- list()

for (ds in datasets) {

  f <- file.path(
    results_dir,
    ds,
    paste0(ds, "_deg_results.csv")
  )

  if (!file.exists(f)) {
    stop(
      "Arquivo DEG não encontrado: ",
      f
    )
  }

  df <- read_csv(
    f,
    show_col_types = FALSE
  )

  sig_genes <- df %>%
    filter(
      adj.P.Val < 0.05,
      abs(logFC) > 0.58
    ) %>%
    pull(gene_id) %>%
    unique()

  deg_lists[[ds]] <- sig_genes

  cat(
    "[INFO]",
    ds,
    "(",
    disease_names[[ds]],
    "):",
    length(sig_genes),
    "DEGs\n"
  )
}

if (length(deg_lists) != 3) {
  stop("Não foi possível carregar os três datasets.")
}

dir.create(
  file.path(results_dir, "intersection"),
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Venn
# ------------------------------------------------------------

p_venn <- ggVennDiagram(
  deg_lists,
  category.names = c(
    "GSE53697 (AD)",
    "GSE64810 (HD)",
    "GSE68719 (PD)"
  )
) +
  labs(
    title = "Interseção de DEGs entre AD, HD e PD"
  )

ggsave(
  file.path(
    results_dir,
    "intersection",
    "venn_degs.png"
  ),
  p_venn,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Genes comuns
# ------------------------------------------------------------

common_genes <- Reduce(
  intersect,
  deg_lists
)

write_csv(
  tibble(
    gene_id = common_genes
  ),
  file.path(
    results_dir,
    "intersection",
    "common_genes.csv"
  )
)

cat(
  "[SUCCESS] Genes comuns:",
  length(common_genes),
  "\n"
)

cat(
  "[SUCCESS] Venn salvo em:",
  file.path(
    results_dir,
    "intersection",
    "venn_degs.png"
  ),
  "\n"
)