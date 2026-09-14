library(dplyr)
library(readr)
library(ggplot2)
library(ggVennDiagram)

results_dir <- "results/01_core_original"
datasets <- c("GSE53697", "GSE64810", "GSE68719")

deg_lists <- list()

for (ds in datasets) {
  f <- file.path(results_dir, ds, paste0(ds, "_deg_results.csv"))
  if (file.exists(f)) {
    df <- read_csv(f, show_col_types = FALSE)
    sig_genes <- df %>% filter(adj.P.Val < 0.05 & abs(logFC) > 0.58) %>% pull(gene_symbol)
    deg_lists[[ds]] <- sig_genes
  }
}

dir.create("results/01_core_original/intersection", recursive = TRUE, showWarnings = FALSE)

p_venn <- ggVennDiagram(deg_lists, category.names = c("GSE53697 (AD)", "GSE64810 (ALS)", "GSE68719 (HD)")) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = "Interseção de DEGs entre Alzheimer, ALS e Huntington")

ggsave("results/01_core_original/intersection/venn_degs.png", p_venn, width = 7, height = 6, dpi = 300)
cat("[SUCCESS] Gráfico de interseção salvo em results/01_core_original/intersection/venn_degs.png\n")