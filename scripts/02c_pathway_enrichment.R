#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript 02c_pathway_enrichment.R <deg_results.csv> <out_prefix>")
}

deg_file   <- args[1]
out_prefix <- args[2]

deg_res <- read_csv(deg_file, show_col_types = FALSE)

# Filtrar DEGs significativos (FDR < 0.05)
sig_genes <- deg_res %>% 
  filter(adj.P.Val < 0.05 & abs(logFC) > 0.58) %>% 
  pull(gene_symbol)

cat("[INFO] Total de genes DEGs significativos para enriquecimento:", length(sig_genes), "\n")

if (length(sig_genes) > 5) {
  # Converter Symbol para Entrez ID para clusterProfiler
  gene_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Enriquecimento GO (Biological Process)
  ego <- enrichGO(gene          = gene_ids$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05)
  
  go_res <- as.data.frame(ego)
  write_csv(go_res, paste0(out_prefix, "_pathways_GO_BP.csv"))
  cat("[SUCCESS] Vias GO salvos em:", paste0(out_prefix, "_pathways_GO_BP.csv"), "\n")
} else {
  cat("[WARN] Poucos genes significativos encontrados com o corte aplicado.\n")
}
