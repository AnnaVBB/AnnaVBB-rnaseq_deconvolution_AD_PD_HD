#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tibble)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Uso: Rscript 02b_exploratory_deg.R <counts.tsv.gz> <metadata.csv> <out_prefix>")
}

tsv_file   <- args[1]
meta_file  <- args[2]
out_prefix <- args[3]

cat("[INFO] Processando dataset:", tsv_file, "\n")

# Criar pasta de saída caso não exista
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

# 1. Carregar dados brutos (.tsv.gz)
counts   <- read.delim(tsv_file, row.names = 1, check.names = FALSE)
metadata <- read_csv(meta_file, show_col_types = FALSE)

# Alinhar amostras
common_samples <- intersect(colnames(counts), metadata$sample_id)
counts   <- counts[, common_samples]
metadata <- metadata %>% 
  filter(sample_id %in% common_samples) %>% 
  slice(match(common_samples, sample_id))

# Definir fator do grupo (Control vs Disease)
metadata$group <- factor(metadata$group, levels = c("Control", unique(metadata$group[metadata$group != "Control"])[1]))
disease_label  <- levels(metadata$group)[2]

cat("[INFO] Amostras:", ncol(counts), "| Grupo Controle vs", disease_label, "\n")

# 2. Análise com limma-voom
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

design <- model.matrix(~ group, data = metadata)
colnames(design) <- c("Intercept", "Disease")

v   <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

res_table <- topTable(fit, coef = "Disease", number = Inf) %>%
  tibble::rownames_to_column(var = "gene_symbol") %>%
  as_tibble()

# Salvar Tabela Completa de DEGs
deg_out_csv <- paste0(out_prefix, "_deg_results.csv")
write_csv(res_table, deg_out_csv)
cat("[SUCCESS] Resultados DEG salvos em:", deg_out_csv, "\n")

# 3. Gerar Gráfico Volcano Plot
pdf(paste0(out_prefix, "_volcano.pdf"), width = 6, height = 5)
p <- ggplot(res_table, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = adj.P.Val < 0.05 & abs(logFC) > 0.58), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("grey60", "red3")) +
  theme_minimal() +
  labs(title = paste("Volcano Plot:", disease_label, "vs Control"),
       subtitle = paste("Dataset:", basename(tsv_file)),
       x = "Log2 Fold Change", y = "-Log10 P-Value", color = "FDR < 0.05 & |LFC| > 0.58")
print(p)
dev.off()

cat("[SUCCESS] Volcano plot salvo em:", paste0(out_prefix, "_volcano.pdf"), "\n")
