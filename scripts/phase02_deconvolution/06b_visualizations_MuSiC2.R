#!/usr/bin/env Rscript
# =============================================================================
# 06b_visualizar_diagnosticos_MuSiC2.R
# =============================================================================
# Gera gráficos de diagnóstico para a etapa 06 SEM executar o MuSiC2 novamente.
# Rodar na raiz do projeto:
#   Rscript scripts/phase02_deconvolution/06b_visualizar_diagnosticos_MuSiC2.R
#
# Entradas (geradas pelo 06):
#   results/deconvolution/MuSiC2/06_proportions_long_all.tsv
#   results/deconvolution/MuSiC2/06_run_summary.tsv
#
# Saída:
#   results/deconvolution/MuSiC2/plots/*.png
#   results/deconvolution/MuSiC2/plots/06b_zero_fraction_by_dataset_celltype.tsv
#   results/deconvolution/MuSiC2/plots/06b_convergence_summary.tsv
#
# ATENÇÃO: figuras são DIAGNÓSTICAS. Estimativas de corridas sem convergência
# completa não devem ser interpretadas como resultados biológicos validados.
# =============================================================================

for (pkg in c("data.table", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Pacote ausente: ", pkg, ". Instale-o no ambiente do projeto.")
  }
}
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

result_root <- file.path("results", "deconvolution", "MuSiC2")
plot_root <- file.path(result_root, "plots")
dir.create(plot_root, showWarnings = FALSE, recursive = TRUE)

p_file <- file.path(result_root, "06_proportions_long_all.tsv")
s_file <- file.path(result_root, "06_run_summary.tsv")
if (!file.exists(p_file) || !file.exists(s_file)) {
  stop("Arquivos consolidados do MuSiC2 ausentes. Execute da raiz do projeto.")
}

x <- fread(p_file, na.strings = c("", "NA"))
s <- fread(s_file, na.strings = c("", "NA"))
req_x <- c("dataset", "disease", "reference_version", "marker_set", "sample_id", "group", "cell_type", "proportion")
req_s <- c("dataset", "reference_version", "marker_set", "cases", "n_id_not_converged", "converged")
if (length(setdiff(req_x, names(x)))) stop("Colunas ausentes em 06_proportions_long_all.tsv")
if (length(setdiff(req_s, names(s)))) stop("Colunas ausentes em 06_run_summary.tsv")
if (anyNA(x$proportion) || any(!is.finite(x$proportion))) stop("Proporções NA/não finitas")

ct_order <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
dataset_order <- c("ROSMAP", "GSE125583", "GSE64810", "GSE68719")
x[, cell_type := factor(cell_type, levels = ct_order)]
x[, dataset := factor(dataset, levels = dataset_order)]

# Usar somente a referência principal no diagnóstico das quatro coortes.
main <- x[reference_version == "main24" & marker_set == "top300"]
if (nrow(main) == 0L) stop("Nenhuma proporção main24/top300 encontrada")
main[, group := factor(group, levels = c("Control", "AD", "HD", "PD"))]

# Cores uniformes: coortes/condições, e uma cor para cada tipo celular.
# Controles sempre cinza; cada doença mantém a mesma cor entre os gráficos.
disease_colors <- c(
  "AD" = "#7B2CBF", "HD" = "#277DA1", "PD" = "#E75480"
)
condition_colors <- c("Control" = "#94A3B8", disease_colors)
cell_colors <- c(
  "Ex" = "#E76F51", "In" = "#D69E00",
  "Ast" = "#7D9D22", "Oli" = "#009E73",
  "Opc" = "#56B4E9", "Mic" = "#0072B2",
  "End" = "#9970AB", "Per" = "#CC79A7"
)
cohort_colors <- c(
  "ROSMAP" = disease_colors[["AD"]],
  "GSE125583" = "#9F5ACD",
  "GSE64810" = disease_colors[["HD"]],
  "GSE68719" = disease_colors[["PD"]]
)
plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", colour = "#141937"),
    plot.subtitle = element_text(colour = "#475569"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.caption = element_text(colour = "#475569", size = 9)
  )

# ---------------------------------------------------------------------------
# 1. Taxa de não convergência POR AMOSTRA DE DOENÇA.
# ---------------------------------------------------------------------------
cv <- copy(s)
cv[, run_label := paste0(as.character(dataset), " / ", reference_version)]
cv[, nonconvergence_rate := n_id_not_converged / cases]
cv[, percentage_label := sprintf("%.1f%% (%d/%d)", 100 * nonconvergence_rate, n_id_not_converged, cases)]
fwrite(cv, file.path(plot_root, "06b_convergence_summary.tsv"), sep = "\t")
cv[, cohort_color := cohort_colors[as.character(dataset)]]
plot_conv <- ggplot(cv, aes(x = reorder(run_label, nonconvergence_rate),
                            y = 100 * nonconvergence_rate, fill = dataset)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = cohort_colors, name = "Coorte", drop = FALSE) +
  geom_text(aes(label = percentage_label), hjust = -0.1, size = 3.6) +
  coord_flip(ylim = c(0, max(55, 100 * max(cv$nonconvergence_rate) * 1.2)), clip = "off") +
  labs(title = "MuSiC2 — casos que não convergiram", x = NULL,
       y = "Percentual das amostras de doença (%)",
       caption = "Convergência das iterações de MuSiC2; não representa falha de todas as amostras.") +
  plot_theme
ggsave(file.path(plot_root, "06b_01_nao_convergencia.png"), plot_conv, width = 11, height = 5.6, dpi = 180)

# ---------------------------------------------------------------------------
# 2. Fração de estimativas exatamente zero por tipo e dataset.
# ---------------------------------------------------------------------------
zero_summary <- main[, .(
  n_samples = .N,
  exact_zero_samples = sum(proportion == 0),
  exact_zero_fraction = mean(proportion == 0),
  near_zero_below_0_1pct = mean(proportion < 0.001),
  mean_proportion = mean(proportion),
  median_proportion = median(proportion)
), by = .(dataset, cell_type)]
fwrite(zero_summary, file.path(plot_root, "06b_zero_fraction_by_dataset_celltype.tsv"), sep = "\t")
zero_summary[, text_color := ifelse(exact_zero_fraction >= 0.65, "white", "#141937")]
plot_zero <- ggplot(zero_summary, aes(x = dataset, y = cell_type,
                                     fill = 100 * exact_zero_fraction)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * exact_zero_fraction),
                colour = text_color), size = 3.6, show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_gradientn(colours = c("#E0F2F1", "#F9C74F", "#F3722C", "#A41638"),
                       limits = c(0, 100), name = "Zeros (%)") +
  labs(title = "MuSiC2 — estimativas exatamente iguais a zero", x = "Coorte", y = "Tipo celular",
       fill = "Zeros (%)", caption = "Referência main24/top300. Zero estimado não prova ausência da célula.") +
  plot_theme + theme(panel.grid = element_blank())
ggsave(file.path(plot_root, "06b_02_mapa_zeros.png"), plot_zero, width = 10.5, height = 6.4, dpi = 180)

# ---------------------------------------------------------------------------
# 3. Perfis composicionais: média por coorte e condição.
# ---------------------------------------------------------------------------
mean_comp <- main[, .(mean_proportion = mean(proportion)),
                  by = .(dataset, disease, group, cell_type)]
mean_comp[, condition_label := paste(as.character(dataset), as.character(group), sep = " · ")]
mean_comp[, condition_label := factor(condition_label,
                                    levels = unique(mean_comp$condition_label))]
plot_mix <- ggplot(mean_comp, aes(x = condition_label, y = 100 * mean_proportion,
                                 fill = cell_type)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = cell_colors, drop = FALSE) +
  coord_flip() +
  labs(title = "MuSiC2 — composição média estimada por coorte e grupo",
       x = NULL, y = "Contribuição estimada (%)", fill = "Tipo celular",
       caption = "Descrição, NÃO inferência. In e Opc são praticamente zerados. Não comparar coortes como amostras equivalentes.") +
  plot_theme
ggsave(file.path(plot_root, "06b_03_composicao_media.png"), plot_mix, width = 12, height = 7, dpi = 180)

# ---------------------------------------------------------------------------
# 4. Distribuições individuais por coorte. Escalas livres por tipo para detectar
# zeros mesmo nas categorias com valores extremamente baixos.
# ---------------------------------------------------------------------------
for (ds in dataset_order) {
  dd <- main[dataset == ds]
  if (nrow(dd) == 0L) next
  plot_box <- ggplot(dd, aes(x = group, y = 100 * proportion, fill = group)) +
    geom_boxplot(outlier.alpha = 0.40, outlier.size = 1.25,
                 colour = "#334155", linewidth = 0.35) +
    scale_fill_manual(values = condition_colors, drop = FALSE) +
    facet_wrap(~cell_type, scales = "free_y", ncol = 4) +
    labs(title = paste0("MuSiC2 — distribuições individuais: ", ds),
         subtitle = "Referência main24 / top300", x = "Grupo", y = "Estimativa (%)",
         caption = "Cada ponto extremo é uma amostra; eixo Y muda entre facetas. Resultados ainda não validados.") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"))
  ggsave(file.path(plot_root, paste0("06b_04_boxplots_", ds, ".png")),
         plot_box, width = 13, height = 7, dpi = 180)
}

# ---------------------------------------------------------------------------
# 5. Sensibilidade ROSMAP — exclusão de 3 doadores sobrepostos.
# ---------------------------------------------------------------------------
ra <- x[dataset == "ROSMAP" & marker_set == "top300" &
          reference_version == "main24", .(sample_id, cell_type, p24 = proportion)]
rb <- x[dataset == "ROSMAP" & marker_set == "top300" &
          reference_version == "no_ROSMAP_overlap21", .(sample_id, cell_type, p21 = proportion)]
if (nrow(ra) > 0L && nrow(rb) > 0L) {
  comp <- merge(ra, rb, by = c("sample_id", "cell_type"))
  if (nrow(comp) == 0L) stop("ROSMAP 24 vs 21 sem amostras alinhadas")
  comp[, difference_pp := 100 * (p21 - p24)]
  difference_long <- comp[, .(
    mean_absolute_difference_pp = mean(abs(difference_pp)),
    max_absolute_difference_pp = max(abs(difference_pp)),
    median_difference_pp = median(difference_pp)
  ), by = cell_type]
  fwrite(difference_long, file.path(plot_root, "06b_ROSMAP_overlap_absolute_differences.tsv"), sep = "\t")
  plot_overlap <- ggplot(difference_long,
                         aes(x = cell_type, y = mean_absolute_difference_pp,
                             fill = cell_type)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = cell_colors, guide = "none") +
    labs(title = "ROSMAP — impacto de remover 3 doadores da referência",
         x = "Tipo celular", y = "Diferença absoluta média (pontos percentuais)",
         caption = "Comparação pareada das mesmas 144 amostras; não demonstra causa da diferença.") +
    plot_theme
  ggsave(file.path(plot_root, "06b_05_ROSMAP_sensibilidade.png"), plot_overlap,
         width = 10, height = 5.5, dpi = 180)
}

cat("\nDiagnósticos 06b concluídos. Arquivos em:\n  ", plot_root, "\n", sep = "")
cat("Não foi necessário executar music2_prop() novamente.\n")
