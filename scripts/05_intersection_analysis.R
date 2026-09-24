#!/usr/bin/env Rscript

# ============================================================
# 05_intersection_analysis.R
#
# Integração transdiagnóstica:
# Alzheimer (AD), Huntington (HD) e Parkinson (PD)
#
# Estratégia:
#   1. Construir universo comum testável
#   2. Identificar interseções estritas de DEGs principais
#   3. Classificar padrões de significância
#   4. Avaliar concordância direcional
#   5. Comparar logFC entre doenças
#   6. Repetir interseção com o critério comparativo do artigo
#
# Critério principal:
#   FDR < 0.05 e |logFC| > 0.58
#
# Critério comparativo:
#   p < 0.05 e |logFC| > 1
#
# IMPORTANTE:
#   ausência de significância != ausência de efeito
#   ausência do universo != gene não significativo
# ============================================================


# ------------------------------------------------------------
# 1. Pacotes
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggVennDiagram)
})


# ------------------------------------------------------------
# 2. Paleta das doenças
# ------------------------------------------------------------

disease_colors <- c(
  "AD" = "#7B2CBF",  # roxo
  "HD" = "#277DA1",  # azul
  "PD" = "#E75480"   # rosa
)

shared_colors <- c(
  "HD+PD" = "#8A6FB0",
  "AD+PD" = "#B249EE",
  "AD+HD" = "#5546DD",
  "AD+HD+PD" = "#F7C75F"
)


# ------------------------------------------------------------
# 3. Argumentos
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    paste0(
      "\nUso:\n",
      "Rscript scripts/05_intersection_analysis.R ",
      "<AD_deg_all.csv> ",
      "<HD_deg_all.csv> ",
      "<PD_deg_all.csv> ",
      "<output_dir>\n"
    )
  )
}

ad_file <- args[1]
hd_file <- args[2]
pd_file <- args[3]
output_dir <- args[4]


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
# 5. Ler resultados
# ------------------------------------------------------------

ad <- read_csv(
  ad_file,
  show_col_types = FALSE
)

hd <- read_csv(
  hd_file,
  show_col_types = FALSE
)

pd <- read_csv(
  pd_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 6. Validar colunas
# ------------------------------------------------------------

required_columns <- c(
  "gene_id",
  "symbol",
  "logFC",
  "P.Value",
  "adj.P.Val",
  "significant_primary",
  "direction_primary",
  "significant_article",
  "direction_article"
)

validate_columns <- function(data, dataset) {

  missing_columns <- setdiff(
    required_columns,
    colnames(data)
  )

  if (length(missing_columns) > 0) {
    stop(
      dataset,
      " possui colunas ausentes: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
}

validate_columns(ad, "AD")
validate_columns(hd, "HD")
validate_columns(pd, "PD")


# ------------------------------------------------------------
# 7. Limpar Ensembl
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

ad <- ad %>%
  mutate(
    gene_id = clean_ensembl(gene_id)
  )

hd <- hd %>%
  mutate(
    gene_id = clean_ensembl(gene_id)
  )

pd <- pd %>%
  mutate(
    gene_id = clean_ensembl(gene_id)
  )


# ------------------------------------------------------------
# 8. Remover genes sem Ensembl
# ------------------------------------------------------------

ad_valid <- ad %>%
  filter(
    !is.na(gene_id)
  )

hd_valid <- hd %>%
  filter(
    !is.na(gene_id)
  )

pd_valid <- pd %>%
  filter(
    !is.na(gene_id)
  )


# ------------------------------------------------------------
# 9. Resolver duplicações de Ensembl
#
# Para a integração precisamos de uma linha por Ensembl.
#
# Em caso de duplicação:
#   1. menor adj.P.Val
#   2. maior |logFC| em caso de empate
#
# Isso NÃO modifica a análise diferencial original.
# ------------------------------------------------------------

collapse_by_ensembl <- function(data) {

  data %>%
    mutate(
      abs_logFC = abs(logFC)
    ) %>%
    arrange(
      gene_id,
      adj.P.Val,
      desc(abs_logFC)
    ) %>%
    group_by(
      gene_id
    ) %>%
    slice_head(
      n = 1
    ) %>%
    ungroup() %>%
    select(
      -abs_logFC
    )
}

ad_unique <- collapse_by_ensembl(ad_valid)
hd_unique <- collapse_by_ensembl(hd_valid)
pd_unique <- collapse_by_ensembl(pd_valid)


# ------------------------------------------------------------
# 10. Universos testáveis
# ------------------------------------------------------------

ad_ids <- unique(ad_unique$gene_id)
hd_ids <- unique(hd_unique$gene_id)
pd_ids <- unique(pd_unique$gene_id)

common_universe <- Reduce(
  intersect,
  list(
    ad_ids,
    hd_ids,
    pd_ids
  )
)

common_n <- length(common_universe)


message("")
message("============================================")
message("INTEGRAÇÃO TRANSDIAGNÓSTICA")
message("============================================")
message("AD genes: ", length(ad_ids))
message("HD genes: ", length(hd_ids))
message("PD genes: ", length(pd_ids))
message("Universo comum AD/HD/PD: ", common_n)
message("============================================")
message("")


# ------------------------------------------------------------
# 11. Salvar universo comum
# ------------------------------------------------------------

write_csv(
  tibble(
    gene_id = sort(common_universe)
  ),
  file.path(
    output_dir,
    "common_testable_universe.csv"
  )
)


# ------------------------------------------------------------
# 12. Restringir ao universo comum
# ------------------------------------------------------------

ad_common <- ad_unique %>%
  filter(
    gene_id %in% common_universe
  )

hd_common <- hd_unique %>%
  filter(
    gene_id %in% common_universe
  )

pd_common <- pd_unique %>%
  filter(
    gene_id %in% common_universe
  )


# ------------------------------------------------------------
# 13. Preparar cada dataset para integração
# ------------------------------------------------------------

prepare_for_merge <- function(data, prefix) {

  data %>%
    select(
      gene_id,
      symbol,
      logFC,
      P.Value,
      adj.P.Val,
      significant_primary,
      direction_primary,
      significant_article,
      direction_article
    ) %>%
    rename(
      !!paste0("symbol_", prefix) := symbol,
      !!paste0("logFC_", prefix) := logFC,
      !!paste0("PValue_", prefix) := P.Value,
      !!paste0("FDR_", prefix) := adj.P.Val,
      !!paste0("sig_primary_", prefix) := significant_primary,
      !!paste0("direction_primary_", prefix) := direction_primary,
      !!paste0("sig_article_", prefix) := significant_article,
      !!paste0("direction_article_", prefix) := direction_article
    )
}

ad_merge <- prepare_for_merge(
  ad_common,
  "AD"
)

hd_merge <- prepare_for_merge(
  hd_common,
  "HD"
)

pd_merge <- prepare_for_merge(
  pd_common,
  "PD"
)


# ------------------------------------------------------------
# 14. Matriz integrada
# ------------------------------------------------------------

integrated <- ad_merge %>%
  inner_join(
    hd_merge,
    by = "gene_id"
  ) %>%
  inner_join(
    pd_merge,
    by = "gene_id"
  )

if (nrow(integrated) != common_n) {
  stop(
    "Erro: matriz integrada possui ",
    nrow(integrated),
    " genes, mas o universo comum possui ",
    common_n,
    "."
  )
}


# ------------------------------------------------------------
# 15. Símbolo consenso
#
# Prioridade:
# PD -> HD -> AD
#
# PD possui cobertura de símbolos mais completa.
# ------------------------------------------------------------

integrated <- integrated %>%
  mutate(
    symbol = coalesce(
      symbol_PD,
      symbol_HD,
      symbol_AD
    )
  )


# ------------------------------------------------------------
# 16. Número de doenças com DEG principal
# ------------------------------------------------------------

integrated <- integrated %>%
  mutate(
    n_sig_primary =
      as.integer(sig_primary_AD) +
      as.integer(sig_primary_HD) +
      as.integer(sig_primary_PD)
  )


# ------------------------------------------------------------
# 17. Padrão de significância principal
# ------------------------------------------------------------

integrated <- integrated %>%
  mutate(
    primary_pattern = case_when(

      sig_primary_AD &
        sig_primary_HD &
        sig_primary_PD ~ "AD+HD+PD",

      sig_primary_AD &
        sig_primary_HD ~ "AD+HD",

      sig_primary_AD &
        sig_primary_PD ~ "AD+PD",

      sig_primary_HD &
        sig_primary_PD ~ "HD+PD",

      sig_primary_AD ~ "AD only",

      sig_primary_HD ~ "HD only",

      sig_primary_PD ~ "PD only",

      TRUE ~ "None"
    )
  )


# ------------------------------------------------------------
# 18. Direção do efeito
#
# Considera apenas o sinal do logFC,
# independentemente da significância.
# ------------------------------------------------------------

effect_direction <- function(x) {

  case_when(
    x > 0 ~ "Up",
    x < 0 ~ "Down",
    TRUE ~ "Zero"
  )
}

integrated <- integrated %>%
  mutate(
    effect_AD = effect_direction(logFC_AD),
    effect_HD = effect_direction(logFC_HD),
    effect_PD = effect_direction(logFC_PD)
  )


# ------------------------------------------------------------
# 19. Concordância direcional nas três doenças
# ------------------------------------------------------------

integrated <- integrated %>%
  mutate(
    directional_pattern = paste(
      effect_AD,
      effect_HD,
      effect_PD,
      sep = "/"
    ),

    concordant_all_three = (
      effect_AD == effect_HD &
        effect_HD == effect_PD &
        effect_AD != "Zero"
    )
  )


# ------------------------------------------------------------
# 20. Interseção estrita HD + PD
# ------------------------------------------------------------

shared_hd_pd <- integrated %>%
  filter(
    sig_primary_HD,
    sig_primary_PD
  ) %>%
  mutate(
    hd_pd_same_direction =
      effect_HD == effect_PD
  ) %>%
  arrange(
    FDR_HD + FDR_PD
  )


# ------------------------------------------------------------
# 21. Compartilhados HD + PD na mesma direção
# ------------------------------------------------------------

shared_hd_pd_same_direction <- shared_hd_pd %>%
  filter(
    hd_pd_same_direction
  )


# ------------------------------------------------------------
# 22. Compartilhados HD + PD em direções opostas
# ------------------------------------------------------------

shared_hd_pd_opposite_direction <- shared_hd_pd %>%
  filter(
    !hd_pd_same_direction
  )


# ------------------------------------------------------------
# 23. Genes concordantes nas três doenças
#
# IMPORTANTE:
# isso NÃO significa DEG compartilhado.
# Significa somente mesmo sinal de logFC nas três doenças.
# ------------------------------------------------------------

concordant_three <- integrated %>%
  filter(
    concordant_all_three
  ) %>%
  arrange(
    desc(n_sig_primary)
  )


# ------------------------------------------------------------
# 24. Interseção pelo critério comparativo do artigo
# ------------------------------------------------------------

integrated <- integrated %>%
  mutate(
    n_sig_article =
      as.integer(sig_article_AD) +
      as.integer(sig_article_HD) +
      as.integer(sig_article_PD),

    article_pattern = case_when(

      sig_article_AD &
        sig_article_HD &
        sig_article_PD ~ "AD+HD+PD",

      sig_article_AD &
        sig_article_HD ~ "AD+HD",

      sig_article_AD &
        sig_article_PD ~ "AD+PD",

      sig_article_HD &
        sig_article_PD ~ "HD+PD",

      sig_article_AD ~ "AD only",

      sig_article_HD ~ "HD only",

      sig_article_PD ~ "PD only",

      TRUE ~ "None"
    )
  )


# ------------------------------------------------------------
# 25. Salvar matriz integrada
# ------------------------------------------------------------

write_csv(
  integrated,
  file.path(
    output_dir,
    "integrated_AD_HD_PD.csv"
  )
)


# ------------------------------------------------------------
# 26. Salvar interseções
# ------------------------------------------------------------

write_csv(
  shared_hd_pd,
  file.path(
    output_dir,
    "shared_primary_HD_PD.csv"
  )
)

write_csv(
  shared_hd_pd_same_direction,
  file.path(
    output_dir,
    "shared_primary_HD_PD_same_direction.csv"
  )
)

write_csv(
  shared_hd_pd_opposite_direction,
  file.path(
    output_dir,
    "shared_primary_HD_PD_opposite_direction.csv"
  )
)

write_csv(
  concordant_three,
  file.path(
    output_dir,
    "directionally_concordant_AD_HD_PD.csv"
  )
)


# ------------------------------------------------------------
# 27. Resumo dos padrões principais
# ------------------------------------------------------------

primary_pattern_summary <- integrated %>%
  count(
    primary_pattern,
    name = "n_genes"
  ) %>%
  arrange(
    desc(n_genes)
  )

write_csv(
  primary_pattern_summary,
  file.path(
    output_dir,
    "primary_pattern_summary.csv"
  )
)


# ------------------------------------------------------------
# 28. Resumo do critério comparativo
# ------------------------------------------------------------

article_pattern_summary <- integrated %>%
  count(
    article_pattern,
    name = "n_genes"
  ) %>%
  arrange(
    desc(n_genes)
  )

write_csv(
  article_pattern_summary,
  file.path(
    output_dir,
    "article_pattern_summary.csv"
  )
)


# ------------------------------------------------------------
# 29. Resumo direcional
# ------------------------------------------------------------

directional_summary <- integrated %>%
  count(
    directional_pattern,
    name = "n_genes"
  ) %>%
  arrange(
    desc(n_genes)
  )

write_csv(
  directional_summary,
  file.path(
    output_dir,
    "directional_pattern_summary.csv"
  )
)


# ------------------------------------------------------------
# 30. Gráfico dos padrões principais
# ------------------------------------------------------------

primary_plot_data <- primary_pattern_summary %>%
  filter(
    primary_pattern != "None"
  ) %>%
  mutate(
    plot_color = case_when(
      primary_pattern == "AD only" ~ disease_colors["AD"],
      primary_pattern == "HD only" ~ disease_colors["HD"],
      primary_pattern == "PD only" ~ disease_colors["PD"],
      primary_pattern == "HD+PD" ~ shared_colors["HD+PD"],
      primary_pattern == "AD+PD" ~ shared_colors["AD+PD"],
      primary_pattern == "AD+HD" ~ shared_colors["AD+HD"],
      primary_pattern == "AD+HD+PD" ~ shared_colors["AD+HD+PD"],
      TRUE ~ "grey60"
    )
  )

if (nrow(primary_plot_data) > 0) {

  p_primary <- ggplot(
    primary_plot_data,
    aes(
      x = reorder(
        primary_pattern,
        n_genes
      ),
      y = n_genes,
      fill = plot_color
    )
  ) +
    geom_col(
      width = 0.75
    ) +
    coord_flip() +
    scale_fill_identity() +
    labs(
      title = "Distribuição dos DEGs entre AD, HD e PD",
      subtitle = paste0(
        "Critério principal no universo comum de ",
        format(common_n, big.mark = "."),
        " genes"
      ),
      x = NULL,
      y = "Número de genes"
    ) +
    theme_bw()

  ggsave(
    file.path(
      plot_dir,
      "primary_deg_patterns.png"
    ),
    p_primary,
    width = 8,
    height = 5,
    dpi = 300
  )
}


# ------------------------------------------------------------
# 31. Correlação de logFC
# ------------------------------------------------------------

cor_matrix <- integrated %>%
  select(
    logFC_AD,
    logFC_HD,
    logFC_PD
  ) %>%
  cor(
    use = "pairwise.complete.obs",
    method = "spearman"
  )

cor_df <- as.data.frame(
  cor_matrix
) %>%
  tibble::rownames_to_column(
    "disease"
  )

write_csv(
  cor_df,
  file.path(
    output_dir,
    "logFC_spearman_correlations.csv"
  )
)


# ------------------------------------------------------------
# 32. Scatter HD vs PD
#
# Cinza  = não significativo nos dois
# Azul   = DEG somente em HD
# Rosa   = DEG somente em PD
# Roxo   = DEG compartilhado HD + PD
# ------------------------------------------------------------
scatter_data <- integrated %>%
  mutate(
    scatter_class = case_when(
      sig_primary_HD & sig_primary_PD ~ "HD+PD",
      sig_primary_HD ~ "HD",
      sig_primary_PD ~ "PD",
      TRUE ~ "Other"
    ),
    scatter_class = factor(
      scatter_class,
      levels = c(
        "Other",
        "HD",
        "PD",
        "HD+PD"
      )
    )
  )

scatter_colors <- c(
  "Other" = "grey80",
  "HD" = unname(disease_colors["HD"]),
  "PD" = unname(disease_colors["PD"]),
  "HD+PD" = unname(shared_colors["HD+PD"])
)

p_hd_pd <- ggplot(
  scatter_data,
  aes(
    x = logFC_HD,
    y = logFC_PD,
    color = scatter_class
  )
) +
  geom_point(
    alpha = 0.55,
    size = 1.3
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = scatter_colors,
    drop = FALSE
  ) +
  labs(
    title = "Concordância de efeito — HD vs PD",
    subtitle = paste0(
      "Universo comum de ",
      format(common_n, big.mark = "."),
      " genes"
    ),
    x = "logFC HD",
    y = "logFC PD",
    color = "Classificação"
  ) +
  theme_bw()

ggsave(
  file.path(
    plot_dir,
    "logFC_HD_vs_PD.png"
  ),
  p_hd_pd,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Função para Venn com cores das doenças
# ------------------------------------------------------------

create_venn_plot <- function(
  pattern_data,
  title,
  subtitle
) {

  # Coordenadas dos círculos
  theta <- seq(
    0,
    2 * pi,
    length.out = 500
  )

  circle_data <- bind_rows(

    tibble(
      x = -0.85 + 1.45 * cos(theta),
      y =  0.55 + 1.45 * sin(theta),
      disease = "AD"
    ),

    tibble(
      x =  0.85 + 1.45 * cos(theta),
      y =  0.55 + 1.45 * sin(theta),
      disease = "HD"
    ),

    tibble(
      x =  0.00 + 1.45 * cos(theta),
      y = -0.85 + 1.45 * sin(theta),
      disease = "PD"
    )
  )

  # Garantir que categorias ausentes tenham valor zero
  get_n <- function(pattern) {

    value <- pattern_data %>%
      filter(
        pattern == !!pattern
      ) %>%
      pull(n_genes)

    if (length(value) == 0) {
      return(0)
    }

    value[1]
  }

  counts <- tibble(
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

  ggplot() +

    # Círculos preenchidos
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

    # Contornos
    geom_path(
      data = circle_data,
      aes(
        x = x,
        y = y,
        group = disease,
        color = disease
      ),
      linewidth = 1.2
    ) +

    # Números
    geom_text(
      data = counts,
      aes(
        x = x,
        y = y,
        label = n
      ),
      size = 6
    ) +

    # Nome AD
    annotate(
      "text",
      x = -1.75,
      y = 2.15,
      label = "AD",
      color = unname(disease_colors["AD"]),
      size = 7,
      fontface = "bold"
    ) +

    # Nome HD
    annotate(
      "text",
      x = 1.75,
      y = 2.15,
      label = "HD",
      color = unname(disease_colors["HD"]),
      size = 7,
      fontface = "bold"
    ) +

    # Nome PD
    annotate(
      "text",
      x = 0,
      y = -2.45,
      label = "PD",
      color = unname(disease_colors["PD"]),
      size = 7,
      fontface = "bold"
    ) +

    scale_fill_manual(
      values = disease_colors
    ) +

    scale_color_manual(
      values = disease_colors
    ) +

    coord_fixed(
      xlim = c(-2.5, 2.5),
      ylim = c(-2.7, 2.5),
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
        size = 18,
        face = "bold"
      ),

      plot.subtitle = element_text(
        size = 12
      ),

      plot.margin = margin(
        20,
        20,
        20,
        20
      )
    )
}

# ------------------------------------------------------------
# 33. Venn — critério principal
# ------------------------------------------------------------

primary_venn_data <- primary_pattern_summary %>%
  rename(
    pattern = primary_pattern
  )

p_venn_primary <- create_venn_plot(
  pattern_data = primary_venn_data,
  title = "Interseção dos DEGs entre AD, HD e PD",
  subtitle = paste0(
    "Critério principal no universo comum de ",
    format(common_n, big.mark = "."),
    " genes"
  )
)

ggsave(
  file.path(
    plot_dir,
    "venn_primary_AD_HD_PD.png"
  ),
  p_venn_primary,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)
# ------------------------------------------------------------
# 34. Diagrama de Venn — critério comparativo do artigo
#
# p < 0.05 e |logFC| > 1
#
# Este gráfico é análise de sensibilidade/comparação.
# Não substitui o critério principal.
# ------------------------------------------------------------
article_venn_data <- article_pattern_summary %>%
  rename(
    pattern = article_pattern
  )

p_venn_article <- create_venn_plot(
  pattern_data = article_venn_data,
  title = "Interseção dos genes — critério comparativo",
  subtitle = paste0(
    "p < 0,05 e |logFC| > 1; universo comum de ",
    format(common_n, big.mark = "."),
    " genes"
  )
)

ggsave(
  file.path(
    plot_dir,
    "venn_article_AD_HD_PD.png"
  ),
  p_venn_article,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)
# ------------------------------------------------------------
# 35. Resumo geral
# ------------------------------------------------------------

summary_table <- tibble(

  common_universe =
    nrow(integrated),

  primary_AD =
    sum(
      integrated$sig_primary_AD
    ),

  primary_HD =
    sum(
      integrated$sig_primary_HD
    ),

  primary_PD =
    sum(
      integrated$sig_primary_PD
    ),

  primary_AD_HD_PD =
    sum(
      integrated$primary_pattern ==
        "AD+HD+PD"
    ),

  primary_AD_HD =
    sum(
      integrated$primary_pattern ==
        "AD+HD"
    ),

  primary_AD_PD =
    sum(
      integrated$primary_pattern ==
        "AD+PD"
    ),

  primary_HD_PD =
    sum(
      integrated$primary_pattern ==
        "HD+PD"
    ),

  primary_AD_only =
    sum(
      integrated$primary_pattern ==
        "AD only"
    ),

  primary_HD_only =
    sum(
      integrated$primary_pattern ==
        "HD only"
    ),

  primary_PD_only =
    sum(
      integrated$primary_pattern ==
        "PD only"
    ),

  shared_HD_PD_same_direction =
    nrow(
      shared_hd_pd_same_direction
    ),

  shared_HD_PD_opposite_direction =
    nrow(
      shared_hd_pd_opposite_direction
    ),

  directionally_concordant_all_three =
    nrow(
      concordant_three
    )
)

write_csv(
  summary_table,
  file.path(
    output_dir,
    "intersection_summary.csv"
  )
)


# ------------------------------------------------------------
# 36. Terminal
# ------------------------------------------------------------

message("")
message("============================================")
message("INTERSEÇÃO FINALIZADA")
message("============================================")

message(
  "Universo comum: ",
  nrow(integrated)
)

message("")
message("DEGs principais no universo comum:")

message(
  "  AD: ",
  sum(integrated$sig_primary_AD)
)

message(
  "  HD: ",
  sum(integrated$sig_primary_HD)
)

message(
  "  PD: ",
  sum(integrated$sig_primary_PD)
)

message("")

message(
  "AD + HD + PD: ",
  sum(
    integrated$primary_pattern ==
      "AD+HD+PD"
  )
)

message(
  "AD + HD: ",
  sum(
    integrated$primary_pattern ==
      "AD+HD"
  )
)

message(
  "AD + PD: ",
  sum(
    integrated$primary_pattern ==
      "AD+PD"
  )
)

message(
  "HD + PD: ",
  sum(
    integrated$primary_pattern ==
      "HD+PD"
  )
)

message(
  "AD only: ",
  sum(
    integrated$primary_pattern ==
      "AD only"
  )
)

message(
  "HD only: ",
  sum(
    integrated$primary_pattern ==
      "HD only"
  )
)

message(
  "PD only: ",
  sum(
    integrated$primary_pattern ==
      "PD only"
  )
)

message("")

message(
  "HD + PD mesma direção: ",
  nrow(
    shared_hd_pd_same_direction
  )
)

message(
  "HD + PD direção oposta: ",
  nrow(
    shared_hd_pd_opposite_direction
  )
)

message("")

message(
  "Mesma direção AD/HD/PD (independente de FDR): ",
  nrow(
    concordant_three
  )
)

message("")
message("Critério comparativo:")

message(
  "  AD + HD + PD: ",
  sum(
    integrated$article_pattern ==
      "AD+HD+PD"
  )
)

message(
  "  AD + HD: ",
  sum(
    integrated$article_pattern ==
      "AD+HD"
  )
)

message(
  "  AD + PD: ",
  sum(
    integrated$article_pattern ==
      "AD+PD"
  )
)

message(
  "  HD + PD: ",
  sum(
    integrated$article_pattern ==
      "HD+PD"
  )
)

message("")
message(
  "Resultados: ",
  output_dir
)

message("============================================")