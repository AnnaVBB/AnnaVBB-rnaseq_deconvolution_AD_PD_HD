#!/usr/bin/env Rscript

# ==============================================================================
# 04_qc_counts_GSE125583.R
#
# FASE 01b — QC, PRÉ-FILTRO PARA EXPRESSÃO DIFERENCIAL E PCA
# DATASET CANDIDATO DE AD: GSE125583
#
# OBJETIVO GERAL
# ------------------------------------------------------------------------------
# Avaliar a qualidade da matriz preparada na etapa 02 e definir o conjunto de
# genes que será utilizado na análise de expressão diferencial (DESeq2).
#
# Esta etapa é a primeira em que o universo gênico é reduzido com base em
# expressão.
#
#
# REGRA DE PRÉ-FILTRO
# ------------------------------------------------------------------------------
# Mantemos a mesma regra metodológica utilizada na análise original:
#
#   count >= 10 em pelo menos o número de amostras do menor grupo
#
# Para o GSE125583:
#
#   AD      = 219 amostras
#   Control = 70 amostras
#
# Portanto:
#
#   manter gene se count >= 10 em pelo menos 70 amostras
#
# Essa regra é aplicada ao conjunto total de 289 amostras.
#
# IMPORTANTE:
# - o filtro NÃO exige que essas 70 amostras pertençam ao grupo Controle;
# - o número 70 é usado apenas por ser o tamanho do menor grupo;
# - o objetivo é remover genes com expressão muito baixa/esparsa antes do DESeq2.
#
#
# O QUE ESTE SCRIPT FAZ
# ------------------------------------------------------------------------------
# 1. Lê:
#      - counts preparados na etapa 02;
#      - metadata preparado na etapa 01.
# 2. Valida novamente counts <-> metadata.
# 3. Calcula QC por amostra:
#      - library size;
#      - genes detectados (>0);
#      - genes com count >= 10;
#      - proporção de zeros.
# 4. Calcula QC por gene:
#      - soma total;
#      - número de amostras com count > 0;
#      - número de amostras com count >= 10;
#      - gene totalmente zerado;
#      - decisão keep_for_DE.
# 5. Aplica o pré-filtro:
#
#      count >= 10 em >= 70 amostras
#
# 6. Salva a matriz filtrada para DESeq2.
# 7. Aplica VST (Variance Stabilizing Transformation) apenas para QC/PCA.
# 8. Executa PCA sobre os genes filtrados.
# 9. Gera gráficos:
#      - library size por grupo;
#      - genes detectados por grupo;
#      - PCA por grupo;
#      - PCA por sexo;
#      - PCA por idade.
# 10. Gera relatório textual com todas as auditorias.
#
#
# O QUE ESTE SCRIPT NÃO FAZ
# ------------------------------------------------------------------------------
# - NÃO executa teste de expressão diferencial.
# - NÃO usa VST para DESeq2.
# - NÃO normaliza a matriz salva para DE.
# - NÃO exclui amostras automaticamente com base na PCA.
# - NÃO redefine os grupos usando Braak.
# - NÃO ajusta por idade ou sexo.
#
# A matriz salva para o DESeq2 permanece em RAW COUNTS inteiros.
#
#
# ENTRADAS
# ------------------------------------------------------------------------------
# results/candidate_GSE125583/counts/counts_GSE125583.tsv
#
# results/candidate_GSE125583/metadata/metadata_GSE125583.tsv
#
#
# SAÍDAS
# ------------------------------------------------------------------------------
# results/candidate_GSE125583/qc/
# ├── counts_GSE125583_DE_filtered.tsv
# ├── gene_filter_summary.tsv
# ├── sample_qc_summary.tsv
# └── 04_qc_counts_GSE125583.txt
#
# results/candidate_GSE125583/plots/qc/
# ├── library_size_GSE125583.png
# ├── detected_genes_GSE125583.png
# ├── PCA_GSE125583_group.png
# ├── PCA_GSE125583_sex.png
# └── PCA_GSE125583_age.png
# ==============================================================================


# ==============================================================================
# 1. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table",
  "DESeq2",
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
    "Pacotes necessários não instalados: ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
})


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

counts_file <- file.path(
  "results",
  "candidate_GSE125583",
  "counts",
  "counts_GSE125583.tsv"
)

metadata_file <- file.path(
  "results",
  "candidate_GSE125583",
  "metadata",
  "metadata_GSE125583.tsv"
)

output_dir <- file.path(
  "results",
  "candidate_GSE125583",
  "qc"
)

plot_dir <- file.path(
  "results",
  "candidate_GSE125583",
  "plots",
  "qc"
)

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

filtered_counts_file <- file.path(
  output_dir,
  "counts_GSE125583_DE_filtered.tsv"
)

gene_filter_file <- file.path(
  output_dir,
  "gene_filter_summary.tsv"
)

sample_qc_file <- file.path(
  output_dir,
  "sample_qc_summary.tsv"
)

audit_file <- file.path(
  output_dir,
  "04_qc_counts_GSE125583.txt"
)

library_size_plot_file <- file.path(
  plot_dir,
  "library_size_GSE125583.png"
)

detected_genes_plot_file <- file.path(
  plot_dir,
  "detected_genes_GSE125583.png"
)

pca_group_plot_file <- file.path(
  plot_dir,
  "PCA_GSE125583_group.png"
)

pca_sex_plot_file <- file.path(
  plot_dir,
  "PCA_GSE125583_sex.png"
)

pca_age_plot_file <- file.path(
  plot_dir,
  "PCA_GSE125583_age.png"
)


# ==============================================================================
# 3. CONSTANTES / REGRAS ANALÍTICAS
# ==============================================================================

EXPECTED_N_GENES <- 39376L
EXPECTED_N_SAMPLES <- 289L
EXPECTED_N_AD <- 219L
EXPECTED_N_CONTROL <- 70L

# Critério de pré-filtro aprovado:
MIN_COUNT <- 10L

# O menor grupo possui 70 amostras.
MIN_SAMPLES <- min(
  EXPECTED_N_AD,
  EXPECTED_N_CONTROL
)

# Paleta usada apenas para visualização.
GROUP_COLORS <- c(
  "Control" = "#7A7A7A",
  "AD" = "#E67E22"
)

SEX_COLORS <- c(
  "F" = "#8E44AD",
  "M" = "#2980B9",
  "Missing" = "#9E9E9E"
)


# ==============================================================================
# 4. FUNÇÕES AUXILIARES
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 Validar arquivo
# ------------------------------------------------------------------------------

validate_file <- function(path) {

  if (!file.exists(path)) {
    stop(
      "Arquivo não encontrado: ",
      path
    )
  }

  info <- file.info(path)

  if (is.na(info$size) || info$size == 0) {
    stop(
      "Arquivo vazio ou tamanho indisponível: ",
      path
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 4.2 Validar GeneID
# ------------------------------------------------------------------------------

validate_gene_id <- function(
  counts,
  dataset_name
) {

  if (!"GeneID" %in% names(counts)) {
    stop(
      dataset_name,
      ": coluna GeneID não encontrada."
    )
  }

  if (names(counts)[1] != "GeneID") {
    stop(
      dataset_name,
      ": GeneID deve ser a primeira coluna."
    )
  }

  gene_id <- as.character(
    counts$GeneID
  )

  missing_gene_id <- (
    is.na(gene_id) |
      trimws(gene_id) == ""
  )

  if (any(missing_gene_id)) {
    stop(
      dataset_name,
      ": existem GeneIDs ausentes."
    )
  }

  if (anyDuplicated(gene_id)) {
    stop(
      dataset_name,
      ": existem GeneIDs duplicados."
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 4.3 Validar metadata
# ------------------------------------------------------------------------------

validate_metadata <- function(metadata) {

  required_columns <- c(
    "sample_id",
    "group",
    "age",
    "sex",
    "braak_score"
  )

  missing_columns <- setdiff(
    required_columns,
    names(metadata)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Metadata sem colunas obrigatórias: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  if (any(is.na(metadata$sample_id))) {
    stop(
      "Existem sample IDs ausentes no metadata."
    )
  }

  if (anyDuplicated(metadata$sample_id)) {
    stop(
      "Existem sample IDs duplicados no metadata."
    )
  }

  if (any(is.na(metadata$group))) {
    stop(
      "Existem grupos ausentes no metadata."
    )
  }

  unexpected_groups <- setdiff(
    unique(
      metadata$group
    ),
    c(
      "Control",
      "AD"
    )
  )

  if (length(unexpected_groups) > 0) {
    stop(
      "Grupos inesperados no metadata: ",
      paste(
        unexpected_groups,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 4.4 Validar correspondência counts <-> metadata
# ------------------------------------------------------------------------------

validate_sample_mapping <- function(
  counts,
  metadata,
  dataset_name
) {

  count_samples <- setdiff(
    names(counts),
    "GeneID"
  )

  metadata_samples <- metadata$sample_id

  if (anyDuplicated(count_samples)) {
    stop(
      dataset_name,
      ": existem amostras duplicadas nos counts."
    )
  }

  only_counts <- setdiff(
    count_samples,
    metadata_samples
  )

  only_metadata <- setdiff(
    metadata_samples,
    count_samples
  )

  if (length(only_counts) > 0) {
    stop(
      dataset_name,
      ": amostras somente nos counts: ",
      paste(
        only_counts,
        collapse = ", "
      )
    )
  }

  if (length(only_metadata) > 0) {
    stop(
      dataset_name,
      ": amostras somente no metadata: ",
      paste(
        only_metadata,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 4.5 Formatação pt-BR
# ------------------------------------------------------------------------------

format_int <- function(x) {

  format(
    x,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


format_num <- function(
  x,
  digits = 2
) {

  format(
    round(
      x,
      digits
    ),
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


# ------------------------------------------------------------------------------
# 4.6 Tema gráfico
# ------------------------------------------------------------------------------

qc_theme <- function() {

  theme_bw(
    base_size = 12
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 15
      ),
      plot.subtitle = element_text(
        size = 11
      ),
      axis.title = element_text(
        face = "bold"
      ),
      legend.title = element_text(
        face = "bold"
      )
    )
}


# ==============================================================================
# 5. VALIDAR ENTRADAS
# ==============================================================================

validate_file(
  counts_file
)

validate_file(
  metadata_file
)


# ==============================================================================
# 6. LER METADATA
# ==============================================================================

metadata <- fread(
  metadata_file,
  na.strings = c(
    "NA",
    ""
  ),
  showProgress = FALSE
)

validate_metadata(
  metadata
)

if (nrow(metadata) != EXPECTED_N_SAMPLES) {
  stop(
    "Número inesperado de amostras no metadata. Esperado: ",
    EXPECTED_N_SAMPLES,
    "; observado: ",
    nrow(metadata),
    "."
  )
}

n_ad <- metadata[
  group == "AD",
  .N
]

n_control <- metadata[
  group == "Control",
  .N
]

if (n_ad != EXPECTED_N_AD) {
  stop(
    "Número inesperado de AD. Esperado: ",
    EXPECTED_N_AD,
    "; observado: ",
    n_ad,
    "."
  )
}

if (n_control != EXPECTED_N_CONTROL) {
  stop(
    "Número inesperado de controles. Esperado: ",
    EXPECTED_N_CONTROL,
    "; observado: ",
    n_control,
    "."
  )
}


# ==============================================================================
# 7. LER MATRIZ DE COUNTS
# ==============================================================================

counts <- fread(
  counts_file,
  showProgress = TRUE
)

counts[
  ,
  GeneID := as.character(
    GeneID
  )
]

validate_gene_id(
  counts,
  "GSE125583"
)

validate_sample_mapping(
  counts,
  metadata,
  "GSE125583"
)

if (nrow(counts) != EXPECTED_N_GENES) {
  stop(
    "Número inesperado de genes. Esperado: ",
    EXPECTED_N_GENES,
    "; observado: ",
    nrow(counts),
    "."
  )
}

sample_ids <- setdiff(
  names(counts),
  "GeneID"
)

if (length(sample_ids) != EXPECTED_N_SAMPLES) {
  stop(
    "Número inesperado de amostras na matriz."
  )
}


# ==============================================================================
# 8. REORDENAR COUNTS PARA A ORDEM DO METADATA
# ==============================================================================
#
# Mesmo que as etapas anteriores já tenham mostrado que a ordem coincide,
# reordenamos explicitamente por segurança.

counts <- counts[
  ,
  c(
    "GeneID",
    metadata$sample_id
  ),
  with = FALSE
]

sample_ids <- metadata$sample_id

if (!identical(
  setdiff(
    names(counts),
    "GeneID"
  ),
  metadata$sample_id
)) {
  stop(
    "Falha ao alinhar counts e metadata."
  )
}


# ==============================================================================
# 9. TRANSFORMAR COUNTS EM MATRIZ NUMÉRICA
# ==============================================================================

count_matrix <- as.matrix(
  counts[
    ,
    ..sample_ids
  ]
)

rownames(
  count_matrix
) <- counts$GeneID


# ==============================================================================
# 10. VALIDAR INTEGRIDADE DA MATRIZ
# ==============================================================================

if (anyNA(count_matrix)) {
  stop(
    "Existem counts ausentes."
  )
}

if (any(!is.finite(count_matrix))) {
  stop(
    "Existem counts não finitos."
  )
}

if (any(count_matrix < 0)) {
  stop(
    "Existem counts negativos."
  )
}

fractional_values <- abs(
  count_matrix -
    round(
      count_matrix
    )
) > 1e-8

if (any(fractional_values)) {
  stop(
    "Existem counts fracionários; esperavam-se raw counts inteiros."
  )
}

# Garantir storage mode inteiro para o DESeq2.
storage.mode(
  count_matrix
) <- "integer"


# ==============================================================================
# 11. QC POR AMOSTRA — ANTES DO FILTRO
# ==============================================================================

library_sizes <- colSums(
  count_matrix
)

detected_genes <- colSums(
  count_matrix > 0
)

genes_count_ge_10 <- colSums(
  count_matrix >= MIN_COUNT
)

zero_fraction_sample <- colMeans(
  count_matrix == 0
)

sample_qc <- data.table(
  sample_id = sample_ids,
  group = metadata$group,
  age = metadata$age,
  sex = metadata$sex,
  braak_score = metadata$braak_score,
  library_size = as.numeric(
    library_sizes[
      sample_ids
    ]
  ),
  detected_genes = as.integer(
    detected_genes[
      sample_ids
    ]
  ),
  genes_count_ge_10 = as.integer(
    genes_count_ge_10[
      sample_ids
    ]
  ),
  zero_fraction = as.numeric(
    zero_fraction_sample[
      sample_ids
    ]
  )
)


# ==============================================================================
# 12. QC POR GENE E DEFINIÇÃO DO PRÉ-FILTRO
# ==============================================================================

gene_total_counts <- rowSums(
  count_matrix
)

gene_n_detected <- rowSums(
  count_matrix > 0
)

gene_n_ge_10 <- rowSums(
  count_matrix >= MIN_COUNT
)

gene_all_zero <- (
  gene_total_counts == 0
)

keep_for_de <- (
  gene_n_ge_10 >= MIN_SAMPLES
)

gene_filter_summary <- data.table(
  GeneID = rownames(
    count_matrix
  ),
  total_count = as.numeric(
    gene_total_counts
  ),
  n_samples_detected = as.integer(
    gene_n_detected
  ),
  n_samples_count_ge_10 = as.integer(
    gene_n_ge_10
  ),
  all_zero = gene_all_zero,
  keep_for_DE = keep_for_de
)


# ==============================================================================
# 13. RESUMO DO FILTRO
# ==============================================================================

n_genes_before <- nrow(
  count_matrix
)

n_all_zero <- sum(
  gene_all_zero
)

n_genes_kept <- sum(
  keep_for_de
)

n_genes_removed <- (
  n_genes_before -
    n_genes_kept
)

if (n_genes_kept == 0) {
  stop(
    "O pré-filtro removeu todos os genes."
  )
}


# ==============================================================================
# 14. CONSTRUIR MATRIZ FILTRADA PARA DESeq2
# ==============================================================================
#
# IMPORTANTE:
# Esta matriz continua em RAW COUNTS.
#
# A única alteração é a remoção de genes que não satisfazem:
#
#   count >= 10 em pelo menos 70 amostras.

filtered_matrix <- count_matrix[
  keep_for_de,
  ,
  drop = FALSE
]

filtered_counts <- data.table(
  GeneID = rownames(
    filtered_matrix
  )
)

filtered_counts <- cbind(
  filtered_counts,
  as.data.table(
    filtered_matrix
  )
)

setnames(
  filtered_counts,
  old = names(
    filtered_counts
  )[-1],
  new = sample_ids
)


# ==============================================================================
# 15. VALIDAR MATRIZ FILTRADA
# ==============================================================================

if (ncol(filtered_matrix) != EXPECTED_N_SAMPLES) {
  stop(
    "O número de amostras mudou após o filtro."
  )
}

if (!identical(
  colnames(filtered_matrix),
  metadata$sample_id
)) {
  stop(
    "A ordem das amostras mudou após o filtro."
  )
}

if (anyNA(filtered_matrix)) {
  stop(
    "A matriz filtrada contém valores ausentes."
  )
}

if (any(filtered_matrix < 0)) {
  stop(
    "A matriz filtrada contém valores negativos."
  )
}

if (any(
  abs(
    filtered_matrix -
      round(
        filtered_matrix
      )
  ) > 1e-8
)) {
  stop(
    "A matriz filtrada contém valores fracionários."
  )
}


# ==============================================================================
# 16. SALVAR RESULTADOS DO PRÉ-FILTRO
# ==============================================================================

fwrite(
  filtered_counts,
  filtered_counts_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  gene_filter_summary,
  gene_filter_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  sample_qc,
  sample_qc_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 17. PREPARAR DESeqDataSet APENAS PARA VST / PCA
# ==============================================================================
#
# O DESeqDataSet criado aqui NÃO é utilizado para testar expressão diferencial.
# Ele serve apenas para obter a transformação VST usada no QC.
#
# Usamos design ~ group apenas para construir o objeto.
# Como blind = TRUE abaixo, a transformação VST não utiliza o desenho
# experimental para estimar a estrutura da PCA.

col_data <- as.data.frame(
  metadata
)

rownames(
  col_data
) <- col_data$sample_id

# Garantir referência e ordem dos níveis.
col_data$group <- factor(
  col_data$group,
  levels = c(
    "Control",
    "AD"
  )
)

dds_qc <- DESeqDataSetFromMatrix(
  countData = filtered_matrix,
  colData = col_data,
  design = ~ group
)


# ==============================================================================
# 18. VST PARA QC
# ==============================================================================
#
# blind = TRUE:
#   a transformação é estimada sem usar o efeito de grupo para "proteger"
#   diferenças entre AD e Control.
#
# Para PCA exploratória de QC, essa é uma escolha conservadora.
#
# Se vst() não puder usar o atalho por alguma razão, varianceStabilizingTransformation
# é usada como fallback.

vsd <- tryCatch(
  {
    vst(
      dds_qc,
      blind = TRUE
    )
  },
  error = function(e) {

    message(
      "vst() não pôde ser usado diretamente; tentando ",
      "varianceStabilizingTransformation()."
    )

    varianceStabilizingTransformation(
      dds_qc,
      blind = TRUE
    )
  }
)


# ==============================================================================
# 19. PCA
# ==============================================================================

pca_fit <- prcomp(
  t(
    assay(
      vsd
    )
  ),
  center = TRUE,
  scale. = FALSE
)

pca_variance <- (
  pca_fit$sdev^2 /
    sum(
      pca_fit$sdev^2
    )
) * 100

pc1_percent <- pca_variance[1]
pc2_percent <- pca_variance[2]

pca_scores <- data.table(
  sample_id = rownames(
    pca_fit$x
  ),
  PC1 = pca_fit$x[
    ,
    1
  ],
  PC2 = pca_fit$x[
    ,
    2
  ]
)

pca_scores <- merge(
  pca_scores,
  metadata[
    ,
    .(
      sample_id,
      group,
      age,
      sex,
      braak_score
    )
  ],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

# Restaurar a ordem do metadata depois do merge.
pca_scores <- pca_scores[
  match(
    metadata$sample_id,
    sample_id
  )
]

if (!identical(
  pca_scores$sample_id,
  metadata$sample_id
)) {
  stop(
    "Falha ao alinhar PCA e metadata."
  )
}

# Para o gráfico por sexo, mantemos explicitamente o valor ausente.
pca_scores[
  ,
  sex_plot := fifelse(
    is.na(sex),
    "Missing",
    sex
  )
]


# ==============================================================================
# 20. GRÁFICO — LIBRARY SIZE
# ==============================================================================

p_library <- ggplot(
  sample_qc,
  aes(
    x = group,
    y = library_size,
    fill = group
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.6
  ) +
  geom_jitter(
    width = 0.18,
    alpha = 0.45,
    size = 1.4
  ) +
  scale_fill_manual(
    values = GROUP_COLORS
  ) +
  scale_y_log10(
    labels = scales::label_number(
      big.mark = ".",
      decimal.mark = ","
    )
  ) +
  labs(
    title = "Library size — GSE125583",
    subtitle = "Raw counts antes do pré-filtro",
    x = NULL,
    y = "Library size (escala log10)",
    fill = "Grupo"
  ) +
  qc_theme() +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = library_size_plot_file,
  plot = p_library,
  width = 7,
  height = 5,
  dpi = 300
)


# ==============================================================================
# 21. GRÁFICO — GENES DETECTADOS POR AMOSTRA
# ==============================================================================

p_detected <- ggplot(
  sample_qc,
  aes(
    x = group,
    y = detected_genes,
    fill = group
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.6
  ) +
  geom_jitter(
    width = 0.18,
    alpha = 0.45,
    size = 1.4
  ) +
  scale_fill_manual(
    values = GROUP_COLORS
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      big.mark = ".",
      decimal.mark = ","
    )
  ) +
  labs(
    title = "Genes detectados por amostra — GSE125583",
    subtitle = "Gene detectado = count > 0",
    x = NULL,
    y = "Número de genes detectados",
    fill = "Grupo"
  ) +
  qc_theme() +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = detected_genes_plot_file,
  plot = p_detected,
  width = 7,
  height = 5,
  dpi = 300
)


# ==============================================================================
# 22. PCA — COLORIDA POR GRUPO
# ==============================================================================

p_pca_group <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = group
  )
) +
  geom_point(
    size = 2.4,
    alpha = 0.75
  ) +
  scale_color_manual(
    values = GROUP_COLORS
  ) +
  labs(
    title = "PCA — GSE125583",
    subtitle = paste0(
      "VST, blind = TRUE; ",
      format_num(
        n_genes_kept,
        digits = 0
      ),
      " genes após pré-filtro"
    ),
    x = paste0(
      "PC1 (",
      format_num(
        pc1_percent,
        digits = 1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      format_num(
        pc2_percent,
        digits = 1
      ),
      "%)"
    ),
    color = "Grupo"
  ) +
  qc_theme()

ggsave(
  filename = pca_group_plot_file,
  plot = p_pca_group,
  width = 7.5,
  height = 6,
  dpi = 300
)


# ==============================================================================
# 23. PCA — COLORIDA POR SEXO
# ==============================================================================

p_pca_sex <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = sex_plot
  )
) +
  geom_point(
    size = 2.4,
    alpha = 0.75
  ) +
  scale_color_manual(
    values = SEX_COLORS
  ) +
  labs(
    title = "PCA — GSE125583 por sexo",
    subtitle = "Exploração de possível estrutura associada ao sexo",
    x = paste0(
      "PC1 (",
      format_num(
        pc1_percent,
        digits = 1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      format_num(
        pc2_percent,
        digits = 1
      ),
      "%)"
    ),
    color = "Sexo"
  ) +
  qc_theme()

ggsave(
  filename = pca_sex_plot_file,
  plot = p_pca_sex,
  width = 7.5,
  height = 6,
  dpi = 300
)


# ==============================================================================
# 24. PCA — COLORIDA POR IDADE
# ==============================================================================

p_pca_age <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = age
  )
) +
  geom_point(
    size = 2.4,
    alpha = 0.8
  ) +
  scale_color_viridis_c(
    option = "viridis",
    na.value = "#9E9E9E"
  ) +
  labs(
    title = "PCA — GSE125583 por idade",
    subtitle = "Exploração de possível estrutura associada à idade",
    x = paste0(
      "PC1 (",
      format_num(
        pc1_percent,
        digits = 1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      format_num(
        pc2_percent,
        digits = 1
      ),
      "%)"
    ),
    color = "Idade"
  ) +
  qc_theme()

ggsave(
  filename = pca_age_plot_file,
  plot = p_pca_age,
  width = 7.5,
  height = 6,
  dpi = 300
)


# ==============================================================================
# 25. RESUMOS DESCRITIVOS POR GRUPO
# ==============================================================================

sample_summary_by_group <- sample_qc[
  ,
  .(
    n = .N,

    library_size_min = min(
      library_size
    ),

    library_size_median = median(
      library_size
    ),

    library_size_max = max(
      library_size
    ),

    detected_genes_min = min(
      detected_genes
    ),

    detected_genes_median = median(
      detected_genes
    ),

    detected_genes_max = max(
      detected_genes
    ),

    zero_fraction_median = median(
      zero_fraction
    )
  ),
  by = group
][
  match(
    c(
      "Control",
      "AD"
    ),
    group
  )
]


# ==============================================================================
# 26. RELATÓRIO TEXTUAL
# ==============================================================================

sink(
  audit_file
)

cat("============================================================\n")
cat("04 — QC E PRÉ-FILTRO GSE125583\n")
cat("============================================================\n\n")


cat("------------------------------------------------------------\n")
cat("1. ENTRADA\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes antes do filtro: ",
  format_int(
    n_genes_before
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras: ",
  format_int(
    ncol(
      count_matrix
    )
  ),
  "\n",
  sep = ""
)

cat(
  "AD: ",
  format_int(
    n_ad
  ),
  "\n",
  sep = ""
)

cat(
  "Control: ",
  format_int(
    n_control
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("2. REGRA DE PRÉ-FILTRO\n")
cat("------------------------------------------------------------\n")

cat(
  "Count mínimo: ",
  MIN_COUNT,
  "\n",
  sep = ""
)

cat(
  "Número mínimo de amostras: ",
  MIN_SAMPLES,
  "\n",
  sep = ""
)

cat(
  "Regra: manter gene se count >= ",
  MIN_COUNT,
  " em pelo menos ",
  MIN_SAMPLES,
  " amostras.\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("3. RESULTADO DO PRÉ-FILTRO\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes totalmente zerados antes do filtro: ",
  format_int(
    n_all_zero
  ),
  "\n",
  sep = ""
)

cat(
  "Genes mantidos para DESeq2: ",
  format_int(
    n_genes_kept
  ),
  "\n",
  sep = ""
)

cat(
  "Genes removidos: ",
  format_int(
    n_genes_removed
  ),
  "\n",
  sep = ""
)

cat(
  "Percentual mantido: ",
  format_num(
    100 *
      n_genes_kept /
      n_genes_before,
    digits = 2
  ),
  "%\n",
  sep = ""
)

cat(
  "Percentual removido: ",
  format_num(
    100 *
      n_genes_removed /
      n_genes_before,
    digits = 2
  ),
  "%\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("4. QC POR AMOSTRA — CONJUNTO COMPLETO\n")
cat("------------------------------------------------------------\n")

cat(
  "Library size — mínimo: ",
  format_int(
    min(
      sample_qc$library_size
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Library size — mediana: ",
  format_int(
    median(
      sample_qc$library_size
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Library size — máximo: ",
  format_int(
    max(
      sample_qc$library_size
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes detectados — mínimo: ",
  format_int(
    min(
      sample_qc$detected_genes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes detectados — mediana: ",
  format_int(
    median(
      sample_qc$detected_genes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes detectados — máximo: ",
  format_int(
    max(
      sample_qc$detected_genes
    )
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("5. QC POR GRUPO\n")
cat("------------------------------------------------------------\n")

print(
  sample_summary_by_group
)

cat("\n")


cat("------------------------------------------------------------\n")
cat("6. PCA — VST\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes utilizados na PCA: ",
  format_int(
    n_genes_kept
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras utilizadas na PCA: ",
  format_int(
    nrow(
      pca_scores
    )
  ),
  "\n",
  sep = ""
)

cat(
  "PC1: ",
  format_num(
    pc1_percent,
    digits = 2
  ),
  "% da variância\n",
  sep = ""
)

cat(
  "PC2: ",
  format_num(
    pc2_percent,
    digits = 2
  ),
  "% da variância\n",
  sep = ""
)

cat(
  "VST: blind = TRUE\n\n"
)


cat("------------------------------------------------------------\n")
cat("7. AMOSTRAS E EXCLUSÕES\n")
cat("------------------------------------------------------------\n")

cat(
  "Amostras antes do QC: ",
  format_int(
    EXPECTED_N_SAMPLES
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras removidas automaticamente: 0\n"
)

cat(
  "Amostras mantidas: ",
  format_int(
    EXPECTED_N_SAMPLES
  ),
  "\n\n",
  sep = ""
)

cat(
  "IMPORTANTE: a PCA e os gráficos de QC são exploratórios. ",
  "Nenhuma amostra é excluída automaticamente com base em posição na PCA, ",
  "library size ou número de genes detectados.\n\n"
)


cat("------------------------------------------------------------\n")
cat("8. ARQUIVOS GERADOS\n")
cat("------------------------------------------------------------\n")

cat(
  "Counts filtrados: ",
  filtered_counts_file,
  "\n",
  sep = ""
)

cat(
  "Resumo por gene: ",
  gene_filter_file,
  "\n",
  sep = ""
)

cat(
  "Resumo por amostra: ",
  sample_qc_file,
  "\n",
  sep = ""
)

cat(
  "Plot library size: ",
  library_size_plot_file,
  "\n",
  sep = ""
)

cat(
  "Plot genes detectados: ",
  detected_genes_plot_file,
  "\n",
  sep = ""
)

cat(
  "PCA por grupo: ",
  pca_group_plot_file,
  "\n",
  sep = ""
)

cat(
  "PCA por sexo: ",
  pca_sex_plot_file,
  "\n",
  sep = ""
)

cat(
  "PCA por idade: ",
  pca_age_plot_file,
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("9. INTERPRETAÇÃO DA ETAPA\n")
cat("------------------------------------------------------------\n")

cat(
  "- O universo bruto de 39.376 genes foi submetido ao pré-filtro de expressão.\n"
)

cat(
  "- O critério aplicado foi count >= 10 em pelo menos 70 amostras.\n"
)

cat(
  "- A matriz filtrada permanece em raw counts inteiros e será a entrada do DESeq2.\n"
)

cat(
  "- A transformação VST foi utilizada exclusivamente para QC/PCA.\n"
)

cat(
  "- Nenhuma amostra foi removida automaticamente.\n"
)

cat(
  "- O número de genes retidos nesta etapa define o universo candidato de genes testáveis para o GSE125583.\n"
)

cat(
  "- A interseção desse universo com os universos filtrados de HD e PD será avaliada na integração posterior.\n"
)

cat("\n")
cat("============================================================\n")
cat("04_qc_counts_GSE125583.R CONCLUÍDO\n")
cat("============================================================\n")

sink()


# ==============================================================================
# 27. RESUMO NO TERMINAL
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("04 — QC GSE125583 CONCLUÍDO\n")
cat("============================================================\n")

cat(
  "Genes antes: ",
  format_int(
    n_genes_before
  ),
  "\n",
  sep = ""
)

cat(
  "Genes após pré-filtro: ",
  format_int(
    n_genes_kept
  ),
  "\n",
  sep = ""
)

cat(
  "Genes removidos: ",
  format_int(
    n_genes_removed
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras: ",
  format_int(
    EXPECTED_N_SAMPLES
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras removidas: 0\n"
)

cat(
  "PC1: ",
  format_num(
    pc1_percent,
    digits = 1
  ),
  "%\n",
  sep = ""
)

cat(
  "PC2: ",
  format_num(
    pc2_percent,
    digits = 1
  ),
  "%\n",
  sep = ""
)

cat(
  "Counts filtrados: ",
  filtered_counts_file,
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
