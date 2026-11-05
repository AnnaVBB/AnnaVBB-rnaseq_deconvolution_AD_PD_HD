# ==============================================================================
# 07_compare_AD_candidates_ROSMAP_GSE125583.R
# ==============================================================================
#
# Fase 01c — comparação formal entre as duas opções de Alzheimer
#
# Objetivo
# --------
# Comparar:
#
#   1) GSE125583 — giro fusiforme
#   2) ROSMAP    — DLPFC
#
# de acordo com os critérios definidos ANTES da inspeção dos resultados do
# ROSMAP:
#
#   A. concordância gene-a-gene dos log2FoldChanges;
#   B. coerência funcional do enriquecimento;
#   C. estabilidade frente ao ajuste por covariáveis.
#
# IMPORTANTE
# ----------
# Este script NÃO cria um cutoff pós-hoc para declarar automaticamente um
# "vencedor". Ele mede os critérios predefinidos e preserva os resultados para
# interpretação conjunta.
#
# A comparação funcional principal é feita por GSEA reexecutada de maneira
# harmonizada nos datasets, usando:
#
#   - a mesma versão instalada de org.Hs.eg.db;
#   - o mesmo clusterProfiler;
#   - o mesmo ranking: Wald statistic do DESeq2;
#   - os mesmos minGSSize/maxGSSize;
#   - o mesmo critério de significância funcional: BH-FDR < 0,05.
#
# Isso evita comparar tabelas de enriquecimento geradas em etapas diferentes
# com parâmetros potencialmente distintos.
#
# Critério principal de DEG preservado:
#
#   BH-FDR < 0,05 e |log2FC| > 0,58
#
# Critério estrito:
#
#   BH-FDR < 0,05 e |log2FC| > 1
#
# ------------------------------------------------------------------------------


# ==============================================================================
# 0. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table",
  "ggplot2",
  "clusterProfiler",
  "org.Hs.eg.db"
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
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})


# ==============================================================================
# 1. PARÂMETROS PREDEFINIDOS
# ==============================================================================

PRIMARY_ALPHA <- 0.05
PRIMARY_LFC <- 0.58

STRICT_ALPHA <- 0.05
STRICT_LFC <- 1.00

FUNCTIONAL_ALPHA <- 0.05

GSEA_MIN_GS_SIZE <- 10L
GSEA_MAX_GS_SIZE <- 500L

TOP_PLOT_TERMS <- 15L


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

gse_dir <- file.path(
  "results",
  "candidate_GSE125583"
)

rosmap_dir <- file.path(
  "results",
  "candidate_ROSMAP"
)

comparison_dir <- file.path(
  rosmap_dir,
  "comparison_GSE125583"
)

gene_dir <- file.path(
  comparison_dir,
  "gene_level"
)

functional_dir <- file.path(
  comparison_dir,
  "functional"
)

plot_dir <- file.path(
  comparison_dir,
  "plots"
)

for (d in c(
  comparison_dir,
  gene_dir,
  functional_dir,
  plot_dir
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------------------------
# Resultados principais
# ------------------------------------------------------------------------------

gse_main_input <- file.path(
  gse_dir,
  "deseq2",
  "GSE125583_deseq2_results.tsv"
)

rosmap_main_input <- file.path(
  rosmap_dir,
  "deseq2",
  "ROSMAP_deseq2_results.tsv"
)

rosmap_adjusted_input <- file.path(
  rosmap_dir,
  "deseq2",
  "ROSMAP_deseq2_results_adjusted.tsv"
)


# ------------------------------------------------------------------------------
# O arquivo ajustado do GSE125583 pode ter sido salvo com nomes diferentes
# durante a Fase 01b. Primeiro tentamos nomes explícitos; se não existir, o
# script procura de forma controlada dentro de candidate_GSE125583/deseq2.
# ------------------------------------------------------------------------------

gse_adjusted_candidates <- c(
  file.path(
    gse_dir,
    "deseq2",
    "GSE125583_deseq2_results_adjusted.tsv"
  ),
  file.path(
    gse_dir,
    "deseq2",
    "GSE125583_deseq2_adjusted_results.tsv"
  )
)


# ------------------------------------------------------------------------------
# Critérios predefinidos da auditoria do ROSMAP, se disponíveis.
# ------------------------------------------------------------------------------

predefined_criteria_input <- file.path(
  rosmap_dir,
  "audit",
  "00_predefined_selection_criteria.tsv"
)


# ==============================================================================
# 3. SAÍDAS
# ==============================================================================

gene_summary_output <- file.path(
  gene_dir,
  "AD_candidate_gene_concordance_summary.tsv"
)

gene_common_output <- file.path(
  gene_dir,
  "GSE125583_ROSMAP_common_gene_effects.tsv"
)

deg_overlap_output <- file.path(
  gene_dir,
  "GSE125583_ROSMAP_primary_DEG_overlap.tsv"
)

adjusted_gene_common_output <- file.path(
  gene_dir,
  "GSE125583_ROSMAP_adjusted_common_gene_effects.tsv"
)

within_dataset_stability_output <- file.path(
  gene_dir,
  "AD_candidates_within_dataset_stability.tsv"
)

functional_summary_output <- file.path(
  functional_dir,
  "AD_candidate_GSEA_concordance_summary.tsv"
)

functional_terms_output <- file.path(
  functional_dir,
  "AD_candidate_GSEA_term_level.tsv"
)

gsea_status_output <- file.path(
  functional_dir,
  "AD_candidate_GSEA_status.tsv"
)

criteria_summary_output <- file.path(
  comparison_dir,
  "AD_candidate_predefined_criteria_summary.tsv"
)

report_output <- file.path(
  comparison_dir,
  "07_compare_ROSMAP_GSE125583.txt"
)

session_info_output <- file.path(
  comparison_dir,
  "07_compare_ROSMAP_GSE125583_sessionInfo.txt"
)


# ==============================================================================
# 4. FUNÇÕES AUXILIARES
# ==============================================================================

assert_file <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Arquivo não encontrado:\n",
      path
    )
  }
  invisible(TRUE)
}


assert_columns <- function(
  dt,
  required,
  object_name
) {
  missing <- setdiff(
    required,
    names(dt)
  )

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


fmt_num <- function(
  x,
  digits = 4
) {
  format(
    round(as.numeric(x), digits),
    nsmall = digits,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


safe_cor <- function(
  x,
  y,
  method = "spearman"
) {
  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 3) {
    return(NA_real_)
  }

  suppressWarnings(
    cor(
      x[ok],
      y[ok],
      method = method
    )
  )
}


safe_jaccard <- function(
  x,
  y
) {
  x <- unique(
    x[
      !is.na(x)
    ]
  )

  y <- unique(
    y[
      !is.na(y)
    ]
  )

  u <- union(
    x,
    y
  )

  if (length(u) == 0) {
    return(NA_real_)
  }

  length(
    intersect(
      x,
      y
    )
  ) /
    length(u)
}


classify_results_if_needed <- function(dt) {

  dt <- copy(dt)

  dt[
    ,
    GeneID := as.character(GeneID)
  ]

  if (!"DEG" %in% names(dt)) {
    dt[
      ,
      DEG :=
        !is.na(padj) &
        padj < PRIMARY_ALPHA &
        abs(log2FoldChange) > PRIMARY_LFC
    ]
  }

  if (!"direction" %in% names(dt)) {
    dt[
      ,
      direction := fifelse(
        DEG &
          log2FoldChange > PRIMARY_LFC,
        "Up",
        fifelse(
          DEG &
            log2FoldChange < -PRIMARY_LFC,
          "Down",
          "Not_DEG"
        )
      )
    ]
  }

  if (!"DEG_strict" %in% names(dt)) {
    if ("DEG_sensitivity" %in% names(dt)) {
      dt[
        ,
        DEG_strict := as.logical(DEG_sensitivity)
      ]
    } else {
      dt[
        ,
        DEG_strict :=
          !is.na(padj) &
          padj < STRICT_ALPHA &
          abs(log2FoldChange) > STRICT_LFC
      ]
    }
  }

  dt
}


# ------------------------------------------------------------------------------
# Localizar opcionalmente um resultado ajustado dentro de uma pasta.
# Só aceita arquivos que contenham as colunas necessárias.
# ------------------------------------------------------------------------------

find_adjusted_result <- function(
  explicit_candidates,
  search_dir,
  dataset_pattern
) {

  direct <- explicit_candidates[
    file.exists(explicit_candidates)
  ]

  if (length(direct) > 0) {
    return(
      normalizePath(
        direct[1],
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  if (!dir.exists(search_dir)) {
    return(NA_character_)
  }

  candidates <- list.files(
    search_dir,
    pattern = "\\.(tsv|txt|csv)(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  candidates <- candidates[
    grepl(
      dataset_pattern,
      basename(candidates),
      ignore.case = TRUE
    ) &
      grepl(
        "adjust",
        basename(candidates),
        ignore.case = TRUE
      ) &
      grepl(
        "deseq|result",
        basename(candidates),
        ignore.case = TRUE
      )
  ]

  if (length(candidates) == 0) {
    return(NA_character_)
  }

  valid <- character()

  for (f in candidates) {

    hdr <- tryCatch(
      names(
        fread(
          f,
          nrows = 3,
          showProgress = FALSE
        )
      ),
      error = function(e) character()
    )

    needed <- c(
      "GeneID",
      "log2FoldChange",
      "stat",
      "pvalue",
      "padj"
    )

    if (all(needed %in% hdr)) {
      valid <- c(valid, f)
    }
  }

  if (length(valid) == 0) {
    return(NA_character_)
  }

  # Preferir nome mais curto e mais diretamente relacionado ao resultado.
  valid <- valid[
    order(
      nchar(valid),
      valid
    )
  ]

  normalizePath(
    valid[1],
    winslash = "/",
    mustWork = TRUE
  )
}


# ------------------------------------------------------------------------------
# Resumo de estabilidade dentro de um mesmo dataset.
# ------------------------------------------------------------------------------

summarize_within_dataset <- function(
  main_dt,
  adjusted_dt,
  dataset_name
) {

  common <- merge(
    main_dt[
      ,
      .(
        GeneID,
        log2FC_main = log2FoldChange,
        DEG_main = DEG,
        direction_main = direction
      )
    ],
    adjusted_dt[
      ,
      .(
        GeneID,
        log2FC_adjusted = log2FoldChange,
        DEG_adjusted = DEG,
        direction_adjusted = direction
      )
    ],
    by = "GeneID",
    all = FALSE,
    sort = FALSE
  )

  shared <- common[
    DEG_main &
      DEG_adjusted
  ]

  union_n <- common[
    DEG_main |
      DEG_adjusted,
    .N
  ]

  data.table(
    dataset = dataset_name,
    genes_compared = nrow(common),
    spearman_log2FC = safe_cor(
      common$log2FC_main,
      common$log2FC_adjusted,
      "spearman"
    ),
    pearson_log2FC = safe_cor(
      common$log2FC_main,
      common$log2FC_adjusted,
      "pearson"
    ),
    primary_main = common[
      DEG_main == TRUE,
      .N
    ],
    primary_adjusted = common[
      DEG_adjusted == TRUE,
      .N
    ],
    shared_primary = nrow(shared),
    shared_same_direction = shared[
      direction_main ==
        direction_adjusted,
      .N
    ],
    shared_opposite_direction = shared[
      direction_main !=
        direction_adjusted,
      .N
    ],
    primary_union = union_n,
    primary_Jaccard = if (
      union_n == 0
    ) {
      NA_real_
    } else {
      nrow(shared) / union_n
    }
  )
}


# ==============================================================================
# 5. LER RESULTADOS DESeq2
# ==============================================================================

assert_file(gse_main_input)
assert_file(rosmap_main_input)
assert_file(rosmap_adjusted_input)

required_result_cols <- c(
  "GeneID",
  "log2FoldChange",
  "stat",
  "pvalue",
  "padj"
)

gse_main <- fread(
  gse_main_input,
  na.strings = c("NA", "")
)

rosmap_main <- fread(
  rosmap_main_input,
  na.strings = c("NA", "")
)

rosmap_adjusted <- fread(
  rosmap_adjusted_input,
  na.strings = c("NA", "")
)

assert_columns(
  gse_main,
  required_result_cols,
  "GSE125583 principal"
)

assert_columns(
  rosmap_main,
  required_result_cols,
  "ROSMAP principal"
)

assert_columns(
  rosmap_adjusted,
  required_result_cols,
  "ROSMAP ajustado"
)

gse_main <- classify_results_if_needed(
  gse_main
)

rosmap_main <- classify_results_if_needed(
  rosmap_main
)

rosmap_adjusted <- classify_results_if_needed(
  rosmap_adjusted
)

if (
  anyDuplicated(gse_main$GeneID) ||
  anyDuplicated(rosmap_main$GeneID) ||
  anyDuplicated(rosmap_adjusted$GeneID)
) {
  stop(
    "GeneID duplicado em pelo menos um resultado DESeq2."
  )
}


# ------------------------------------------------------------------------------
# Resultado ajustado do GSE125583 — opcional, mas desejável para o critério C.
# ------------------------------------------------------------------------------

gse_adjusted_input <- find_adjusted_result(
  explicit_candidates = gse_adjusted_candidates,
  search_dir = file.path(
    gse_dir,
    "deseq2"
  ),
  dataset_pattern = "GSE125583"
)

gse_adjusted_available <- !is.na(
  gse_adjusted_input
)

if (gse_adjusted_available) {

  gse_adjusted <- fread(
    gse_adjusted_input,
    na.strings = c("NA", "")
  )

  assert_columns(
    gse_adjusted,
    required_result_cols,
    "GSE125583 ajustado"
  )

  gse_adjusted <- classify_results_if_needed(
    gse_adjusted
  )

  if (anyDuplicated(gse_adjusted$GeneID)) {
    stop(
      "GeneID duplicado no resultado ajustado do GSE125583."
    )
  }

} else {

  gse_adjusted <- NULL
}


# ==============================================================================
# 6. UNIVERSOS TESTÁVEIS
# ==============================================================================

gse_universe <- unique(
  gse_main$GeneID
)

rosmap_universe <- unique(
  rosmap_main$GeneID
)

common_universe <- intersect(
  gse_universe,
  rosmap_universe
)

universe_union <- union(
  gse_universe,
  rosmap_universe
)

universe_summary <- data.table(
  metric = c(
    "GSE125583_tested_genes",
    "ROSMAP_tested_genes",
    "common_tested_genes",
    "union_tested_genes",
    "testable_universe_Jaccard",
    "GSE125583_covered_by_common",
    "ROSMAP_covered_by_common"
  ),
  value = c(
    length(gse_universe),
    length(rosmap_universe),
    length(common_universe),
    length(universe_union),
    length(common_universe) /
      length(universe_union),
    length(common_universe) /
      length(gse_universe),
    length(common_universe) /
      length(rosmap_universe)
  )
)

fwrite(
  universe_summary,
  file.path(
    gene_dir,
    "AD_candidate_testable_universe_summary.tsv"
  ),
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 7. CONCORDÂNCIA GENE-A-GENE — MODELO PRINCIPAL
# ==============================================================================

common_main <- merge(
  gse_main[
    GeneID %in% common_universe,
    .(
      GeneID,
      symbol_GSE125583 = if (
        "symbol" %in% names(gse_main)
      ) {
        symbol
      } else {
        NA_character_
      },
      log2FC_GSE125583 = log2FoldChange,
      stat_GSE125583 = stat,
      pvalue_GSE125583 = pvalue,
      padj_GSE125583 = padj,
      DEG_GSE125583 = DEG,
      direction_GSE125583 = direction,
      DEG_strict_GSE125583 = DEG_strict
    )
  ],
  rosmap_main[
    GeneID %in% common_universe,
    .(
      GeneID,
      symbol_ROSMAP = if (
        "symbol" %in% names(rosmap_main)
      ) {
        symbol
      } else {
        NA_character_
      },
      log2FC_ROSMAP = log2FoldChange,
      stat_ROSMAP = stat,
      pvalue_ROSMAP = pvalue,
      padj_ROSMAP = padj,
      DEG_ROSMAP = DEG,
      direction_ROSMAP = direction,
      DEG_strict_ROSMAP = DEG_strict
    )
  ],
  by = "GeneID",
  all = FALSE,
  sort = FALSE
)

common_main[
  ,
  symbol := fifelse(
    !is.na(symbol_GSE125583) &
      symbol_GSE125583 != "",
    symbol_GSE125583,
    fifelse(
      !is.na(symbol_ROSMAP) &
        symbol_ROSMAP != "",
      symbol_ROSMAP,
      GeneID
    )
  )
]

common_main[
  ,
  primary_overlap_class := fifelse(
    DEG_GSE125583 &
      DEG_ROSMAP,
    "Shared_primary",
    fifelse(
      DEG_GSE125583 &
        !DEG_ROSMAP,
      "GSE125583_only",
      fifelse(
        !DEG_GSE125583 &
          DEG_ROSMAP,
        "ROSMAP_only",
        "Neither"
      )
    )
  )
]

common_main[
  ,
  same_effect_sign :=
    sign(log2FC_GSE125583) ==
    sign(log2FC_ROSMAP)
]

common_main[
  ,
  informative_effect :=
    abs(log2FC_GSE125583) >
      PRIMARY_LFC |
    abs(log2FC_ROSMAP) >
      PRIMARY_LFC
]

common_main[
  ,
  same_sign_informative :=
    informative_effect &
    same_effect_sign
]

common_main[
  ,
  opposite_sign_informative :=
    informative_effect &
    !same_effect_sign
]

fwrite(
  common_main,
  gene_common_output,
  sep = "\t",
  na = "NA"
)


# ------------------------------------------------------------------------------
# Métricas globais
# ------------------------------------------------------------------------------

rho_main <- safe_cor(
  common_main$log2FC_GSE125583,
  common_main$log2FC_ROSMAP,
  "spearman"
)

pearson_main <- safe_cor(
  common_main$log2FC_GSE125583,
  common_main$log2FC_ROSMAP,
  "pearson"
)

informative_dt <- common_main[
  informative_effect == TRUE
]

rho_informative <- safe_cor(
  informative_dt$log2FC_GSE125583,
  informative_dt$log2FC_ROSMAP,
  "spearman"
)

pearson_informative <- safe_cor(
  informative_dt$log2FC_GSE125583,
  informative_dt$log2FC_ROSMAP,
  "pearson"
)


# ------------------------------------------------------------------------------
# Sobreposição dos DEGs no MESMO universo comum.
# Isso evita atribuir ao método diferenças causadas apenas por genes que não
# foram testados em uma das coortes.
# ------------------------------------------------------------------------------

gse_primary_common <- common_main[
  DEG_GSE125583 == TRUE,
  GeneID
]

rosmap_primary_common <- common_main[
  DEG_ROSMAP == TRUE,
  GeneID
]

shared_primary <- intersect(
  gse_primary_common,
  rosmap_primary_common
)

primary_union <- union(
  gse_primary_common,
  rosmap_primary_common
)

shared_dt <- common_main[
  GeneID %in% shared_primary
]

n_shared_same_direction <- shared_dt[
  direction_GSE125583 ==
    direction_ROSMAP,
  .N
]

n_shared_opposite_direction <- shared_dt[
  direction_GSE125583 !=
    direction_ROSMAP,
  .N
]

primary_jaccard <- safe_jaccard(
  gse_primary_common,
  rosmap_primary_common
)

deg_overlap_table <- common_main[
  DEG_GSE125583 |
    DEG_ROSMAP,
  .(
    GeneID,
    symbol,
    log2FC_GSE125583,
    padj_GSE125583,
    direction_GSE125583,
    log2FC_ROSMAP,
    padj_ROSMAP,
    direction_ROSMAP,
    primary_overlap_class
  )
][
  order(
    primary_overlap_class,
    GeneID
  )
]

fwrite(
  deg_overlap_table,
  deg_overlap_output,
  sep = "\t",
  na = "NA"
)


# ------------------------------------------------------------------------------
# Correlação na união dos genes DEG em pelo menos uma coorte.
# ------------------------------------------------------------------------------

primary_union_dt <- common_main[
  GeneID %in% primary_union
]

rho_primary_union <- safe_cor(
  primary_union_dt$log2FC_GSE125583,
  primary_union_dt$log2FC_ROSMAP,
  "spearman"
)

pearson_primary_union <- safe_cor(
  primary_union_dt$log2FC_GSE125583,
  primary_union_dt$log2FC_ROSMAP,
  "pearson"
)


# ==============================================================================
# 8. COMPARAÇÃO GSE125583 PRINCIPAL × ROSMAP AJUSTADO
# ==============================================================================
#
# Esta comparação não substitui a análise principal principal×principal.
# Ela testa se a concordância externa do ROSMAP melhora ou piora quando o forte
# efeito de sequencingBatch e as demais covariáveis são explicitamente modelados.
# ==============================================================================

common_rosmap_adjusted <- merge(
  gse_main[
    GeneID %in% common_universe,
    .(
      GeneID,
      log2FC_GSE125583 = log2FoldChange,
      DEG_GSE125583 = DEG,
      direction_GSE125583 = direction
    )
  ],
  rosmap_adjusted[
    GeneID %in% common_universe,
    .(
      GeneID,
      log2FC_ROSMAP_adjusted = log2FoldChange,
      DEG_ROSMAP_adjusted = DEG,
      direction_ROSMAP_adjusted = direction
    )
  ],
  by = "GeneID",
  all = FALSE,
  sort = FALSE
)

fwrite(
  common_rosmap_adjusted,
  adjusted_gene_common_output,
  sep = "\t",
  na = "NA"
)

rho_gse_vs_rosmap_adjusted <- safe_cor(
  common_rosmap_adjusted$log2FC_GSE125583,
  common_rosmap_adjusted$log2FC_ROSMAP_adjusted,
  "spearman"
)

pearson_gse_vs_rosmap_adjusted <- safe_cor(
  common_rosmap_adjusted$log2FC_GSE125583,
  common_rosmap_adjusted$log2FC_ROSMAP_adjusted,
  "pearson"
)

gse_primary_for_adjusted <- common_rosmap_adjusted[
  DEG_GSE125583 == TRUE,
  GeneID
]

rosmap_adjusted_primary_common <- common_rosmap_adjusted[
  DEG_ROSMAP_adjusted == TRUE,
  GeneID
]

shared_primary_adjusted <- intersect(
  gse_primary_for_adjusted,
  rosmap_adjusted_primary_common
)

jaccard_primary_adjusted <- safe_jaccard(
  gse_primary_for_adjusted,
  rosmap_adjusted_primary_common
)

shared_adjusted_dt <- common_rosmap_adjusted[
  GeneID %in% shared_primary_adjusted
]

shared_adjusted_same_direction <- shared_adjusted_dt[
  direction_GSE125583 ==
    direction_ROSMAP_adjusted,
  .N
]

shared_adjusted_opposite_direction <- shared_adjusted_dt[
  direction_GSE125583 !=
    direction_ROSMAP_adjusted,
  .N
]


# ==============================================================================
# 9. RESUMO DA CONCORDÂNCIA GENE-A-GENE
# ==============================================================================

gene_summary <- data.table(
  comparison = c(
    "GSE125583_main_vs_ROSMAP_main",
    "GSE125583_main_vs_ROSMAP_adjusted"
  ),
  common_tested_genes = c(
    nrow(common_main),
    nrow(common_rosmap_adjusted)
  ),
  spearman_log2FC_all = c(
    rho_main,
    rho_gse_vs_rosmap_adjusted
  ),
  pearson_log2FC_all = c(
    pearson_main,
    pearson_gse_vs_rosmap_adjusted
  ),
  informative_genes_absLFC_gt_0.58_either = c(
    nrow(informative_dt),
    NA_integer_
  ),
  spearman_log2FC_informative = c(
    rho_informative,
    NA_real_
  ),
  pearson_log2FC_informative = c(
    pearson_informative,
    NA_real_
  ),
  primary_GSE125583_common = c(
    length(gse_primary_common),
    length(gse_primary_for_adjusted)
  ),
  primary_ROSMAP_common = c(
    length(rosmap_primary_common),
    length(rosmap_adjusted_primary_common)
  ),
  shared_primary = c(
    length(shared_primary),
    length(shared_primary_adjusted)
  ),
  primary_Jaccard = c(
    primary_jaccard,
    jaccard_primary_adjusted
  ),
  shared_same_direction = c(
    n_shared_same_direction,
    shared_adjusted_same_direction
  ),
  shared_opposite_direction = c(
    n_shared_opposite_direction,
    shared_adjusted_opposite_direction
  ),
  spearman_primary_union = c(
    rho_primary_union,
    NA_real_
  ),
  pearson_primary_union = c(
    pearson_primary_union,
    NA_real_
  )
)

fwrite(
  gene_summary,
  gene_summary_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 10. ESTABILIDADE DENTRO DE CADA COORTE
# ==============================================================================

rosmap_stability <- summarize_within_dataset(
  rosmap_main,
  rosmap_adjusted,
  "ROSMAP"
)

if (gse_adjusted_available) {

  gse_stability <- summarize_within_dataset(
    gse_main,
    gse_adjusted,
    "GSE125583"
  )

} else {

  gse_stability <- data.table(
    dataset = "GSE125583",
    genes_compared = NA_integer_,
    spearman_log2FC = NA_real_,
    pearson_log2FC = NA_real_,
    primary_main = sum(
      gse_main$DEG,
      na.rm = TRUE
    ),
    primary_adjusted = NA_integer_,
    shared_primary = NA_integer_,
    shared_same_direction = NA_integer_,
    shared_opposite_direction = NA_integer_,
    primary_union = NA_integer_,
    primary_Jaccard = NA_real_
  )
}

within_dataset_stability <- rbind(
  gse_stability,
  rosmap_stability,
  fill = TRUE
)

fwrite(
  within_dataset_stability,
  within_dataset_stability_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 11. FIGURAS DE CONCORDÂNCIA GENE-A-GENE
# ==============================================================================

gene_scatter_main_output <- file.path(
  plot_dir,
  "log2FC_GSE125583_vs_ROSMAP_main.png"
)

gene_scatter_adjusted_output <- file.path(
  plot_dir,
  "log2FC_GSE125583_vs_ROSMAP_adjusted.png"
)


plot_main_dt <- copy(
  common_main
)

plot_main_dt[
  ,
  plot_class := factor(
    primary_overlap_class,
    levels = c(
      "Neither",
      "GSE125583_only",
      "ROSMAP_only",
      "Shared_primary"
    )
  )
]

p_main <- ggplot(
  plot_main_dt,
  aes(
    x = log2FC_GSE125583,
    y = log2FC_ROSMAP,
    color = plot_class
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
    alpha = 0.45,
    size = 1.4
  ) +
  labs(
    title = "Alzheimer — concordância dos efeitos",
    subtitle = paste0(
      "GSE125583 principal vs ROSMAP principal | Spearman rho = ",
      fmt_num(rho_main, 3)
    ),
    x = "log2FC — GSE125583",
    y = "log2FC — ROSMAP",
    color = "DEG principal"
  ) +
  theme_bw(
    base_size = 12
  )

ggsave(
  gene_scatter_main_output,
  p_main,
  width = 8,
  height = 7,
  dpi = 300
)


plot_adjusted_dt <- copy(
  common_rosmap_adjusted
)

plot_adjusted_dt[
  ,
  plot_class := factor(
    fifelse(
      DEG_GSE125583 &
        DEG_ROSMAP_adjusted,
      "Shared_primary",
      fifelse(
        DEG_GSE125583 &
          !DEG_ROSMAP_adjusted,
        "GSE125583_only",
        fifelse(
          !DEG_GSE125583 &
            DEG_ROSMAP_adjusted,
          "ROSMAP_adjusted_only",
          "Neither"
        )
      )
    ),
    levels = c(
      "Neither",
      "GSE125583_only",
      "ROSMAP_adjusted_only",
      "Shared_primary"
    )
  )
]

p_adjusted <- ggplot(
  plot_adjusted_dt,
  aes(
    x = log2FC_GSE125583,
    y = log2FC_ROSMAP_adjusted,
    color = plot_class
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
    alpha = 0.45,
    size = 1.4
  ) +
  labs(
    title = "Alzheimer — concordância dos efeitos",
    subtitle = paste0(
      "GSE125583 principal vs ROSMAP ajustado | Spearman rho = ",
      fmt_num(
        rho_gse_vs_rosmap_adjusted,
        3
      )
    ),
    x = "log2FC — GSE125583",
    y = "log2FC — ROSMAP ajustado",
    color = "DEG principal"
  ) +
  theme_bw(
    base_size = 12
  )

ggsave(
  gene_scatter_adjusted_output,
  p_adjusted,
  width = 8,
  height = 7,
  dpi = 300
)


# ==============================================================================
# 12. GSEA HARMONIZADA
# ==============================================================================
#
# Reexecutamos GSEA nesta etapa para garantir comparação "apples-to-apples".
# ------------------------------------------------------------------------------

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
    stop(
      "GeneID duplicado no ranking GSEA."
    )
  }

  ranked <- ranked[
    order(-stat)
  ]

  gene_list <- ranked$stat
  names(gene_list) <- ranked$GeneID

  gene_list
}


cp_to_dt <- function(x) {

  if (is.null(x)) {
    return(data.table())
  }

  z <- tryCatch(
    as.data.frame(x),
    error = function(e) NULL
  )

  if (
    is.null(z) ||
    nrow(z) == 0
  ) {
    return(data.table())
  }

  as.data.table(z)
}


run_gsea <- function(
  res_dt,
  dataset_model,
  database
) {

  gene_list <- make_ranked_list(
    res_dt
  )

  if (database == "GO_BP") {

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

  } else if (database == "KEGG") {

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

  } else {

    stop(
      "Database GSEA não reconhecida: ",
      database
    )
  }

  if (inherits(obj, "error")) {

    return(
      list(
        result = data.table(),
        status = data.table(
          dataset_model = dataset_model,
          database = database,
          n_ranked_genes = length(gene_list),
          status = "error",
          message = conditionMessage(obj)
        )
      )
    )
  }

  dt <- cp_to_dt(obj)

  if (nrow(dt) > 0) {

    dt[
      ,
      dataset_model := dataset_model
    ]

    dt[
      ,
      database := database
    ]
  }

  list(
    result = dt,
    status = data.table(
      dataset_model = dataset_model,
      database = database,
      n_ranked_genes = length(gene_list),
      status = "completed",
      message = NA_character_
    )
  )
}


gsea_inputs <- list(
  GSE125583_main = gse_main,
  ROSMAP_main = rosmap_main,
  ROSMAP_adjusted = rosmap_adjusted
)

if (gse_adjusted_available) {
  gsea_inputs[["GSE125583_adjusted"]] <- gse_adjusted
}


gsea_results <- list()
gsea_status <- list()

cat(
  "\nExecutando GSEA harmonizada para comparação entre coortes...\n"
)

for (dataset_model in names(gsea_inputs)) {

  for (database in c(
    "GO_BP",
    "KEGG"
  )) {

    cat(
      "  ",
      dataset_model,
      " / ",
      database,
      "\n",
      sep = ""
    )

    run <- run_gsea(
      gsea_inputs[[dataset_model]],
      dataset_model,
      database
    )

    key <- paste(
      dataset_model,
      database,
      sep = "__"
    )

    gsea_results[[key]] <- run$result
    gsea_status[[length(gsea_status) + 1L]] <- run$status

    output_file <- file.path(
      functional_dir,
      paste0(
        "GSEA_",
        dataset_model,
        "_",
        database,
        ".tsv"
      )
    )

    if (nrow(run$result) == 0) {
      fwrite(
        data.table(
          note = "Nenhum resultado disponível."
        ),
        output_file,
        sep = "\t"
      )
    } else {
      fwrite(
        run$result,
        output_file,
        sep = "\t",
        na = "NA"
      )
    }
  }
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
# 13. COMPARAÇÃO FUNCIONAL POR GSEA
# ==============================================================================

get_sig_ids <- function(dt) {

  if (
    nrow(dt) == 0 ||
    !all(
      c(
        "ID",
        "p.adjust"
      ) %in% names(dt)
    )
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


compare_gsea_pair <- function(
  x_dt,
  y_dt,
  x_name,
  y_name,
  database
) {

  if (
    nrow(x_dt) == 0 ||
    nrow(y_dt) == 0
  ) {

    return(
      list(
        summary = data.table(
          comparison = paste0(
            x_name,
            "_vs_",
            y_name
          ),
          database = database,
          terms_tested_both = NA_integer_,
          spearman_NES_all_common_terms = NA_real_,
          significant_x = length(
            get_sig_ids(x_dt)
          ),
          significant_y = length(
            get_sig_ids(y_dt)
          ),
          shared_significant = NA_integer_,
          union_significant = NA_integer_,
          Jaccard = NA_real_,
          shared_same_direction = NA_integer_,
          shared_opposite_direction = NA_integer_,
          spearman_NES_union_significant = NA_real_
        ),
        terms = data.table()
      )
    )
  }

  x_sub <- x_dt[
    ,
    .(
      ID = as.character(ID),
      Description_x = Description,
      NES_x = NES,
      padj_x = p.adjust
    )
  ]

  y_sub <- y_dt[
    ,
    .(
      ID = as.character(ID),
      Description_y = Description,
      NES_y = NES,
      padj_y = p.adjust
    )
  ]

  merged <- merge(
    x_sub,
    y_sub,
    by = "ID",
    all = FALSE,
    sort = FALSE
  )

  merged[
    ,
    Description := fifelse(
      !is.na(Description_x),
      Description_x,
      Description_y
    )
  ]

  merged[
    ,
    significant_x :=
      !is.na(padj_x) &
      padj_x < FUNCTIONAL_ALPHA
  ]

  merged[
    ,
    significant_y :=
      !is.na(padj_y) &
      padj_y < FUNCTIONAL_ALPHA
  ]

  merged[
    ,
    shared_significant :=
      significant_x &
      significant_y
  ]

  merged[
    ,
    same_NES_direction :=
      sign(NES_x) ==
      sign(NES_y)
  ]

  merged[
    ,
    significant_either :=
      significant_x |
      significant_y
  ]

  sig_x_ids <- get_sig_ids(
    x_dt
  )

  sig_y_ids <- get_sig_ids(
    y_dt
  )

  shared_ids <- intersect(
    sig_x_ids,
    sig_y_ids
  )

  union_ids <- union(
    sig_x_ids,
    sig_y_ids
  )

  shared_rows <- merged[
    ID %in% shared_ids
  ]

  union_rows <- merged[
    ID %in% union_ids
  ]

  summary <- data.table(
    comparison = paste0(
      x_name,
      "_vs_",
      y_name
    ),
    database = database,
    terms_tested_both = nrow(merged),
    spearman_NES_all_common_terms = safe_cor(
      merged$NES_x,
      merged$NES_y,
      "spearman"
    ),
    significant_x = length(
      sig_x_ids
    ),
    significant_y = length(
      sig_y_ids
    ),
    shared_significant = length(
      shared_ids
    ),
    union_significant = length(
      union_ids
    ),
    Jaccard = safe_jaccard(
      sig_x_ids,
      sig_y_ids
    ),
    shared_same_direction = shared_rows[
      same_NES_direction == TRUE,
      .N
    ],
    shared_opposite_direction = shared_rows[
      same_NES_direction == FALSE,
      .N
    ],
    spearman_NES_union_significant = safe_cor(
      union_rows$NES_x,
      union_rows$NES_y,
      "spearman"
    )
  )

  terms <- merged[
    significant_either == TRUE,
    .(
      comparison = paste0(
        x_name,
        "_vs_",
        y_name
      ),
      database = database,
      ID,
      Description,
      NES_x,
      padj_x,
      NES_y,
      padj_y,
      significant_x,
      significant_y,
      shared_significant,
      same_NES_direction
    )
  ]

  list(
    summary = summary,
    terms = terms
  )
}


functional_comparisons <- list(
  c(
    "GSE125583_main",
    "ROSMAP_main"
  ),
  c(
    "GSE125583_main",
    "ROSMAP_adjusted"
  )
)

if (gse_adjusted_available) {
  functional_comparisons <- c(
    functional_comparisons,
    list(
      c(
        "GSE125583_adjusted",
        "ROSMAP_adjusted"
      )
    )
  )
}


functional_summary_list <- list()
functional_terms_list <- list()

for (pair in functional_comparisons) {

  x_name <- pair[1]
  y_name <- pair[2]

  for (database in c(
    "GO_BP",
    "KEGG"
  )) {

    x_key <- paste(
      x_name,
      database,
      sep = "__"
    )

    y_key <- paste(
      y_name,
      database,
      sep = "__"
    )

    cmp <- compare_gsea_pair(
      gsea_results[[x_key]],
      gsea_results[[y_key]],
      x_name,
      y_name,
      database
    )

    functional_summary_list[[length(functional_summary_list) + 1L]] <- cmp$summary

    if (nrow(cmp$terms) > 0) {
      functional_terms_list[[length(functional_terms_list) + 1L]] <- cmp$terms
    }
  }
}

functional_summary <- rbindlist(
  functional_summary_list,
  fill = TRUE
)

functional_terms <- if (
  length(
    functional_terms_list
  ) > 0
) {
  rbindlist(
    functional_terms_list,
    fill = TRUE
  )
} else {
  data.table()
}

fwrite(
  functional_summary,
  functional_summary_output,
  sep = "\t",
  na = "NA"
)

if (nrow(functional_terms) > 0) {

  fwrite(
    functional_terms,
    functional_terms_output,
    sep = "\t",
    na = "NA"
  )

} else {

  fwrite(
    data.table(
      note = "Nenhum termo significativo disponível para comparação."
    ),
    functional_terms_output,
    sep = "\t"
  )
}


# ==============================================================================
# 14. FIGURAS DE CONCORDÂNCIA FUNCIONAL
# ==============================================================================

make_gsea_scatter <- function(
  terms_dt,
  comparison_name,
  database_name,
  output_file
) {

  dt <- terms_dt[
    comparison ==
      comparison_name &
      database ==
      database_name
  ]

  if (nrow(dt) == 0) {
    return(
      invisible(NULL)
    )
  }

  dt[
    ,
    significance_class := factor(
      fifelse(
        significant_x &
          significant_y,
        "Significativo em ambos",
        fifelse(
          significant_x,
          "Somente GSE125583",
          "Somente ROSMAP"
        )
      ),
      levels = c(
        "Somente GSE125583",
        "Somente ROSMAP",
        "Significativo em ambos"
      )
    )
  ]

  p <- ggplot(
    dt[
      is.finite(NES_x) &
        is.finite(NES_y)
    ],
    aes(
      x = NES_x,
      y = NES_y,
      color = significance_class
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
      alpha = 0.7,
      size = 2
    ) +
    labs(
      title = paste0(
        "Alzheimer — concordância funcional | ",
        database_name
      ),
      subtitle = gsub(
        "_",
        " ",
        comparison_name
      ),
      x = "NES — GSE125583",
      y = "NES — ROSMAP",
      color = "Significância"
    ) +
    theme_bw(
      base_size = 12
    )

  ggsave(
    output_file,
    p,
    width = 8,
    height = 7,
    dpi = 300
  )
}


for (comparison_name in unique(
  functional_summary$comparison
)) {

  for (database_name in c(
    "GO_BP",
    "KEGG"
  )) {

    output_file <- file.path(
      plot_dir,
      paste0(
        "GSEA_NES_",
        comparison_name,
        "_",
        database_name,
        ".png"
      )
    )

    make_gsea_scatter(
      functional_terms,
      comparison_name,
      database_name,
      output_file
    )
  }
}


# ==============================================================================
# 15. RESUMO DOS CRITÉRIOS PREDEFINIDOS
# ==============================================================================
#
# Não há classificação automática "pass/fail".
# O arquivo reúne métricas observadas para os três critérios definidos antes.
# ==============================================================================

criteria_summary <- rbindlist(
  list(

    data.table(
      criterion = "A_log2FC_concordance",
      comparison = "GSE125583_main_vs_ROSMAP_main",
      metric = c(
        "common_tested_genes",
        "spearman_log2FC_all",
        "pearson_log2FC_all",
        "primary_Jaccard",
        "shared_primary",
        "shared_same_direction",
        "shared_opposite_direction"
      ),
      value = as.character(
        c(
          nrow(common_main),
          rho_main,
          pearson_main,
          primary_jaccard,
          length(shared_primary),
          n_shared_same_direction,
          n_shared_opposite_direction
        )
      )
    ),

    data.table(
      criterion = "A_log2FC_concordance",
      comparison = "GSE125583_main_vs_ROSMAP_adjusted",
      metric = c(
        "common_tested_genes",
        "spearman_log2FC_all",
        "pearson_log2FC_all",
        "primary_Jaccard",
        "shared_primary",
        "shared_same_direction",
        "shared_opposite_direction"
      ),
      value = as.character(
        c(
          nrow(common_rosmap_adjusted),
          rho_gse_vs_rosmap_adjusted,
          pearson_gse_vs_rosmap_adjusted,
          jaccard_primary_adjusted,
          length(shared_primary_adjusted),
          shared_adjusted_same_direction,
          shared_adjusted_opposite_direction
        )
      )
    ),

    functional_summary[
      ,
      .(
        criterion = "B_functional_coherence_GSEA",
        comparison,
        metric = paste0(
          database,
          "__",
          c(
            rep(
              "Jaccard",
              .N
            )
          )
        ),
        value = as.character(
          Jaccard
        )
      )
    ],

    within_dataset_stability[
      ,
      .(
        criterion = "C_covariate_stability",
        comparison = dataset,
        metric = "spearman_log2FC_main_vs_adjusted",
        value = as.character(
          spearman_log2FC
        )
      )
    ],

    within_dataset_stability[
      ,
      .(
        criterion = "C_covariate_stability",
        comparison = dataset,
        metric = "primary_Jaccard_main_vs_adjusted",
        value = as.character(
          primary_Jaccard
        )
      )
    ]
  ),
  fill = TRUE
)

fwrite(
  criteria_summary,
  criteria_summary_output,
  sep = "\t",
  na = "NA"
)


# ------------------------------------------------------------------------------
# Preservar uma cópia dos critérios predefinidos originais ao lado da
# comparação, caso o arquivo da auditoria exista.
# ------------------------------------------------------------------------------

if (file.exists(predefined_criteria_input)) {

  predefined_criteria <- fread(
    predefined_criteria_input,
    na.strings = c(
      "NA",
      ""
    )
  )

  fwrite(
    predefined_criteria,
    file.path(
      comparison_dir,
      "00_predefined_selection_criteria_copy.tsv"
    ),
    sep = "\t",
    na = "NA"
  )
}


# ==============================================================================
# 16. RELATÓRIO TEXTUAL
# ==============================================================================

sink(report_output)

cat(
  "============================================================\n"
)
cat(
  "FASE 01c — 07_compare_AD_candidates_ROSMAP_GSE125583.R\n"
)
cat(
  "COMPARAÇÃO FORMAL ROSMAP × GSE125583\n"
)
cat(
  "============================================================\n\n"
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "1. PRINCÍPIO DA COMPARAÇÃO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "A comparação segue critérios definidos antes da inspeção dos resultados\n"
)
cat(
  "do ROSMAP:\n\n"
)
cat(
  "  A. concordância dos log2FC entre as duas coortes;\n"
)
cat(
  "  B. coerência do enriquecimento funcional;\n"
)
cat(
  "  C. estabilidade dos efeitos frente às covariáveis.\n\n"
)
cat(
  "Nenhum cutoff pós-hoc é aplicado para escolher automaticamente uma coorte.\n\n"
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "2. UNIVERSO TESTÁVEL\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "GSE125583: ",
  fmt_int(length(gse_universe)),
  " genes\n",
  sep = ""
)

cat(
  "ROSMAP: ",
  fmt_int(length(rosmap_universe)),
  " genes\n",
  sep = ""
)

cat(
  "Interseção: ",
  fmt_int(length(common_universe)),
  " genes\n",
  sep = ""
)

cat(
  "União: ",
  fmt_int(length(universe_union)),
  " genes\n",
  sep = ""
)

cat(
  "Jaccard do universo: ",
  fmt_num(
    length(common_universe) /
      length(universe_union),
    4
  ),
  "\n\n",
  sep = ""
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "3. CRITÉRIO A — CONCORDÂNCIA GENE-A-GENE\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "GSE125583 principal × ROSMAP principal\n"
)

cat(
  "  genes comparados: ",
  fmt_int(nrow(common_main)),
  "\n",
  sep = ""
)

cat(
  "  Spearman log2FC: ",
  fmt_num(rho_main, 4),
  "\n",
  sep = ""
)

cat(
  "  Pearson log2FC: ",
  fmt_num(pearson_main, 4),
  "\n",
  sep = ""
)

cat(
  "  genes informativos (|log2FC| > 0,58 em pelo menos uma coorte): ",
  fmt_int(nrow(informative_dt)),
  "\n",
  sep = ""
)

cat(
  "  Spearman nos genes informativos: ",
  fmt_num(rho_informative, 4),
  "\n",
  sep = ""
)

cat(
  "  DEGs principais GSE125583 no universo comum: ",
  fmt_int(length(gse_primary_common)),
  "\n",
  sep = ""
)

cat(
  "  DEGs principais ROSMAP no universo comum: ",
  fmt_int(length(rosmap_primary_common)),
  "\n",
  sep = ""
)

cat(
  "  compartilhados: ",
  fmt_int(length(shared_primary)),
  "\n",
  sep = ""
)

cat(
  "  Jaccard dos DEGs: ",
  fmt_num(primary_jaccard, 4),
  "\n",
  sep = ""
)

cat(
  "  compartilhados mesma direção: ",
  fmt_int(n_shared_same_direction),
  "\n",
  sep = ""
)

cat(
  "  compartilhados direção oposta: ",
  fmt_int(n_shared_opposite_direction),
  "\n\n",
  sep = ""
)


cat(
  "GSE125583 principal × ROSMAP ajustado\n"
)

cat(
  "  Spearman log2FC: ",
  fmt_num(
    rho_gse_vs_rosmap_adjusted,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "  Pearson log2FC: ",
  fmt_num(
    pearson_gse_vs_rosmap_adjusted,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "  compartilhados: ",
  fmt_int(
    length(
      shared_primary_adjusted
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  Jaccard dos DEGs: ",
  fmt_num(
    jaccard_primary_adjusted,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "  compartilhados mesma direção: ",
  fmt_int(
    shared_adjusted_same_direction
  ),
  "\n",
  sep = ""
)

cat(
  "  compartilhados direção oposta: ",
  fmt_int(
    shared_adjusted_opposite_direction
  ),
  "\n\n",
  sep = ""
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "4. CRITÉRIO B — COERÊNCIA FUNCIONAL POR GSEA\n"
)
cat(
  "------------------------------------------------------------\n"
)

for (i in seq_len(nrow(functional_summary))) {

  row <- functional_summary[i]

  cat(
    row$comparison,
    " | ",
    row$database,
    "\n",
    sep = ""
  )

  cat(
    "  termos comuns testados: ",
    ifelse(
      is.na(row$terms_tested_both),
      "NA",
      fmt_int(row$terms_tested_both)
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Spearman NES — termos comuns: ",
    ifelse(
      is.na(row$spearman_NES_all_common_terms),
      "NA",
      fmt_num(
        row$spearman_NES_all_common_terms,
        4
      )
    ),
    "\n",
    sep = ""
  )

  cat(
    "  significativos GSE125583: ",
    fmt_int(row$significant_x),
    "\n",
    sep = ""
  )

  cat(
    "  significativos ROSMAP: ",
    fmt_int(row$significant_y),
    "\n",
    sep = ""
  )

  cat(
    "  compartilhados significativos: ",
    ifelse(
      is.na(row$shared_significant),
      "NA",
      fmt_int(row$shared_significant)
    ),
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

  cat(
    "  compartilhados mesma direção NES: ",
    ifelse(
      is.na(row$shared_same_direction),
      "NA",
      fmt_int(row$shared_same_direction)
    ),
    "\n",
    sep = ""
  )

  cat(
    "  compartilhados direção oposta NES: ",
    ifelse(
      is.na(row$shared_opposite_direction),
      "NA",
      fmt_int(row$shared_opposite_direction)
    ),
    "\n\n",
    sep = ""
  )
}


cat(
  "------------------------------------------------------------\n"
)
cat(
  "5. CRITÉRIO C — ESTABILIDADE COM COVARIÁVEIS\n"
)
cat(
  "------------------------------------------------------------\n"
)

for (i in seq_len(nrow(within_dataset_stability))) {

  row <- within_dataset_stability[i]

  cat(
    row$dataset,
    "\n",
    sep = ""
  )

  cat(
    "  Spearman main × adjusted: ",
    ifelse(
      is.na(row$spearman_log2FC),
      "NA",
      fmt_num(
        row$spearman_log2FC,
        4
      )
    ),
    "\n",
    sep = ""
  )

  cat(
    "  primary main: ",
    ifelse(
      is.na(row$primary_main),
      "NA",
      fmt_int(row$primary_main)
    ),
    "\n",
    sep = ""
  )

  cat(
    "  primary adjusted: ",
    ifelse(
      is.na(row$primary_adjusted),
      "NA",
      fmt_int(row$primary_adjusted)
    ),
    "\n",
    sep = ""
  )

  cat(
    "  shared primary: ",
    ifelse(
      is.na(row$shared_primary),
      "NA",
      fmt_int(row$shared_primary)
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Jaccard: ",
    ifelse(
      is.na(row$primary_Jaccard),
      "NA",
      fmt_num(
        row$primary_Jaccard,
        4
      )
    ),
    "\n\n",
    sep = ""
  )
}


cat(
  "------------------------------------------------------------\n"
)
cat(
  "6. CONTEXTO ANATÔMICO — NÃO É UM TESTE ESTATÍSTICO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "ROSMAP utiliza DLPFC, região diretamente alinhada à referência de\n"
)
cat(
  "deconvolução Mathys 2019 e anatomicamente mais próxima do córtex\n"
)
cat(
  "pré-frontal BA9 utilizado nos datasets atuais de PD e HD.\n"
)
cat(
  "GSE125583 utiliza giro fusiforme e permanece como coorte independente\n"
)
cat(
  "de Alzheimer com maior tamanho amostral.\n\n"
)
cat(
  "A região cerebral é interpretada em conjunto com os três critérios acima;\n"
)
cat(
  "não recebe um escore numérico inventado após os resultados.\n\n"
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "7. ARQUIVOS PRINCIPAIS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Concordância gene-a-gene:\n  ",
  gene_summary_output,
  "\n\n",
  sep = ""
)

cat(
  "Estabilidade dentro das coortes:\n  ",
  within_dataset_stability_output,
  "\n\n",
  sep = ""
)

cat(
  "Concordância funcional GSEA:\n  ",
  functional_summary_output,
  "\n\n",
  sep = ""
)

cat(
  "Termos GSEA:\n  ",
  functional_terms_output,
  "\n\n",
  sep = ""
)

cat(
  "Resumo dos critérios predefinidos:\n  ",
  criteria_summary_output,
  "\n\n",
  sep = ""
)

cat(
  "Scatter principal:\n  ",
  gene_scatter_main_output,
  "\n",
  sep = ""
)

cat(
  "Scatter ROSMAP ajustado:\n  ",
  gene_scatter_adjusted_output,
  "\n\n",
  sep = ""
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "8. INTERPRETAÇÃO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "- A escolha entre ROSMAP e GSE125583 NÃO é feita pelo número absoluto de DEGs.\n"
)
cat(
  "- A concordância gene-a-gene deve ser interpretada conjuntamente com a\n"
)
cat(
  "  coerência funcional e a estabilidade frente às covariáveis.\n"
)
cat(
  "- A GSEA foi reexecutada de modo harmonizado para reduzir diferenças de\n"
)
cat(
  "  implementação entre as fases.\n"
)
cat(
  "- O universo comum é usado nas comparações de DEGs para evitar confundir\n"
)
cat(
  "  diferenças de cobertura gênica com diferenças biológicas.\n"
)
cat(
  "- A compatibilidade anatômica permanece um argumento independente e deve\n"
)
cat(
  "  ser considerada na decisão final.\n\n"
)

if (gse_adjusted_available) {

  cat(
    "Resultado ajustado do GSE125583 localizado em:\n  ",
    gse_adjusted_input,
    "\n\n",
    sep = ""
  )

} else {

  cat(
    "ATENÇÃO: resultado ajustado do GSE125583 não foi localizado automaticamente.\n"
  )
  cat(
    "A estabilidade do GSE125583 ficará como NA neste relatório; isso NÃO afeta\n"
  )
  cat(
    "a comparação principal GSE125583 × ROSMAP nem a GSEA harmonizada principal.\n\n"
  )
}

cat(
  "============================================================\n"
)
cat(
  "07_compare_AD_candidates_ROSMAP_GSE125583.R CONCLUÍDO\n"
)
cat(
  "============================================================\n"
)

sink()


# ==============================================================================
# 17. SESSION INFO
# ==============================================================================

sink(session_info_output)
print(sessionInfo())
sink()


# ==============================================================================
# 18. RESUMO NO TERMINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "COMPARAÇÃO ROSMAP × GSE125583 CONCLUÍDA\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Universo comum: ",
  fmt_int(length(common_universe)),
  " genes\n",
  sep = ""
)

cat(
  "Spearman log2FC — principal × principal: ",
  fmt_num(rho_main, 4),
  "\n",
  sep = ""
)

cat(
  "Jaccard DEGs — principal × principal: ",
  fmt_num(primary_jaccard, 4),
  "\n",
  sep = ""
)

cat(
  "Spearman log2FC — GSE principal × ROSMAP ajustado: ",
  fmt_num(
    rho_gse_vs_rosmap_adjusted,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "Jaccard DEGs — GSE principal × ROSMAP ajustado: ",
  fmt_num(
    jaccard_primary_adjusted,
    4
  ),
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
  "Depois desta execução, revisar especialmente:\n"
)
cat(
  "  1) AD_candidate_gene_concordance_summary.tsv\n"
)
cat(
  "  2) AD_candidate_GSEA_concordance_summary.tsv\n"
)
cat(
  "  3) AD_candidates_within_dataset_stability.tsv\n"
)
