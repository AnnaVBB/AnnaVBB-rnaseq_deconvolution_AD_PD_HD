# ==============================================================================
# 08_finalize_AD_HD_PD_integration_ROSMAP_validation.R
# ==============================================================================
#
# FASE 01c — INTEGRAÇÃO FINAL AD–HD–PD + VALIDAÇÃO ROSMAP
#
# Decisão metodológica consolidada após a comparação formal das duas coortes
# de Alzheimer:
#
#   AD principal  = GSE125583
#   AD validação  = ROSMAP DLPFC
#   HD principal  = GSE64810
#   PD principal  = GSE68719
#
# O objetivo deste script é:
#
#   1) refazer a integração transdiagnóstica principal AD–HD–PD usando
#      GSE125583 como coorte principal de Alzheimer;
#
#   2) garantir que primary/strict sejam comparados dentro do MESMO universo
#      testável comum entre AD, HD e PD;
#
#   3) preservar a análise "article-like nominal" em listas completas por
#      dataset, exclusivamente para comparação histórica com a Fase 01;
#
#   4) quantificar interseções e concordância de direção dos efeitos;
#
#   5) gerar tabelas dos genes compartilhados AD–HD, AD–PD, HD–PD e AD–HD–PD;
#
#   6) usar ROSMAP como validação independente da componente de Alzheimer das
#      assinaturas transdiagnósticas, SEM substituir GSE125583 na integração
#      principal;
#
#   7) documentar explicitamente quais genes/signaturas são reproduzidos no
#      ROSMAP em direção e/ou significância.
#
# ------------------------------------------------------------------------------
# CRITÉRIOS PREDEFINIDOS
# ------------------------------------------------------------------------------
#
# Principal:
#   BH-FDR < 0,05 e |log2FC| > 0,58
#
# Estrito:
#   BH-FDR < 0,05 e |log2FC| > 1
#
# Comparativo nominal ("article-like"):
#   p < 0,05 e |log2FC| > 1
#
# IMPORTANTE:
#   - primary e strict usam o universo triplo comum testável;
#   - article-like nominal usa as listas COMPLETAS de cada dataset, preservando
#     a lógica da comparação histórica com a Fase 01;
#   - o ROSMAP NÃO entra na definição da assinatura principal AD–HD–PD;
#     ele entra apenas como validação independente de AD.
#
# ==============================================================================


# ==============================================================================
# 0. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table",
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
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})


# ==============================================================================
# 1. PARÂMETROS
# ==============================================================================

ALPHA_FDR <- 0.05

LFC_PRIMARY <- 0.58
LFC_STRICT <- 1.00

P_NOMINAL <- 0.05
LFC_NOMINAL <- 1.00

# ------------------------------------------------------------------------------
# Paleta visual padronizada do projeto
#
# Mantém exatamente o padrão utilizado na integração anterior:
#   AD = roxo
#   HD = azul
#   PD = rosa
# ------------------------------------------------------------------------------
COLOR_AD <- "#7B2CBF"
COLOR_HD <- "#2783A5"
COLOR_PD <- "#ED4C78"

# Valor esperado da interseção tripla com base nas etapas anteriores.
# NÃO é usado para filtrar; serve apenas como sanity check.
EXPECTED_COMMON_UNIVERSE <- 20576L


# ==============================================================================
# 2. CAMINHOS DE ENTRADA
# ==============================================================================

# ------------------------------------------------------------------------------
# Alzheimer principal — GSE125583
# ------------------------------------------------------------------------------

ad_input <- file.path(
  "results",
  "candidate_GSE125583",
  "deseq2",
  "GSE125583_deseq2_results.tsv"
)


# ------------------------------------------------------------------------------
# Huntington e Parkinson — resultados principais da Fase 01
#
# Esses arquivos historicamente usam "gene_id" em vez de "GeneID".
# A função standardize_result() normaliza isso automaticamente.
# ------------------------------------------------------------------------------

hd_input <- file.path(
  "results",
  "deseq2",
  "HD_deseq2_results.tsv"
)

pd_input <- file.path(
  "results",
  "deseq2",
  "PD_deseq2_results.tsv"
)


# ------------------------------------------------------------------------------
# ROSMAP — validação independente de Alzheimer
# ------------------------------------------------------------------------------

rosmap_main_input <- file.path(
  "results",
  "candidate_ROSMAP",
  "deseq2",
  "ROSMAP_deseq2_results.tsv"
)

rosmap_adjusted_input <- file.path(
  "results",
  "candidate_ROSMAP",
  "deseq2",
  "ROSMAP_deseq2_results_adjusted.tsv"
)


# ==============================================================================
# 3. CAMINHOS DE SAÍDA
# ==============================================================================

output_dir <- file.path(
  "results",
  "candidate_ROSMAP",
  "final_integration_AD_HD_PD"
)

tables_dir <- file.path(
  output_dir,
  "tables"
)

validation_dir <- file.path(
  output_dir,
  "ROSMAP_validation"
)

plot_dir <- file.path(
  output_dir,
  "plots"
)

for (d in c(
  output_dir,
  tables_dir,
  validation_dir,
  plot_dir
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------------------------
# Saídas principais
# ------------------------------------------------------------------------------

universe_summary_output <- file.path(
  tables_dir,
  "final_testable_universe_summary.tsv"
)

integration_summary_output <- file.path(
  tables_dir,
  "final_AD_HD_PD_integration_summary.tsv"
)

pairwise_summary_output <- file.path(
  tables_dir,
  "final_pairwise_overlap_summary.tsv"
)

threeway_regions_output <- file.path(
  tables_dir,
  "final_threeway_regions.tsv"
)

all_common_effects_output <- file.path(
  tables_dir,
  "final_common_universe_effects_AD_HD_PD.tsv"
)

ad_hd_shared_output <- file.path(
  tables_dir,
  "final_primary_shared_AD_HD.tsv"
)

ad_pd_shared_output <- file.path(
  tables_dir,
  "final_primary_shared_AD_PD.tsv"
)

hd_pd_shared_output <- file.path(
  tables_dir,
  "final_primary_shared_HD_PD.tsv"
)

triple_shared_output <- file.path(
  tables_dir,
  "final_primary_shared_AD_HD_PD.tsv"
)

transdiagnostic_ad_output <- file.path(
  tables_dir,
  "final_AD_primary_shared_with_HD_or_PD.tsv"
)

validation_summary_output <- file.path(
  validation_dir,
  "ROSMAP_validation_summary.tsv"
)

validation_genes_output <- file.path(
  validation_dir,
  "ROSMAP_validation_transdiagnostic_AD_genes.tsv"
)

report_output <- file.path(
  output_dir,
  "08_finalize_AD_HD_PD_integration_ROSMAP_validation.txt"
)

session_info_output <- file.path(
  output_dir,
  "08_sessionInfo.txt"
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


fmt_int <- function(x) {

  if (length(x) == 0 || is.na(x)) {
    return("NA")
  }

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

  if (length(x) == 0 || is.na(x)) {
    return("NA")
  }

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


# ------------------------------------------------------------------------------
# Escolher coluna de símbolo, quando disponível.
# ------------------------------------------------------------------------------

extract_symbol <- function(dt) {

  candidates <- c(
    "symbol",
    "Symbol",
    "gene_symbol",
    "GeneSymbol",
    "plot_label"
  )

  hit <- candidates[
    candidates %in% names(dt)
  ]

  if (length(hit) == 0) {
    return(
      rep(
        NA_character_,
        nrow(dt)
      )
    )
  }

  as.character(
    dt[[hit[1]]]
  )
}


# ------------------------------------------------------------------------------
# Padronizar resultados DESeq2 de diferentes fases.
#
# Saída interna:
#
#   GeneID
#   symbol
#   baseMean
#   log2FoldChange
#   lfcSE
#   stat
#   pvalue
#   padj
#   DEG_primary
#   DEG_strict
#   DEG_nominal
#   direction
# ------------------------------------------------------------------------------

standardize_result <- function(
  dt,
  dataset_name
) {

  dt <- copy(dt)

  # GeneID pode estar como GeneID ou gene_id.
  if (!"GeneID" %in% names(dt)) {

    if ("gene_id" %in% names(dt)) {

      setnames(
        dt,
        "gene_id",
        "GeneID"
      )

    } else {

      stop(
        dataset_name,
        ": não foi encontrada coluna GeneID/gene_id."
      )
    }
  }

  required <- c(
    "GeneID",
    "log2FoldChange",
    "pvalue",
    "padj"
  )

  missing <- setdiff(
    required,
    names(dt)
  )

  if (length(missing) > 0) {
    stop(
      dataset_name,
      ": colunas obrigatórias ausentes: ",
      paste(missing, collapse = ", ")
    )
  }

  dt[
    ,
    GeneID := as.character(GeneID)
  ]

  if (anyNA(dt$GeneID) || any(trimws(dt$GeneID) == "")) {
    stop(
      dataset_name,
      ": GeneID ausente."
    )
  }

  if (anyDuplicated(dt$GeneID)) {
    stop(
      dataset_name,
      ": GeneID duplicado."
    )
  }

  # Colunas opcionais são criadas se não existirem.
  if (!"baseMean" %in% names(dt)) {
    dt[
      ,
      baseMean := NA_real_
    ]
  }

  if (!"lfcSE" %in% names(dt)) {
    dt[
      ,
      lfcSE := NA_real_
    ]
  }

  if (!"stat" %in% names(dt)) {
    dt[
      ,
      stat := NA_real_
    ]
  }

  dt[
    ,
    symbol := extract_symbol(dt)
  ]

  # Critério principal.
  dt[
    ,
    DEG_primary :=
      !is.na(padj) &
      padj < ALPHA_FDR &
      abs(log2FoldChange) > LFC_PRIMARY
  ]

  # Critério estrito.
  dt[
    ,
    DEG_strict :=
      !is.na(padj) &
      padj < ALPHA_FDR &
      abs(log2FoldChange) > LFC_STRICT
  ]

  # Critério nominal utilizado apenas na comparação histórica.
  dt[
    ,
    DEG_nominal :=
      !is.na(pvalue) &
      pvalue < P_NOMINAL &
      abs(log2FoldChange) > LFC_NOMINAL
  ]

  # Direção sempre recalculada a partir do sinal do log2FC.
  dt[
    ,
    direction := fifelse(
      log2FoldChange > 0,
      "Up",
      fifelse(
        log2FoldChange < 0,
        "Down",
        "Zero"
      )
    )
  ]

  dt[
    ,
    dataset := dataset_name
  ]

  # Retornar apenas colunas padronizadas.
  dt[
    ,
    .(
      GeneID,
      symbol,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      DEG_primary,
      DEG_strict,
      DEG_nominal,
      direction,
      dataset
    )
  ]
}


# ------------------------------------------------------------------------------
# Selecionar genes DE por critério.
# ------------------------------------------------------------------------------

get_gene_set <- function(
  dt,
  criterion,
  universe = NULL
) {

  if (!criterion %in% c(
    "primary",
    "strict",
    "nominal"
  )) {
    stop(
      "Critério desconhecido: ",
      criterion
    )
  }

  work <- dt

  if (!is.null(universe)) {
    work <- work[
      GeneID %in% universe
    ]
  }

  flag <- switch(
    criterion,
    primary = "DEG_primary",
    strict = "DEG_strict",
    nominal = "DEG_nominal"
  )

  unique(
    work[
      get(flag) == TRUE,
      GeneID
    ]
  )
}


# ------------------------------------------------------------------------------
# Resumo de sobreposição de dois conjuntos.
# ------------------------------------------------------------------------------

pairwise_set_summary <- function(
  x_set,
  y_set,
  x_name,
  y_name,
  criterion
) {

  shared <- intersect(
    x_set,
    y_set
  )

  data.table(
    criterion = criterion,
    comparison = paste0(
      x_name,
      "_vs_",
      y_name
    ),
    n_x = length(x_set),
    n_y = length(y_set),
    shared = length(shared),
    union = length(
      union(
        x_set,
        y_set
      )
    ),
    Jaccard = safe_jaccard(
      x_set,
      y_set
    )
  )
}


# ------------------------------------------------------------------------------
# Calcular as sete regiões exclusivas de um Venn de três conjuntos.
# ------------------------------------------------------------------------------

threeway_regions <- function(
  A,
  B,
  C,
  criterion
) {

  A_only <- setdiff(
    A,
    union(
      B,
      C
    )
  )

  B_only <- setdiff(
    B,
    union(
      A,
      C
    )
  )

  C_only <- setdiff(
    C,
    union(
      A,
      B
    )
  )

  AB_only <- setdiff(
    intersect(
      A,
      B
    ),
    C
  )

  AC_only <- setdiff(
    intersect(
      A,
      C
    ),
    B
  )

  BC_only <- setdiff(
    intersect(
      B,
      C
    ),
    A
  )

  ABC <- Reduce(
    intersect,
    list(
      A,
      B,
      C
    )
  )

  data.table(
    criterion = criterion,
    region = c(
      "AD_only",
      "HD_only",
      "PD_only",
      "AD_HD_only",
      "AD_PD_only",
      "HD_PD_only",
      "AD_HD_PD"
    ),
    N = c(
      length(A_only),
      length(B_only),
      length(C_only),
      length(AB_only),
      length(AC_only),
      length(BC_only),
      length(ABC)
    )
  )
}


# ------------------------------------------------------------------------------
# Venn simples em ggplot2.
#
# Os círculos NÃO são proporcionais às áreas. O objetivo é documentar
# visualmente as sete regiões e seus números sem adicionar dependência externa.
# ------------------------------------------------------------------------------

circle_df <- function(
  cx,
  cy,
  r,
  label
) {

  theta <- seq(
    0,
    2 * pi,
    length.out = 400
  )

  data.table(
    x = cx + r * cos(theta),
    y = cy + r * sin(theta),
    set = label
  )
}


plot_threeway_venn <- function(
  region_dt,
  title,
  subtitle,
  output_file
) {

  # --------------------------------------------------------------------------
  # Esta função reproduz o padrão visual já estabelecido na Fase 01b:
  #
  #   AD = roxo  (#7B2CBF)
  #   HD = azul  (#2783A5)
  #   PD = rosa  (#ED4C78)
  #
  # O fundo é explicitamente branco para evitar PNG transparente/escuro.
  # --------------------------------------------------------------------------

  circle_points <- function(
    cx,
    cy,
    radius = 2.35,
    n = 500,
    group
  ) {

    theta <- seq(
      0,
      2 * pi,
      length.out = n
    )

    data.table(
      x = cx + radius * cos(theta),
      y = cy + radius * sin(theta),
      group = group
    )
  }

  get_n <- function(region_name) {

    value <- region_dt[
      region == region_name,
      N
    ]

    if (length(value) == 0) {
      return(0L)
    }

    as.integer(value[1])
  }

  # Mesma geometria utilizada no script 08 da Fase 01b.
  circles <- rbindlist(
    list(
      circle_points(
        3.90,
        5.45,
        group = "AD"
      ),
      circle_points(
        6.60,
        5.45,
        group = "HD"
      ),
      circle_points(
        5.25,
        3.15,
        group = "PD"
      )
    )
  )

  # Posições das sete regiões do Venn.
  label_positions <- data.table(
    region = c(
      "AD_only",
      "HD_only",
      "PD_only",
      "AD_HD_only",
      "AD_PD_only",
      "HD_PD_only",
      "AD_HD_PD"
    ),
    x = c(
      3.25,
      7.25,
      5.25,
      5.25,
      4.10,
      6.40,
      5.25
    ),
    y = c(
      5.75,
      5.75,
      2.20,
      6.35,
      4.05,
      4.05,
      4.95
    ),
    n = c(
      get_n("AD_only"),
      get_n("HD_only"),
      get_n("PD_only"),
      get_n("AD_HD_only"),
      get_n("AD_PD_only"),
      get_n("HD_PD_only"),
      get_n("AD_HD_PD")
    )
  )

  # Rótulos das doenças.
  set_labels <- data.table(
    label = c(
      "AD",
      "HD",
      "PD"
    ),
    x = c(
      2.35,
      8.15,
      5.25
    ),
    y = c(
      8.15,
      8.15,
      0.10
    ),
    color = c(
      COLOR_AD,
      COLOR_HD,
      COLOR_PD
    )
  )

  p <- ggplot() +

    # AD — roxo.
    geom_polygon(
      data = circles[
        group == "AD"
      ],
      aes(
        x = x,
        y = y
      ),
      fill = COLOR_AD,
      alpha = 0.25,
      color = COLOR_AD,
      linewidth = 1.25
    ) +

    # HD — azul.
    geom_polygon(
      data = circles[
        group == "HD"
      ],
      aes(
        x = x,
        y = y
      ),
      fill = COLOR_HD,
      alpha = 0.25,
      color = COLOR_HD,
      linewidth = 1.25
    ) +

    # PD — rosa.
    geom_polygon(
      data = circles[
        group == "PD"
      ],
      aes(
        x = x,
        y = y
      ),
      fill = COLOR_PD,
      alpha = 0.25,
      color = COLOR_PD,
      linewidth = 1.25
    ) +

    # Números das sete regiões, sem caixas, como no padrão anterior.
    geom_text(
      data = label_positions,
      aes(
        x = x,
        y = y,
        label = n
      ),
      size = 7,
      color = "black"
    ) +

    # Nome de cada doença na cor correspondente.
    geom_text(
      data = set_labels,
      aes(
        x = x,
        y = y,
        label = label,
        color = color
      ),
      size = 10,
      fontface = "bold",
      show.legend = FALSE
    ) +

    scale_color_identity() +

    coord_fixed(
      xlim = c(
        0.7,
        9.8
      ),
      ylim = c(
        -0.4,
        9.5
      ),
      clip = "off"
    ) +

    labs(
      title = title,
      subtitle = subtitle
    ) +

    theme_void(
      base_size = 16
    ) +

    theme(
      plot.background = element_rect(
        fill = "white",
        color = NA
      ),
      panel.background = element_rect(
        fill = "white",
        color = NA
      ),
      plot.title = element_text(
        face = "bold",
        size = 22,
        hjust = 0.5,
        margin = margin(
          b = 5
        ),
        color = "black"
      ),
      plot.subtitle = element_text(
        size = 14,
        hjust = 0.5,
        margin = margin(
          b = 15
        ),
        color = "black"
      ),
      plot.margin = margin(
        20,
        20,
        20,
        20
      )
    )

  ggsave(
    filename = output_file,
    plot = p,
    width = 11,
    height = 8.5,
    dpi = 300,
    bg = "white"
  )

  invisible(p)
}


# ------------------------------------------------------------------------------
# Obter símbolo preferencial de várias colunas.
# ------------------------------------------------------------------------------

coalesce_symbol <- function(...) {

  values <- list(...)

  out <- rep(
    NA_character_,
    length(values[[1]])
  )

  for (v in values) {

    v <- as.character(v)

    replace_idx <- (
      is.na(out) |
        out == ""
    ) &
      !is.na(v) &
      v != ""

    out[
      replace_idx
    ] <- v[
      replace_idx
    ]
  }

  out
}


# ------------------------------------------------------------------------------
# Construir tabela compartilhada entre dois datasets.
# ------------------------------------------------------------------------------

build_pair_shared_table <- function(
  merged,
  set_genes,
  x_prefix,
  y_prefix
) {

  x_lfc <- paste0(
    "log2FC_",
    x_prefix
  )

  y_lfc <- paste0(
    "log2FC_",
    y_prefix
  )

  x_dir <- paste0(
    "direction_",
    x_prefix
  )

  y_dir <- paste0(
    "direction_",
    y_prefix
  )

  out <- merged[
    GeneID %in% set_genes
  ]

  out[
    ,
    same_direction :=
      get(x_dir) ==
      get(y_dir)
  ]

  out[
    ,
    opposite_direction :=
      get(x_dir) !=
      get(y_dir)
  ]

  out[
    ,
    abs_effect_difference :=
      abs(
        get(x_lfc) -
          get(y_lfc)
      )
  ]

  out
}


# ------------------------------------------------------------------------------
# Resumo direcional de um conjunto compartilhado.
# ------------------------------------------------------------------------------

direction_summary <- function(
  dt,
  x_direction_col,
  y_direction_col
) {

  if (nrow(dt) == 0) {
    return(
      list(
        same = 0L,
        opposite = 0L
      )
    )
  }

  same <- sum(
    dt[[x_direction_col]] ==
      dt[[y_direction_col]],
    na.rm = TRUE
  )

  opposite <- sum(
    dt[[x_direction_col]] !=
      dt[[y_direction_col]],
    na.rm = TRUE
  )

  list(
    same = same,
    opposite = opposite
  )
}


# ------------------------------------------------------------------------------
# Scatter de log2FC para dois datasets no universo comum.
# ------------------------------------------------------------------------------

plot_pairwise_effects <- function(
  dt,
  x_col,
  y_col,
  x_primary_col,
  y_primary_col,
  x_name,
  y_name,
  x_color,
  y_color,
  output_file
) {

  plot_dt <- copy(dt)

  rho <- safe_cor(
    plot_dt[[x_col]],
    plot_dt[[y_col]],
    "spearman"
  )

  # Classificação visual:
  #   - genes não-primary em nenhum dataset ficam cinza;
  #   - genes primary exclusivos recebem a cor da respectiva doença;
  #   - genes primary compartilhados e concordantes ficam roxos;
  #   - compartilhados com direção oposta ficam vermelhos.
  plot_dt[
    ,
    plot_class := fifelse(
      get(x_primary_col) == TRUE &
        get(y_primary_col) == TRUE &
        sign(get(x_col)) == sign(get(y_col)),
      "Shared primary — same direction",
      fifelse(
        get(x_primary_col) == TRUE &
          get(y_primary_col) == TRUE &
          sign(get(x_col)) != sign(get(y_col)),
        "Shared primary — opposite direction",
        fifelse(
          get(x_primary_col) == TRUE,
          paste0(
            x_name,
            " primary only"
          ),
          fifelse(
            get(y_primary_col) == TRUE,
            paste0(
              y_name,
              " primary only"
            ),
            "Neither"
          )
        )
      )
    )
  ]

  class_levels <- c(
    "Neither",
    paste0(
      x_name,
      " primary only"
    ),
    paste0(
      y_name,
      " primary only"
    ),
    "Shared primary — same direction",
    "Shared primary — opposite direction"
  )

  plot_dt[
    ,
    plot_class := factor(
      plot_class,
      levels = class_levels
    )
  ]

  color_values <- c(
    "Neither" = "#BDBDBD",
    setNames(
      x_color,
      paste0(
        x_name,
        " primary only"
      )
    ),
    setNames(
      y_color,
      paste0(
        y_name,
        " primary only"
      )
    ),
    "Shared primary — same direction" = "#444444",
    "Shared primary — opposite direction" = "#D73027"
  )

  alpha_values <- c(
    "Neither" = 0.18,
    setNames(
      0.58,
      paste0(
        x_name,
        " primary only"
      )
    ),
    setNames(
      0.58,
      paste0(
        y_name,
        " primary only"
      )
    ),
    "Shared primary — same direction" = 0.82,
    "Shared primary — opposite direction" = 0.95
  )

  size_values <- c(
    "Neither" = 0.75,
    setNames(
      1.05,
      paste0(
        x_name,
        " primary only"
      )
    ),
    setNames(
      1.05,
      paste0(
        y_name,
        " primary only"
      )
    ),
    "Shared primary — same direction" = 1.45,
    "Shared primary — opposite direction" = 1.70
  )

  p <- ggplot(
    plot_dt,
    aes(
      x = .data[[x_col]],
      y = .data[[y_col]]
    )
  ) +

    geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      color = "black"
    ) +

    geom_vline(
      xintercept = 0,
      linewidth = 0.3,
      color = "black"
    ) +

    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2,
      linewidth = 0.55,
      color = "black"
    ) +

    geom_point(
      aes(
        color = plot_class,
        alpha = plot_class,
        size = plot_class
      )
    ) +

    scale_color_manual(
      values = color_values,
      drop = FALSE
    ) +

    scale_alpha_manual(
      values = alpha_values,
      drop = FALSE
    ) +

    scale_size_manual(
      values = size_values,
      drop = FALSE
    ) +

    guides(
      alpha = "none",
      size = "none"
    ) +

    labs(
      title = paste0(
        x_name,
        " × ",
        y_name,
        " — efeitos no universo comum"
      ),
      subtitle = paste0(
        "Spearman rho = ",
        fmt_num(
          rho,
          3
        ),
        " | cores = DEGs do critério principal"
      ),
      x = paste0(
        "log2FC — ",
        x_name
      ),
      y = paste0(
        "log2FC — ",
        y_name
      ),
      color = "Critério principal"
    ) +

    theme_bw(
      base_size = 12
    ) +

    theme(
      plot.background = element_rect(
        fill = "white",
        color = NA
      ),
      panel.background = element_rect(
        fill = "white",
        color = NA
      ),
      plot.title = element_text(
        face = "bold"
      ),
      legend.position = "right"
    )

  ggsave(
    output_file,
    p,
    width = 8.6,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}


# ==============================================================================
# 5. LER E PADRONIZAR RESULTADOS
# ==============================================================================

for (f in c(
  ad_input,
  hd_input,
  pd_input,
  rosmap_main_input,
  rosmap_adjusted_input
)) {
  assert_file(f)
}

cat(
  "\n============================================================\n"
)
cat(
  "08 — INTEGRAÇÃO FINAL AD–HD–PD + VALIDAÇÃO ROSMAP\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Lendo resultados...\n"
)

AD <- standardize_result(
  fread(
    ad_input,
    na.strings = c(
      "NA",
      ""
    )
  ),
  "AD_GSE125583"
)

HD <- standardize_result(
  fread(
    hd_input,
    na.strings = c(
      "NA",
      ""
    )
  ),
  "HD_GSE64810"
)

PD <- standardize_result(
  fread(
    pd_input,
    na.strings = c(
      "NA",
      ""
    )
  ),
  "PD_GSE68719"
)

ROSMAP_main <- standardize_result(
  fread(
    rosmap_main_input,
    na.strings = c(
      "NA",
      ""
    )
  ),
  "AD_ROSMAP_main"
)

ROSMAP_adjusted <- standardize_result(
  fread(
    rosmap_adjusted_input,
    na.strings = c(
      "NA",
      ""
    )
  ),
  "AD_ROSMAP_adjusted"
)


# ==============================================================================
# 6. UNIVERSOS TESTÁVEIS
# ==============================================================================

AD_universe <- AD$GeneID
HD_universe <- HD$GeneID
PD_universe <- PD$GeneID

AD_HD_universe <- intersect(
  AD_universe,
  HD_universe
)

AD_PD_universe <- intersect(
  AD_universe,
  PD_universe
)

HD_PD_universe <- intersect(
  HD_universe,
  PD_universe
)

common_universe <- Reduce(
  intersect,
  list(
    AD_universe,
    HD_universe,
    PD_universe
  )
)

union_universe <- Reduce(
  union,
  list(
    AD_universe,
    HD_universe,
    PD_universe
  )
)


universe_summary <- data.table(
  metric = c(
    "AD_GSE125583",
    "HD_GSE64810",
    "PD_GSE68719",
    "AD_HD_intersection",
    "AD_PD_intersection",
    "HD_PD_intersection",
    "AD_HD_PD_intersection",
    "AD_HD_PD_union",
    "triple_universe_Jaccard_vs_union"
  ),
  value = c(
    length(AD_universe),
    length(HD_universe),
    length(PD_universe),
    length(AD_HD_universe),
    length(AD_PD_universe),
    length(HD_PD_universe),
    length(common_universe),
    length(union_universe),
    length(common_universe) /
      length(union_universe)
  )
)

fwrite(
  universe_summary,
  universe_summary_output,
  sep = "\t",
  na = "NA"
)


cat(
  "\nUniversos testáveis:\n"
)

cat(
  "  AD: ",
  fmt_int(length(AD_universe)),
  "\n",
  sep = ""
)

cat(
  "  HD: ",
  fmt_int(length(HD_universe)),
  "\n",
  sep = ""
)

cat(
  "  PD: ",
  fmt_int(length(PD_universe)),
  "\n",
  sep = ""
)

cat(
  "  Interseção tripla: ",
  fmt_int(length(common_universe)),
  "\n\n",
  sep = ""
)


if (length(common_universe) != EXPECTED_COMMON_UNIVERSE) {

  warning(
    paste0(
      "Interseção tripla atual = ",
      length(common_universe),
      "; valor esperado da fase anterior = ",
      EXPECTED_COMMON_UNIVERSE,
      ". O script prosseguirá usando o universo observado."
    )
  )
}


# ==============================================================================
# 7. MATRIZ DE EFEITOS NO UNIVERSO TRIPLO COMUM
# ==============================================================================

AD_common <- AD[
  GeneID %in% common_universe,
  .(
    GeneID,
    symbol_AD = symbol,
    log2FC_AD = log2FoldChange,
    stat_AD = stat,
    pvalue_AD = pvalue,
    padj_AD = padj,
    primary_AD = DEG_primary,
    strict_AD = DEG_strict,
    direction_AD = direction
  )
]

HD_common <- HD[
  GeneID %in% common_universe,
  .(
    GeneID,
    symbol_HD = symbol,
    log2FC_HD = log2FoldChange,
    stat_HD = stat,
    pvalue_HD = pvalue,
    padj_HD = padj,
    primary_HD = DEG_primary,
    strict_HD = DEG_strict,
    direction_HD = direction
  )
]

PD_common <- PD[
  GeneID %in% common_universe,
  .(
    GeneID,
    symbol_PD = symbol,
    log2FC_PD = log2FoldChange,
    stat_PD = stat,
    pvalue_PD = pvalue,
    padj_PD = padj,
    primary_PD = DEG_primary,
    strict_PD = DEG_strict,
    direction_PD = direction
  )
]

common_effects <- merge(
  AD_common,
  HD_common,
  by = "GeneID",
  all = FALSE,
  sort = FALSE
)

common_effects <- merge(
  common_effects,
  PD_common,
  by = "GeneID",
  all = FALSE,
  sort = FALSE
)

common_effects[
  ,
  symbol := coalesce_symbol(
    symbol_AD,
    symbol_HD,
    symbol_PD
  )
]

setcolorder(
  common_effects,
  c(
    "GeneID",
    "symbol",
    setdiff(
      names(common_effects),
      c(
        "GeneID",
        "symbol"
      )
    )
  )
)

if (nrow(common_effects) != length(common_universe)) {
  stop(
    "Falha ao construir matriz de efeitos do universo comum."
  )
}

fwrite(
  common_effects,
  all_common_effects_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 8. CONJUNTOS DE DEGs — PRIMARY E STRICT NO UNIVERSO COMUM
# ==============================================================================

AD_primary <- get_gene_set(
  AD,
  "primary",
  common_universe
)

HD_primary <- get_gene_set(
  HD,
  "primary",
  common_universe
)

PD_primary <- get_gene_set(
  PD,
  "primary",
  common_universe
)

AD_strict <- get_gene_set(
  AD,
  "strict",
  common_universe
)

HD_strict <- get_gene_set(
  HD,
  "strict",
  common_universe
)

PD_strict <- get_gene_set(
  PD,
  "strict",
  common_universe
)


# ==============================================================================
# 9. ARTICLE-LIKE NOMINAL — LISTAS COMPLETAS
# ==============================================================================
#
# Essa análise NÃO usa o universo triplo comum de propósito.
# Ela preserva a comparação histórica:
#
#   p < 0,05 e |log2FC| > 1
#
# em cada lista completa do dataset.
# ==============================================================================

AD_nominal_full <- get_gene_set(
  AD,
  "nominal",
  universe = NULL
)

HD_nominal_full <- get_gene_set(
  HD,
  "nominal",
  universe = NULL
)

PD_nominal_full <- get_gene_set(
  PD,
  "nominal",
  universe = NULL
)


# ==============================================================================
# 10. RESUMOS DE VENN / REGIÕES EXCLUSIVAS
# ==============================================================================

regions_primary <- threeway_regions(
  AD_primary,
  HD_primary,
  PD_primary,
  "primary_common_universe"
)

regions_strict <- threeway_regions(
  AD_strict,
  HD_strict,
  PD_strict,
  "strict_common_universe"
)

regions_nominal <- threeway_regions(
  AD_nominal_full,
  HD_nominal_full,
  PD_nominal_full,
  "article_like_nominal_full_lists"
)

regions_all <- rbindlist(
  list(
    regions_primary,
    regions_strict,
    regions_nominal
  )
)

fwrite(
  regions_all,
  threeway_regions_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 11. FIGURAS VENN
# ==============================================================================

plot_threeway_venn(
  regions_primary,
  "Interseção dos DEGs — critério principal",
  paste0(
    "Universo testável comum; FDR < 0,05 e |log2FC| > 0,58; ",
    "listas no universo comum: AD = ",
    length(AD_primary),
    ", HD = ",
    length(HD_primary),
    ", PD = ",
    length(PD_primary)
  ),
  file.path(
    plot_dir,
    "venn_final_primary_AD_HD_PD.png"
  )
)

plot_threeway_venn(
  regions_strict,
  "Interseção dos DEGs — análise de sensibilidade",
  paste0(
    "Universo testável comum; FDR < 0,05 e |log2FC| > 1; ",
    "listas no universo comum: AD = ",
    length(AD_strict),
    ", HD = ",
    length(HD_strict),
    ", PD = ",
    length(PD_strict)
  ),
  file.path(
    plot_dir,
    "venn_final_strict_AD_HD_PD.png"
  )
)

plot_threeway_venn(
  regions_nominal,
  "Interseção dos genes — critério comparativo nominal",
  paste0(
    "p < 0,05 e |log2FC| > 1; listas completas: AD = ",
    length(AD_nominal_full),
    ", HD = ",
    length(HD_nominal_full),
    ", PD = ",
    length(PD_nominal_full)
  ),
  file.path(
    plot_dir,
    "venn_final_article_like_nominal_AD_HD_PD.png"
  )
)


# ==============================================================================
# 12. SOBREPOSIÇÕES PAIRWISE — PRIMARY
# ==============================================================================

pairwise_primary_summary <- rbindlist(
  list(
    pairwise_set_summary(
      AD_primary,
      HD_primary,
      "AD",
      "HD",
      "primary_common_universe"
    ),
    pairwise_set_summary(
      AD_primary,
      PD_primary,
      "AD",
      "PD",
      "primary_common_universe"
    ),
    pairwise_set_summary(
      HD_primary,
      PD_primary,
      "HD",
      "PD",
      "primary_common_universe"
    )
  )
)


# ------------------------------------------------------------------------------
# Tabelas dos genes compartilhados.
# ------------------------------------------------------------------------------

AD_HD_shared_genes <- intersect(
  AD_primary,
  HD_primary
)

AD_PD_shared_genes <- intersect(
  AD_primary,
  PD_primary
)

HD_PD_shared_genes <- intersect(
  HD_primary,
  PD_primary
)

triple_shared_genes <- Reduce(
  intersect,
  list(
    AD_primary,
    HD_primary,
    PD_primary
  )
)


AD_HD_shared <- build_pair_shared_table(
  common_effects,
  AD_HD_shared_genes,
  "AD",
  "HD"
)

AD_PD_shared <- build_pair_shared_table(
  common_effects,
  AD_PD_shared_genes,
  "AD",
  "PD"
)

HD_PD_shared <- build_pair_shared_table(
  common_effects,
  HD_PD_shared_genes,
  "HD",
  "PD"
)


# ------------------------------------------------------------------------------
# Para a interseção tripla, classificar se os três efeitos têm a mesma direção.
# ------------------------------------------------------------------------------

triple_shared <- common_effects[
  GeneID %in% triple_shared_genes
]

triple_shared[
  ,
  all_three_same_direction :=
    direction_AD ==
    direction_HD &
    direction_AD ==
    direction_PD
]

triple_shared[
  ,
  direction_pattern := paste(
    direction_AD,
    direction_HD,
    direction_PD,
    sep = "_"
  )
]


fwrite(
  AD_HD_shared,
  ad_hd_shared_output,
  sep = "\t",
  na = "NA"
)

fwrite(
  AD_PD_shared,
  ad_pd_shared_output,
  sep = "\t",
  na = "NA"
)

fwrite(
  HD_PD_shared,
  hd_pd_shared_output,
  sep = "\t",
  na = "NA"
)

fwrite(
  triple_shared,
  triple_shared_output,
  sep = "\t",
  na = "NA"
)


# ------------------------------------------------------------------------------
# Direção pairwise.
# ------------------------------------------------------------------------------

AD_HD_direction <- direction_summary(
  AD_HD_shared,
  "direction_AD",
  "direction_HD"
)

AD_PD_direction <- direction_summary(
  AD_PD_shared,
  "direction_AD",
  "direction_PD"
)

HD_PD_direction <- direction_summary(
  HD_PD_shared,
  "direction_HD",
  "direction_PD"
)


pairwise_primary_summary[
  comparison == "AD_vs_HD",
  `:=`(
    shared_same_direction = AD_HD_direction$same,
    shared_opposite_direction = AD_HD_direction$opposite
  )
]

pairwise_primary_summary[
  comparison == "AD_vs_PD",
  `:=`(
    shared_same_direction = AD_PD_direction$same,
    shared_opposite_direction = AD_PD_direction$opposite
  )
]

pairwise_primary_summary[
  comparison == "HD_vs_PD",
  `:=`(
    shared_same_direction = HD_PD_direction$same,
    shared_opposite_direction = HD_PD_direction$opposite
  )
]


# ==============================================================================
# 13. CORRELAÇÃO DOS EFEITOS NO UNIVERSO COMUM
# ==============================================================================

effect_correlations <- data.table(
  comparison = c(
    "AD_vs_HD",
    "AD_vs_PD",
    "HD_vs_PD"
  ),
  spearman_all_common = c(
    safe_cor(
      common_effects$log2FC_AD,
      common_effects$log2FC_HD,
      "spearman"
    ),
    safe_cor(
      common_effects$log2FC_AD,
      common_effects$log2FC_PD,
      "spearman"
    ),
    safe_cor(
      common_effects$log2FC_HD,
      common_effects$log2FC_PD,
      "spearman"
    )
  ),
  pearson_all_common = c(
    safe_cor(
      common_effects$log2FC_AD,
      common_effects$log2FC_HD,
      "pearson"
    ),
    safe_cor(
      common_effects$log2FC_AD,
      common_effects$log2FC_PD,
      "pearson"
    ),
    safe_cor(
      common_effects$log2FC_HD,
      common_effects$log2FC_PD,
      "pearson"
    )
  )
)

pairwise_primary_summary <- merge(
  pairwise_primary_summary,
  effect_correlations,
  by = "comparison",
  all.x = TRUE,
  sort = FALSE
)

fwrite(
  pairwise_primary_summary,
  pairwise_summary_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 14. SCATTERS PAIRWISE
# ==============================================================================

plot_pairwise_effects(
  common_effects,
  "log2FC_AD",
  "log2FC_HD",
  "primary_AD",
  "primary_HD",
  "AD",
  "HD",
  COLOR_AD,
  COLOR_HD,
  file.path(
    plot_dir,
    "log2FC_AD_vs_HD_common_universe.png"
  )
)

plot_pairwise_effects(
  common_effects,
  "log2FC_AD",
  "log2FC_PD",
  "primary_AD",
  "primary_PD",
  "AD",
  "PD",
  COLOR_AD,
  COLOR_PD,
  file.path(
    plot_dir,
    "log2FC_AD_vs_PD_common_universe.png"
  )
)

plot_pairwise_effects(
  common_effects,
  "log2FC_HD",
  "log2FC_PD",
  "primary_HD",
  "primary_PD",
  "HD",
  "PD",
  COLOR_HD,
  COLOR_PD,
  file.path(
    plot_dir,
    "log2FC_HD_vs_PD_common_universe.png"
  )
)


# ==============================================================================
# 15. ASSINATURA AD COMPARTILHADA COM PELO MENOS UMA OUTRA DOENÇA
# ==============================================================================
#
# Esses genes formam o conjunto mais diretamente relevante para validar a
# componente AD das assinaturas transdiagnósticas no ROSMAP.
#
# Um gene entra se:
#
#   DEG principal em AD
#   E
#   DEG principal em HD OU PD.
# ==============================================================================

AD_shared_any_genes <- intersect(
  AD_primary,
  union(
    HD_primary,
    PD_primary
  )
)

AD_shared_any <- common_effects[
  GeneID %in% AD_shared_any_genes
]

AD_shared_any[
  ,
  shared_with_HD := primary_HD
]

AD_shared_any[
  ,
  shared_with_PD := primary_PD
]

AD_shared_any[
  ,
  shared_with_both_HD_PD :=
    primary_HD &
    primary_PD
]

AD_shared_any[
  ,
  AD_HD_same_direction :=
    fifelse(
      shared_with_HD,
      direction_AD ==
        direction_HD,
      NA
    )
]

AD_shared_any[
  ,
  AD_PD_same_direction :=
    fifelse(
      shared_with_PD,
      direction_AD ==
        direction_PD,
      NA
    )
]

fwrite(
  AD_shared_any,
  transdiagnostic_ad_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 16. PREPARAR ROSMAP PARA VALIDAÇÃO
# ==============================================================================

ROSMAP_main_small <- ROSMAP_main[
  ,
  .(
    GeneID,
    symbol_ROSMAP_main = symbol,
    log2FC_ROSMAP_main = log2FoldChange,
    pvalue_ROSMAP_main = pvalue,
    padj_ROSMAP_main = padj,
    primary_ROSMAP_main = DEG_primary,
    strict_ROSMAP_main = DEG_strict,
    direction_ROSMAP_main = direction
  )
]

ROSMAP_adjusted_small <- ROSMAP_adjusted[
  ,
  .(
    GeneID,
    symbol_ROSMAP_adjusted = symbol,
    log2FC_ROSMAP_adjusted = log2FoldChange,
    pvalue_ROSMAP_adjusted = pvalue,
    padj_ROSMAP_adjusted = padj,
    primary_ROSMAP_adjusted = DEG_primary,
    strict_ROSMAP_adjusted = DEG_strict,
    direction_ROSMAP_adjusted = direction
  )
]


validation_genes <- merge(
  AD_shared_any,
  ROSMAP_main_small,
  by = "GeneID",
  all.x = TRUE,
  sort = FALSE
)

validation_genes <- merge(
  validation_genes,
  ROSMAP_adjusted_small,
  by = "GeneID",
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------------------------
# Concordância de direção GSE125583 ↔ ROSMAP.
# ------------------------------------------------------------------------------

validation_genes[
  ,
  same_direction_AD_ROSMAP_main :=
    !is.na(
      direction_ROSMAP_main
    ) &
    direction_AD ==
      direction_ROSMAP_main
]

validation_genes[
  ,
  same_direction_AD_ROSMAP_adjusted :=
    !is.na(
      direction_ROSMAP_adjusted
    ) &
    direction_AD ==
      direction_ROSMAP_adjusted
]


# ------------------------------------------------------------------------------
# Para genes que já eram concordantes entre AD e outra doença, testar se ROSMAP
# preserva a MESMA direção transdiagnóstica.
# ------------------------------------------------------------------------------

validation_genes[
  ,
  AD_HD_same_direction_validated_ROSMAP_main :=
    fifelse(
      shared_with_HD &
        AD_HD_same_direction == TRUE &
        !is.na(
          direction_ROSMAP_main
        ),
      direction_ROSMAP_main ==
        direction_AD,
      NA
    )
]

validation_genes[
  ,
  AD_PD_same_direction_validated_ROSMAP_main :=
    fifelse(
      shared_with_PD &
        AD_PD_same_direction == TRUE &
        !is.na(
          direction_ROSMAP_main
        ),
      direction_ROSMAP_main ==
        direction_AD,
      NA
    )
]

validation_genes[
  ,
  AD_HD_same_direction_validated_ROSMAP_adjusted :=
    fifelse(
      shared_with_HD &
        AD_HD_same_direction == TRUE &
        !is.na(
          direction_ROSMAP_adjusted
        ),
      direction_ROSMAP_adjusted ==
        direction_AD,
      NA
    )
]

validation_genes[
  ,
  AD_PD_same_direction_validated_ROSMAP_adjusted :=
    fifelse(
      shared_with_PD &
        AD_PD_same_direction == TRUE &
        !is.na(
          direction_ROSMAP_adjusted
        ),
      direction_ROSMAP_adjusted ==
        direction_AD,
      NA
    )
]


fwrite(
  validation_genes,
  validation_genes_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 17. RESUMO DA VALIDAÇÃO ROSMAP
# ==============================================================================

summarize_validation_set <- function(
  dt,
  set_name,
  subset_fun
) {

  keep <- subset_fun(
    dt
  )

  if (
    length(keep) != nrow(dt)
  ) {
    stop(
      "A função de subset da validação não retornou um vetor lógico ",
      "com o mesmo número de linhas da tabela."
    )
  }

  keep[
    is.na(keep)
  ] <- FALSE

  x <- dt[
    keep == TRUE
  ]

  if (nrow(x) == 0) {

    return(
      data.table(
        signature = set_name,
        N = 0L,
        N_present_ROSMAP_main = 0L,
        N_present_ROSMAP_adjusted = 0L,
        Spearman_AD_vs_ROSMAP_main = NA_real_,
        Spearman_AD_vs_ROSMAP_adjusted = NA_real_,
        same_direction_ROSMAP_main = 0L,
        opposite_direction_ROSMAP_main = 0L,
        primary_significant_ROSMAP_main = 0L,
        same_direction_ROSMAP_adjusted = 0L,
        opposite_direction_ROSMAP_adjusted = 0L,
        primary_significant_ROSMAP_adjusted = 0L
      )
    )
  }

  present_main <- !is.na(
    x$log2FC_ROSMAP_main
  )

  present_adjusted <- !is.na(
    x$log2FC_ROSMAP_adjusted
  )

  same_main <- sum(
    x$same_direction_AD_ROSMAP_main,
    na.rm = TRUE
  )

  same_adjusted <- sum(
    x$same_direction_AD_ROSMAP_adjusted,
    na.rm = TRUE
  )

  data.table(
    signature = set_name,
    N = nrow(x),
    N_present_ROSMAP_main = sum(
      present_main
    ),
    N_present_ROSMAP_adjusted = sum(
      present_adjusted
    ),
    Spearman_AD_vs_ROSMAP_main = safe_cor(
      x$log2FC_AD,
      x$log2FC_ROSMAP_main,
      "spearman"
    ),
    Spearman_AD_vs_ROSMAP_adjusted = safe_cor(
      x$log2FC_AD,
      x$log2FC_ROSMAP_adjusted,
      "spearman"
    ),
    same_direction_ROSMAP_main = same_main,
    opposite_direction_ROSMAP_main = sum(
      present_main
    ) -
      same_main,
    primary_significant_ROSMAP_main = sum(
      x$primary_ROSMAP_main == TRUE,
      na.rm = TRUE
    ),
    same_direction_ROSMAP_adjusted = same_adjusted,
    opposite_direction_ROSMAP_adjusted = sum(
      present_adjusted
    ) -
      same_adjusted,
    primary_significant_ROSMAP_adjusted = sum(
      x$primary_ROSMAP_adjusted == TRUE,
      na.rm = TRUE
    )
  )
}


validation_summary <- rbindlist(
  list(

    summarize_validation_set(
      validation_genes,
      "AD_shared_with_HD_or_PD",
      function(x) {
        rep(
          TRUE,
          nrow(x)
        )
      }
    ),

    summarize_validation_set(
      validation_genes,
      "AD_HD_shared_primary",
      function(x) {
        x$shared_with_HD == TRUE
      }
    ),

    summarize_validation_set(
      validation_genes,
      "AD_PD_shared_primary",
      function(x) {
        x$shared_with_PD == TRUE
      }
    ),

    summarize_validation_set(
      validation_genes,
      "AD_HD_PD_shared_primary",
      function(x) {
        x$shared_with_both_HD_PD == TRUE
      }
    ),

    summarize_validation_set(
      validation_genes,
      "AD_HD_same_direction_primary",
      function(x) {
        x$shared_with_HD == TRUE &
          x$AD_HD_same_direction == TRUE
      }
    ),

    summarize_validation_set(
      validation_genes,
      "AD_PD_same_direction_primary",
      function(x) {
        x$shared_with_PD == TRUE &
          x$AD_PD_same_direction == TRUE
      }
    )
  ),
  fill = TRUE
)

fwrite(
  validation_summary,
  validation_summary_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 18. FIGURA DE VALIDAÇÃO ROSMAP DA ASSINATURA AD TRANSDIAGNÓSTICA
# ==============================================================================

validation_plot_main_output <- file.path(
  plot_dir,
  "ROSMAP_validation_transdiagnostic_AD_main.png"
)

validation_plot_adjusted_output <- file.path(
  plot_dir,
  "ROSMAP_validation_transdiagnostic_AD_adjusted.png"
)


if (nrow(validation_genes) > 0) {

  validation_genes[
    ,
    validation_class_main := fifelse(
      primary_ROSMAP_main == TRUE &
        same_direction_AD_ROSMAP_main == TRUE,
      "Same direction + ROSMAP primary",
      fifelse(
        same_direction_AD_ROSMAP_main == TRUE,
        "Same direction",
        "Opposite direction"
      )
    )
  ]

  validation_genes[
    ,
    validation_class_adjusted := fifelse(
      primary_ROSMAP_adjusted == TRUE &
        same_direction_AD_ROSMAP_adjusted == TRUE,
      "Same direction + ROSMAP primary",
      fifelse(
        same_direction_AD_ROSMAP_adjusted == TRUE,
        "Same direction",
        "Opposite direction"
      )
    )
  ]


  p_val_main <- ggplot(
    validation_genes[
      is.finite(log2FC_AD) &
        is.finite(log2FC_ROSMAP_main)
    ],
    aes(
      x = log2FC_AD,
      y = log2FC_ROSMAP_main,
      color = validation_class_main
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
      alpha = 0.78,
      size = 2
    ) +
    scale_color_manual(
      values = c(
        "Opposite direction" = "#E76F51",
        "Same direction" = "#2A9D8F",
        "Same direction + ROSMAP primary" = "#457B9D"
      )
    ) +
    labs(
      title = "Validação ROSMAP da assinatura AD transdiagnóstica",
      subtitle = paste0(
        "ROSMAP principal | genes AD compartilhados com HD e/ou PD | Spearman rho = ",
        fmt_num(
          safe_cor(
            validation_genes$log2FC_AD,
            validation_genes$log2FC_ROSMAP_main,
            "spearman"
          ),
          3
        )
      ),
      x = "log2FC AD — GSE125583",
      y = "log2FC AD — ROSMAP principal",
      color = "Validação"
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      plot.background = element_rect(
        fill = "white",
        color = NA
      )
    )

  ggsave(
    validation_plot_main_output,
    p_val_main,
    width = 8,
    height = 7,
    dpi = 300,
    bg = "white"
  )


  p_val_adjusted <- ggplot(
    validation_genes[
      is.finite(log2FC_AD) &
        is.finite(log2FC_ROSMAP_adjusted)
    ],
    aes(
      x = log2FC_AD,
      y = log2FC_ROSMAP_adjusted,
      color = validation_class_adjusted
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
      alpha = 0.78,
      size = 2
    ) +
    scale_color_manual(
      values = c(
        "Opposite direction" = "#E76F51",
        "Same direction" = "#2A9D8F",
        "Same direction + ROSMAP primary" = "#457B9D"
      )
    ) +
    labs(
      title = "Validação ROSMAP da assinatura AD transdiagnóstica",
      subtitle = paste0(
        "ROSMAP ajustado | genes AD compartilhados com HD e/ou PD | Spearman rho = ",
        fmt_num(
          safe_cor(
            validation_genes$log2FC_AD,
            validation_genes$log2FC_ROSMAP_adjusted,
            "spearman"
          ),
          3
        )
      ),
      x = "log2FC AD — GSE125583",
      y = "log2FC AD — ROSMAP ajustado",
      color = "Validação"
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      plot.background = element_rect(
        fill = "white",
        color = NA
      )
    )

  ggsave(
    validation_plot_adjusted_output,
    p_val_adjusted,
    width = 8,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}


# ==============================================================================
# 19. RESUMO GERAL DA INTEGRAÇÃO
# ==============================================================================

integration_summary <- data.table(
  metric = c(
    "common_testable_universe",
    "AD_primary",
    "HD_primary",
    "PD_primary",
    "AD_strict",
    "HD_strict",
    "PD_strict",
    "AD_nominal_full",
    "HD_nominal_full",
    "PD_nominal_full",
    "AD_HD_primary_shared",
    "AD_PD_primary_shared",
    "HD_PD_primary_shared",
    "AD_HD_PD_primary_shared",
    "AD_HD_primary_same_direction",
    "AD_HD_primary_opposite_direction",
    "AD_PD_primary_same_direction",
    "AD_PD_primary_opposite_direction",
    "HD_PD_primary_same_direction",
    "HD_PD_primary_opposite_direction",
    "AD_HD_PD_primary_all_same_direction",
    "AD_primary_shared_with_HD_or_PD"
  ),
  value = c(
    length(common_universe),
    length(AD_primary),
    length(HD_primary),
    length(PD_primary),
    length(AD_strict),
    length(HD_strict),
    length(PD_strict),
    length(AD_nominal_full),
    length(HD_nominal_full),
    length(PD_nominal_full),
    length(AD_HD_shared_genes),
    length(AD_PD_shared_genes),
    length(HD_PD_shared_genes),
    length(triple_shared_genes),
    AD_HD_direction$same,
    AD_HD_direction$opposite,
    AD_PD_direction$same,
    AD_PD_direction$opposite,
    HD_PD_direction$same,
    HD_PD_direction$opposite,
    sum(
      triple_shared$all_three_same_direction,
      na.rm = TRUE
    ),
    length(AD_shared_any_genes)
  )
)

fwrite(
  integration_summary,
  integration_summary_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 20. RELATÓRIO TEXTUAL
# ==============================================================================

sink(
  report_output
)

cat(
  "============================================================\n"
)
cat(
  "FASE 01c — 08_finalize_AD_HD_PD_integration_ROSMAP_validation.R\n"
)
cat(
  "INTEGRAÇÃO FINAL AD–HD–PD + VALIDAÇÃO ROSMAP\n"
)
cat(
  "============================================================\n\n"
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "1. DESENHO FINAL\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "AD principal: GSE125583\n"
)
cat(
  "AD validação independente: ROSMAP DLPFC\n"
)
cat(
  "HD principal: GSE64810\n"
)
cat(
  "PD principal: GSE68719\n\n"
)

cat(
  "O ROSMAP NÃO redefine a assinatura transdiagnóstica principal.\n"
)
cat(
  "Ele é usado para testar a reprodutibilidade da componente AD das\n"
)
cat(
  "assinaturas compartilhadas.\n\n"
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "2. UNIVERSOS TESTÁVEIS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "AD: ",
  fmt_int(length(AD_universe)),
  "\n",
  sep = ""
)

cat(
  "HD: ",
  fmt_int(length(HD_universe)),
  "\n",
  sep = ""
)

cat(
  "PD: ",
  fmt_int(length(PD_universe)),
  "\n",
  sep = ""
)

cat(
  "AD ∩ HD: ",
  fmt_int(length(AD_HD_universe)),
  "\n",
  sep = ""
)

cat(
  "AD ∩ PD: ",
  fmt_int(length(AD_PD_universe)),
  "\n",
  sep = ""
)

cat(
  "HD ∩ PD: ",
  fmt_int(length(HD_PD_universe)),
  "\n",
  sep = ""
)

cat(
  "AD ∩ HD ∩ PD: ",
  fmt_int(length(common_universe)),
  "\n\n",
  sep = ""
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "3. CRITÉRIOS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Principal:\n"
)
cat(
  "  BH-FDR < 0,05 e |log2FC| > 0,58\n\n"
)

cat(
  "Estrito:\n"
)
cat(
  "  BH-FDR < 0,05 e |log2FC| > 1\n\n"
)

cat(
  "Comparativo nominal:\n"
)
cat(
  "  p < 0,05 e |log2FC| > 1\n\n"
)

cat(
  "Primary e strict foram calculados no universo triplo comum.\n"
)
cat(
  "O nominal foi calculado nas listas completas por dataset para preservar\n"
)
cat(
  "comparabilidade histórica com a Fase 01.\n\n"
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "4. NÚMERO DE DEGs\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Principal — universo comum:\n"
)

cat(
  "  AD: ",
  fmt_int(length(AD_primary)),
  "\n",
  sep = ""
)

cat(
  "  HD: ",
  fmt_int(length(HD_primary)),
  "\n",
  sep = ""
)

cat(
  "  PD: ",
  fmt_int(length(PD_primary)),
  "\n\n",
  sep = ""
)

cat(
  "Estrito — universo comum:\n"
)

cat(
  "  AD: ",
  fmt_int(length(AD_strict)),
  "\n",
  sep = ""
)

cat(
  "  HD: ",
  fmt_int(length(HD_strict)),
  "\n",
  sep = ""
)

cat(
  "  PD: ",
  fmt_int(length(PD_strict)),
  "\n\n",
  sep = ""
)

cat(
  "Nominal — listas completas:\n"
)

cat(
  "  AD: ",
  fmt_int(length(AD_nominal_full)),
  "\n",
  sep = ""
)

cat(
  "  HD: ",
  fmt_int(length(HD_nominal_full)),
  "\n",
  sep = ""
)

cat(
  "  PD: ",
  fmt_int(length(PD_nominal_full)),
  "\n\n",
  sep = ""
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "5. SOBREPOSIÇÃO PRINCIPAL\n"
)
cat(
  "------------------------------------------------------------\n"
)

for (i in seq_len(nrow(pairwise_primary_summary))) {

  row <- pairwise_primary_summary[i]

  cat(
    row$comparison,
    "\n",
    sep = ""
  )

  cat(
    "  compartilhados: ",
    fmt_int(row$shared),
    "\n",
    sep = ""
  )

  cat(
    "  Jaccard: ",
    fmt_num(row$Jaccard, 4),
    "\n",
    sep = ""
  )

  cat(
    "  mesma direção: ",
    fmt_int(row$shared_same_direction),
    "\n",
    sep = ""
  )

  cat(
    "  direção oposta: ",
    fmt_int(row$shared_opposite_direction),
    "\n",
    sep = ""
  )

  cat(
    "  Spearman de todos os efeitos no universo comum: ",
    fmt_num(row$spearman_all_common, 4),
    "\n\n",
    sep = ""
  )
}


cat(
  "Interseção tripla principal:\n"
)
cat(
  "  AD ∩ HD ∩ PD = ",
  fmt_int(length(triple_shared_genes)),
  "\n",
  sep = ""
)

cat(
  "  todos os três na mesma direção = ",
  fmt_int(
    sum(
      triple_shared$all_three_same_direction,
      na.rm = TRUE
    )
  ),
  "\n\n",
  sep = ""
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "6. VALIDAÇÃO ROSMAP DA COMPONENTE AD TRANSDIAGNÓSTICA\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Genes AD principais compartilhados com HD e/ou PD: ",
  fmt_int(length(AD_shared_any_genes)),
  "\n\n",
  sep = ""
)

for (i in seq_len(nrow(validation_summary))) {

  row <- validation_summary[i]

  cat(
    row$signature,
    "\n",
    sep = ""
  )

  cat(
    "  N: ",
    fmt_int(row$N),
    "\n",
    sep = ""
  )

  cat(
    "  Spearman GSE125583 × ROSMAP main: ",
    fmt_num(
      row$Spearman_AD_vs_ROSMAP_main,
      4
    ),
    "\n",
    sep = ""
  )

  cat(
    "  mesma direção no ROSMAP main: ",
    fmt_int(
      row$same_direction_ROSMAP_main
    ),
    " / ",
    fmt_int(
      row$N_present_ROSMAP_main
    ),
    "\n",
    sep = ""
  )

  cat(
    "  ROSMAP main também primary: ",
    fmt_int(
      row$primary_significant_ROSMAP_main
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Spearman GSE125583 × ROSMAP adjusted: ",
    fmt_num(
      row$Spearman_AD_vs_ROSMAP_adjusted,
      4
    ),
    "\n",
    sep = ""
  )

  cat(
    "  mesma direção no ROSMAP adjusted: ",
    fmt_int(
      row$same_direction_ROSMAP_adjusted
    ),
    " / ",
    fmt_int(
      row$N_present_ROSMAP_adjusted
    ),
    "\n",
    sep = ""
  )

  cat(
    "  ROSMAP adjusted também primary: ",
    fmt_int(
      row$primary_significant_ROSMAP_adjusted
    ),
    "\n\n",
    sep = ""
  )
}


cat(
  "------------------------------------------------------------\n"
)
cat(
  "7. INTERPRETAÇÃO METODOLÓGICA\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "- GSE125583 permanece a coorte principal de AD na integração AD–HD–PD.\n"
)
cat(
  "- ROSMAP é tratado como validação independente de AD e não como quarto\n"
)
cat(
  "  conjunto na definição das interseções principais.\n"
)
cat(
  "- A significância no ROSMAP NÃO é exigida para validar direção, porque o\n"
)
cat(
  "  ROSMAP possui menor estabilidade gene-a-gene e menor poder para o critério\n"
)
cat(
  "  principal que GSE125583.\n"
)
cat(
  "- Portanto, a validação ROSMAP é documentada separadamente em termos de:\n"
)
cat(
  "  presença do GeneID, concordância de direção, correlação dos efeitos e\n"
)
cat(
  "  significância quando presente.\n"
)
cat(
  "- Isso evita descartar uma assinatura transdiagnóstica apenas porque um gene\n"
)
cat(
  "  não cruza o cutoff binário em uma coorte independente.\n\n"
)


cat(
  "------------------------------------------------------------\n"
)
cat(
  "8. ARQUIVOS PRINCIPAIS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Resumo dos universos:\n  ",
  universe_summary_output,
  "\n\n",
  sep = ""
)

cat(
  "Resumo da integração:\n  ",
  integration_summary_output,
  "\n\n",
  sep = ""
)

cat(
  "Resumo pairwise:\n  ",
  pairwise_summary_output,
  "\n\n",
  sep = ""
)

cat(
  "Regiões dos Venns:\n  ",
  threeway_regions_output,
  "\n\n",
  sep = ""
)

cat(
  "Genes compartilhados AD–HD:\n  ",
  ad_hd_shared_output,
  "\n\n",
  sep = ""
)

cat(
  "Genes compartilhados AD–PD:\n  ",
  ad_pd_shared_output,
  "\n\n",
  sep = ""
)

cat(
  "Genes compartilhados HD–PD:\n  ",
  hd_pd_shared_output,
  "\n\n",
  sep = ""
)

cat(
  "Genes compartilhados AD–HD–PD:\n  ",
  triple_shared_output,
  "\n\n",
  sep = ""
)

cat(
  "Assinatura AD compartilhada com HD/PD:\n  ",
  transdiagnostic_ad_output,
  "\n\n",
  sep = ""
)

cat(
  "Validação ROSMAP gene a gene:\n  ",
  validation_genes_output,
  "\n\n",
  sep = ""
)

cat(
  "Resumo da validação ROSMAP:\n  ",
  validation_summary_output,
  "\n\n",
  sep = ""
)


cat(
  "============================================================\n"
)
cat(
  "08_finalize_AD_HD_PD_integration_ROSMAP_validation.R CONCLUÍDO\n"
)
cat(
  "============================================================\n"
)

sink()


# ==============================================================================
# 21. SESSION INFO
# ==============================================================================

sink(
  session_info_output
)

print(
  sessionInfo()
)

sink()


# ==============================================================================
# 22. RESUMO NO TERMINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "INTEGRAÇÃO FINAL AD–HD–PD CONCLUÍDA\n"
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
  "DEGs primary:\n"
)

cat(
  "  AD = ",
  fmt_int(length(AD_primary)),
  "\n",
  sep = ""
)

cat(
  "  HD = ",
  fmt_int(length(HD_primary)),
  "\n",
  sep = ""
)

cat(
  "  PD = ",
  fmt_int(length(PD_primary)),
  "\n\n",
  sep = ""
)

cat(
  "Compartilhados primary:\n"
)

cat(
  "  AD–HD = ",
  fmt_int(length(AD_HD_shared_genes)),
  "\n",
  sep = ""
)

cat(
  "  AD–PD = ",
  fmt_int(length(AD_PD_shared_genes)),
  "\n",
  sep = ""
)

cat(
  "  HD–PD = ",
  fmt_int(length(HD_PD_shared_genes)),
  "\n",
  sep = ""
)

cat(
  "  AD–HD–PD = ",
  fmt_int(length(triple_shared_genes)),
  "\n\n",
  sep = ""
)

cat(
  "Genes AD compartilhados com HD e/ou PD para validação ROSMAP: ",
  fmt_int(length(AD_shared_any_genes)),
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
  "Após revisar esta integração, a próxima etapa pode ser a preparação da\n"
)
cat(
  "referência Mathys 2019 e a deconvolução celular.\n"
)
