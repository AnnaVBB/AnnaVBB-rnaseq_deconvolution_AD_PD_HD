# ============================================================================== 
# 08_integrate_AD_HD_PD.R | FASE 01c | VERSAO COMPARATIVA REVISAO 02 (2026-10-20)
# ============================================================================== 
# Esta versao preserva as FUNCOES de Venn, scatter e estatisticas da antiga 01c,
# mas NAO elege previamente ROSMAP ou GSE125583 como AD principal.
#
# Saidas: duas integracoes exploratorias INDEPENDENTES AD-HD-PD,
# cada uma com a outra coorte AD como validacao externa. Inclui criterios
# anatomicos, metadados, robustez de DE e concordancia entre estimadores.
#
# Importante: este script NAO executa MuSiC2, dtangle nem Bisque. Somente
# le suas proporcoes ja estimadas, quando presentes. A decisao AD principal
# depende da revisao dos TRES metodos, nunca apenas do MuSiC2.
#
# Executar na RAIZ do repositorio:
#   Rscript scripts/phase01c_ROSMAP/08_integrate_AD_HD_PD.R
#
# Requer data.table e ggplot2. Nao sobrescreve integracoes historicas.
# ============================================================================== 

required_packages <- c("data.table", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Pacotes nao instalados: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# 0. PARAMETROS GERAIS
ALPHA_FDR <- 0.05
LFC_PRIMARY <- 0.58
LFC_STRICT <- 1.00
P_NOMINAL <- 0.05
LFC_NOMINAL <- 1.00
COLOR_AD <- "#7B2CBF"
COLOR_HD <- "#2783A5"
COLOR_PD <- "#ED4C78"
CANDIDATES <- c("ROSMAP", "GSE125583")
CELL_TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
METHODS <- c("MuSiC2", "dtangle", "Bisque")
REFERENCE_VERSION <- "main24"
MARKER_SET <- "top300"

# Os caminhos abaixo sao os nomes EFETIVOS produzidos por scripts 05 anteriores.
DE_FILES <- c(
  ROSMAP = "results/candidate_ROSMAP/deseq2/ROSMAP_deseq2_results.tsv",
  GSE125583 = "results/candidate_GSE125583/deseq2/GSE125583_deseq2_results.tsv",
  HD = "results/deseq2/HD_deseq2_results.tsv", # GSE64810
  PD = "results/deseq2/PD_deseq2_results.tsv"  # GSE68719
)
DE_ADJUSTED_FILES <- c(
  ROSMAP = "results/candidate_ROSMAP/deseq2/ROSMAP_deseq2_results_adjusted.tsv",
  GSE125583 = "results/candidate_GSE125583/deseq2/GSE125583_deseq2_adjusted_results.tsv"
)

# O arquivo MuSiC2 e produzido por 06_run_MuSiC2.R. Os caminhos de dtangle
# e Bisque NAO foram confirmados. Substituir "" pelos nomes verdadeiros
# quando esses algoritmos forem executados e exportados ao MESMO esquema.
PROPORTION_FILES <- c(
  MuSiC2 = "results/deconvolution/MuSiC2/06_proportions_long_all.tsv",
  dtangle = "",
  Bisque = ""
)
MUSIC2_RUN_SUMMARY <- "results/deconvolution/MuSiC2/06_run_summary.tsv"

COMPARISON_DIR <- "results/phase01c_comparacao_coortes_AD"
OUTPUT_SUBDIR <- "integration_candidate_AD_HD_PD"
dir.create(COMPARISON_DIR, recursive = TRUE, showWarnings = FALSE)

# NA = nenhuma escolha definitiva. Para declarar uma coorte como principal,
# preencher com ROSMAP ou GSE125583 e fornecer justificativa cientifica.
# O script RECUSA a declaracao se nao houver 3 metodos completos e comparaveis.
AD_PRIMARY <- NA_character_
AD_DECISION_JUSTIFICATION <- ""

# Valores extraidos dos relatorios 01b/01c (NAO de resultados dos algoritmos).
COHORT_INFO <- data.table(
  dataset = CANDIDATES,
  brain_region = c("DLPFC", "giro fusiforme"),
  n_AD_report = c(77L, 219L),
  n_Control_report = c(67L, 70L),
  covariates = c("sequencingBatch; msex; RIN; pmi; grupo",
                 "idade; sexo; grupo"),
  metadata_limitations = c(
    "Idade >90 censurada em parte da coorte; forte efeito de lote; metadados tecnicos disponiveis",
    "Uma amostra sem covariaveis completas no modelo ajustado; regiao fora do cortex prefrontal"
  ),
  anatomic_note = c("Maior compatibilidade anatomica com BA9 de HD/PD",
                    "Comparacao regional menos direta com BA9 de HD/PD")
)

# ==============================================================================
# 1. FUNCOES AUXILIARES ORIGINAIS, COM SEGURANCA PARA VALORES NA/NAO FINITOS
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
      is.finite(padj) &
      padj < ALPHA_FDR &
      is.finite(log2FoldChange) &
      abs(log2FoldChange) > LFC_PRIMARY
  ]

  # Critério estrito.
  dt[
    ,
    DEG_strict :=
      is.finite(padj) &
      padj < ALPHA_FDR &
      is.finite(log2FoldChange) &
      abs(log2FoldChange) > LFC_STRICT
  ]

  # Critério nominal utilizado apenas na comparação histórica.
  dt[
    ,
    DEG_nominal :=
      is.finite(pvalue) &
      pvalue < P_NOMINAL &
      is.finite(log2FoldChange) &
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
# 2. LEITURA DOS RESULTADOS DESEQ2
# ==============================================================================
read_de <- function(path, label, mandatory = TRUE) {
  if (!file.exists(path)) {
    if (mandatory) stop(label, ": resultado DESeq2 nao encontrado: ", path)
    warning(label, ": resultado ajustado ausente: ", path,
            "; metricas de robustez/validacao ajustada serao NA.")
    return(NULL)
  }
  x <- fread(path, na.strings = c("NA", ""))
  x <- standardize_result(x, label)
  if (!nrow(x)) stop(label, ": tabela vazia: ", path)
  for (nm in c("log2FoldChange", "pvalue", "padj")) {
    x[, (nm) := suppressWarnings(as.numeric(get(nm)))]
  }
  for (nm in c("pvalue", "padj")) {
    if (any(is.finite(x[[nm]]) & (x[[nm]] < 0 | x[[nm]] > 1))) {
      stop(label, ": valores de ", nm, " fora do intervalo [0,1].")
    }
  }
  x[, `:=`(
    DEG_primary = is.finite(padj) & is.finite(log2FoldChange) &
      padj < ALPHA_FDR & abs(log2FoldChange) > LFC_PRIMARY,
    DEG_strict = is.finite(padj) & is.finite(log2FoldChange) &
      padj < ALPHA_FDR & abs(log2FoldChange) > LFC_STRICT,
    DEG_nominal = is.finite(pvalue) & is.finite(log2FoldChange) &
      pvalue < P_NOMINAL & abs(log2FoldChange) > LFC_NOMINAL,
    testable = is.finite(pvalue) & is.finite(log2FoldChange),
    direction = fifelse(is.finite(log2FoldChange),
                        fifelse(log2FoldChange > 0, "Up",
                                fifelse(log2FoldChange < 0, "Down", "Zero")),
                        NA_character_)
  )]
  x
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), path, sep = "\t", na = "NA")
  invisible(path)
}

require_columns <- function(x, required, description) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(description, ": colunas obrigatorias ausentes: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

if (any(grepl("GSE88719", c(DE_FILES, names(DE_FILES)), fixed = TRUE))) {
  stop("Identificador inesperado GSE88719; PD correto e GSE68719.")
}

DE <- list(
  ROSMAP = read_de(DE_FILES[["ROSMAP"]], "AD_ROSMAP"),
  GSE125583 = read_de(DE_FILES[["GSE125583"]], "AD_GSE125583"),
  HD = read_de(DE_FILES[["HD"]], "HD_GSE64810"),
  PD = read_de(DE_FILES[["PD"]], "PD_GSE68719")
)
DE_ADJUSTED <- list(
  ROSMAP = read_de(DE_ADJUSTED_FILES[["ROSMAP"]], "AD_ROSMAP_adjusted", FALSE),
  GSE125583 = read_de(DE_ADJUSTED_FILES[["GSE125583"]],
                      "AD_GSE125583_adjusted", FALSE)
)

# O conjunto de genes testaveis deve ser auditado de duas maneiras:
#  (a) universo historico = GeneIDs presentes nas tabelas, como no script antigo;
#  (b) universo de efeitos = GeneIDs com pvalue e log2FC finitos nos tres.
# O calculo principal mantem (a) para compatibilidade; (b) e relatado a parte.

# ==============================================================================
# 3. ESTABILIDADE DE EXPRES SAO DIFERENCIAL (PRINCIPAL VS AJUSTADO)
# ==============================================================================
compare_DE_models <- function(main, adjusted, candidate) {
  if (is.null(adjusted)) {
    return(data.table(dataset = candidate,
                      n_shared_testable = NA_integer_,
                      spearman_log2FC = NA_real_,
                      main_primary = sum(main$DEG_primary),
                      adjusted_primary = NA_integer_,
                      n_shared_primary = NA_integer_,
                      jaccard_primary = NA_real_,
                      concordant_shared_direction = NA_real_,
                      status = "PENDENTE: modelo ajustado nao encontrado"))
  }
  a <- main[, .(GeneID, lfc_main = log2FoldChange,
                primary_main = DEG_primary, testable_main = testable)]
  b <- adjusted[, .(GeneID, lfc_adjusted = log2FoldChange,
                    primary_adjusted = DEG_primary, testable_adjusted = testable)]
  w <- merge(a, b, by = "GeneID", all = FALSE)
  w <- w[testable_main & testable_adjusted]
  # data.table moderno: nao usar DT[coluna_logica] sem condicao explicita.
  # Isso tentaria encontrar "primary_main" no ambiente externo e falharia.
  m <- w[primary_main == TRUE, GeneID]
  a_set <- w[primary_adjusted == TRUE, GeneID]
  shared <- intersect(m, a_set)
  both <- w[GeneID %in% shared]
  data.table(dataset = candidate,
             n_shared_testable = nrow(w),
             spearman_log2FC = safe_cor(w$lfc_main, w$lfc_adjusted),
             main_primary = length(m),
             adjusted_primary = length(a_set),
             n_shared_primary = length(shared),
             jaccard_primary = safe_jaccard(m, a_set),
             concordant_shared_direction = if (nrow(both)) {
               mean(sign(both$lfc_main) == sign(both$lfc_adjusted))
             } else NA_real_,
             status = "Calculado a partir dos dois arquivos DESeq2")
}

DE_STABILITY <- rbindlist(lapply(CANDIDATES, function(id) {
  compare_DE_models(DE[[id]], DE_ADJUSTED[[id]], id)
}), fill = TRUE)
write_tsv(DE_STABILITY, file.path(COMPARISON_DIR, "08_DE_model_stability.tsv"))

# Concordancia global de efeitos entre as duas coortes AD (contexto adicional).
# Nenhuma delas e declarada principal a partir deste indicador isolado.
A_ROS <- DE$ROSMAP[, .(GeneID, log2FC_ROSMAP = log2FoldChange,
                       DEG_ROSMAP = DEG_primary, tested_ROSMAP = testable)]
A_GSE <- DE$GSE125583[, .(GeneID, log2FC_GSE125583 = log2FoldChange,
                          DEG_GSE125583 = DEG_primary, tested_GSE125583 = testable)]
AD_cross <- merge(A_ROS, A_GSE, by = "GeneID", all = FALSE)
AD_cross <- AD_cross[tested_ROSMAP & tested_GSE125583]
AD_cross_summary <- data.table(
  n_testable_both = nrow(AD_cross),
  spearman_log2FC = safe_cor(AD_cross$log2FC_ROSMAP,
                            AD_cross$log2FC_GSE125583),
  n_same_direction = sum(sign(AD_cross$log2FC_ROSMAP) ==
                         sign(AD_cross$log2FC_GSE125583), na.rm = TRUE),
  fraction_same_direction = if (nrow(AD_cross)) {
    mean(sign(AD_cross$log2FC_ROSMAP) == sign(AD_cross$log2FC_GSE125583))
  } else NA_real_,
  n_primary_shared = sum(AD_cross$DEG_ROSMAP & AD_cross$DEG_GSE125583)
)
write_tsv(AD_cross_summary,
          file.path(COMPARISON_DIR, "08_between_AD_cohorts_DE_concordance.tsv"))

# ==============================================================================
# 4. ANALISE TRANS DIAGNOSTICA COMPLETA PARA UM CANDIDATO AD
# ==============================================================================
run_candidate <- function(candidate, validation_id) {
  AD <- DE[[candidate]]
  HD <- DE$HD
  PD <- DE$PD
  VAL <- DE[[validation_id]]
  VAL_ADJ <- DE_ADJUSTED[[validation_id]]

  output_dir <- file.path("results", paste0("candidate_", candidate), OUTPUT_SUBDIR)
  tables_dir <- file.path(output_dir, "tables")
  plot_dir <- file.path(output_dir, "plots")
  validation_dir <- file.path(output_dir, paste0("validation_", validation_id))
  for (p in c(output_dir, tables_dir, plot_dir, validation_dir)) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
  cat("\n>>> INTEGRACAO EXPLORATORIA: AD = ", candidate,
      "; AD validacao = ", validation_id, "\n", sep = "")

  # 4.1. Universo historico e universo estatisticamente testavel
  AD_universe <- AD$GeneID
  HD_universe <- HD$GeneID
  PD_universe <- PD$GeneID
  AD_HD_universe <- intersect(AD_universe, HD_universe)
  AD_PD_universe <- intersect(AD_universe, PD_universe)
  HD_PD_universe <- intersect(HD_universe, PD_universe)
  common_universe <- Reduce(intersect, list(AD_universe, HD_universe, PD_universe))
  common_effect_testable <- Reduce(intersect, list(
    AD[testable == TRUE, GeneID],
    HD[testable == TRUE, GeneID],
    PD[testable == TRUE, GeneID]
  ))
  union_universe <- Reduce(union, list(AD_universe, HD_universe, PD_universe))

  if (!length(common_universe)) stop(candidate, ": universo comum vazio.")
  if (!length(common_effect_testable)) {
    warning(candidate, ": nenhum gene testavel nos tres conjuntos.")
  }
  universe_summary <- data.table(
    metric = c("AD", "HD_GSE64810", "PD_GSE68719", "AD_HD_intersection",
               "AD_PD_intersection", "HD_PD_intersection", "AD_HD_PD_intersection",
               "AD_HD_PD_union", "AD_HD_PD_finite_effects_intersection",
               "triple_universe_Jaccard_vs_union"),
    value = c(length(AD_universe), length(HD_universe), length(PD_universe),
              length(AD_HD_universe), length(AD_PD_universe),
              length(HD_PD_universe), length(common_universe),
              length(union_universe), length(common_effect_testable),
              if (length(union_universe)) length(common_universe) / length(union_universe) else NA_real_)
  )
  universe_summary[, AD_dataset := candidate]
  write_tsv(universe_summary, file.path(tables_dir, "final_testable_universe_summary.tsv"))

  # 4.2. Montar matriz conjunta de efeitos e flags na intersecao tripla
  AD_common <- AD[GeneID %in% common_universe,
    .(GeneID, symbol_AD = symbol, log2FC_AD = log2FoldChange,
      stat_AD = stat, pvalue_AD = pvalue, padj_AD = padj,
      primary_AD = DEG_primary, strict_AD = DEG_strict, direction_AD = direction)]
  HD_common <- HD[GeneID %in% common_universe,
    .(GeneID, symbol_HD = symbol, log2FC_HD = log2FoldChange,
      stat_HD = stat, pvalue_HD = pvalue, padj_HD = padj,
      primary_HD = DEG_primary, strict_HD = DEG_strict, direction_HD = direction)]
  PD_common <- PD[GeneID %in% common_universe,
    .(GeneID, symbol_PD = symbol, log2FC_PD = log2FoldChange,
      stat_PD = stat, pvalue_PD = pvalue, padj_PD = padj,
      primary_PD = DEG_primary, strict_PD = DEG_strict, direction_PD = direction)]
  common_effects <- merge(merge(AD_common, HD_common, by = "GeneID", sort = FALSE),
                          PD_common, by = "GeneID", sort = FALSE)
  if (nrow(common_effects) != length(common_universe)) {
    stop(candidate, ": GeneIDs duplicados/ausentes apos o merge.")
  }
  common_effects[, symbol := coalesce_symbol(symbol_AD, symbol_HD, symbol_PD)]
  setcolorder(common_effects, c("GeneID", "symbol", setdiff(names(common_effects),
                                                            c("GeneID", "symbol"))))
  write_tsv(common_effects, file.path(tables_dir,
                                     "final_common_universe_effects_AD_HD_PD.tsv"))

  # 4.3. Mesmo criterio do script historico: primary e strict no universo triplo;
  #      criterio article-like nominal em listas COMPLETAS por coorte.
  AD_primary <- get_gene_set(AD, "primary", common_universe)
  HD_primary <- get_gene_set(HD, "primary", common_universe)
  PD_primary <- get_gene_set(PD, "primary", common_universe)
  AD_strict <- get_gene_set(AD, "strict", common_universe)
  HD_strict <- get_gene_set(HD, "strict", common_universe)
  PD_strict <- get_gene_set(PD, "strict", common_universe)
  AD_nominal_full <- get_gene_set(AD, "nominal")
  HD_nominal_full <- get_gene_set(HD, "nominal")
  PD_nominal_full <- get_gene_set(PD, "nominal")

  regions_primary <- threeway_regions(AD_primary, HD_primary, PD_primary,
                                      "primary_common_universe")
  regions_strict <- threeway_regions(AD_strict, HD_strict, PD_strict,
                                     "strict_common_universe")
  regions_nominal <- threeway_regions(AD_nominal_full, HD_nominal_full,
                                      PD_nominal_full, "article_like_nominal_full_lists")
  regions_all <- rbindlist(list(regions_primary, regions_strict, regions_nominal))
  write_tsv(regions_all, file.path(tables_dir, "final_threeway_regions.tsv"))

  # 4.4. Preservar os tres diagramas Venn do script original.
  plot_threeway_venn(regions_primary,
    paste0("AD (", candidate, ") - HD - PD | DEGs principais"),
    paste0("Universo triplo comum; FDR < 0,05 e |log2FC| > 0,58; ",
           "AD=", length(AD_primary), ", HD=", length(HD_primary),
           ", PD=", length(PD_primary)),
    file.path(plot_dir, "venn_final_primary_AD_HD_PD.png"))
  plot_threeway_venn(regions_strict,
    paste0("AD (", candidate, ") - HD - PD | DEGs estritos"),
    paste0("Universo triplo comum; FDR < 0,05 e |log2FC| > 1; ",
           "AD=", length(AD_strict), ", HD=", length(HD_strict),
           ", PD=", length(PD_strict)),
    file.path(plot_dir, "venn_final_strict_AD_HD_PD.png"))
  plot_threeway_venn(regions_nominal,
    paste0("AD (", candidate, ") - HD - PD | criterio nominal"),
    paste0("Listas completas; p < 0,05 e |log2FC| > 1; ",
           "AD=", length(AD_nominal_full), ", HD=", length(HD_nominal_full),
           ", PD=", length(PD_nominal_full)),
    file.path(plot_dir, "venn_final_article_like_nominal_AD_HD_PD.png"))

  # 4.5. Sobreposicao pairwise, direcao, genes compartilhados e intersecao tripla
  AD_HD_shared_genes <- intersect(AD_primary, HD_primary)
  AD_PD_shared_genes <- intersect(AD_primary, PD_primary)
  HD_PD_shared_genes <- intersect(HD_primary, PD_primary)
  triple_shared_genes <- Reduce(intersect, list(AD_primary, HD_primary, PD_primary))
  AD_HD_shared <- build_pair_shared_table(common_effects, AD_HD_shared_genes, "AD", "HD")
  AD_PD_shared <- build_pair_shared_table(common_effects, AD_PD_shared_genes, "AD", "PD")
  HD_PD_shared <- build_pair_shared_table(common_effects, HD_PD_shared_genes, "HD", "PD")
  triple_shared <- copy(common_effects[GeneID %in% triple_shared_genes])
  triple_shared[, all_three_same_direction :=
                 !is.na(direction_AD) & !is.na(direction_HD) &
                 !is.na(direction_PD) & direction_AD == direction_HD &
                 direction_AD == direction_PD]
  triple_shared[, direction_pattern := paste(direction_AD, direction_HD,
                                              direction_PD, sep = "_")]
  write_tsv(AD_HD_shared, file.path(tables_dir, "final_primary_shared_AD_HD.tsv"))
  write_tsv(AD_PD_shared, file.path(tables_dir, "final_primary_shared_AD_PD.tsv"))
  write_tsv(HD_PD_shared, file.path(tables_dir, "final_primary_shared_HD_PD.tsv"))
  write_tsv(triple_shared, file.path(tables_dir, "final_primary_shared_AD_HD_PD.tsv"))

  AD_HD_direction <- direction_summary(AD_HD_shared, "direction_AD", "direction_HD")
  AD_PD_direction <- direction_summary(AD_PD_shared, "direction_AD", "direction_PD")
  HD_PD_direction <- direction_summary(HD_PD_shared, "direction_HD", "direction_PD")
  pairwise_primary_summary <- rbindlist(list(
    pairwise_set_summary(AD_primary, HD_primary, "AD", "HD", "primary_common_universe"),
    pairwise_set_summary(AD_primary, PD_primary, "AD", "PD", "primary_common_universe"),
    pairwise_set_summary(HD_primary, PD_primary, "HD", "PD", "primary_common_universe")
  ))
  pairwise_primary_summary[, `:=`(
    shared_same_direction = c(AD_HD_direction$same, AD_PD_direction$same,
                              HD_PD_direction$same),
    shared_opposite_direction = c(AD_HD_direction$opposite,
                                  AD_PD_direction$opposite, HD_PD_direction$opposite),
    spearman_all_common = c(
      safe_cor(common_effects$log2FC_AD, common_effects$log2FC_HD),
      safe_cor(common_effects$log2FC_AD, common_effects$log2FC_PD),
      safe_cor(common_effects$log2FC_HD, common_effects$log2FC_PD)),
    pearson_all_common = c(
      safe_cor(common_effects$log2FC_AD, common_effects$log2FC_HD, "pearson"),
      safe_cor(common_effects$log2FC_AD, common_effects$log2FC_PD, "pearson"),
      safe_cor(common_effects$log2FC_HD, common_effects$log2FC_PD, "pearson"))
  )]
  write_tsv(pairwise_primary_summary,
            file.path(tables_dir, "final_pairwise_overlap_summary.tsv"))

  # 4.6. Preservar os 3 scatterplots de efeitos originais.
  plot_pairwise_effects(common_effects, "log2FC_AD", "log2FC_HD",
                        "primary_AD", "primary_HD", "AD", "HD", COLOR_AD,
                        COLOR_HD, file.path(plot_dir,
                                          "log2FC_AD_vs_HD_common_universe.png"))
  plot_pairwise_effects(common_effects, "log2FC_AD", "log2FC_PD",
                        "primary_AD", "primary_PD", "AD", "PD", COLOR_AD,
                        COLOR_PD, file.path(plot_dir,
                                          "log2FC_AD_vs_PD_common_universe.png"))
  plot_pairwise_effects(common_effects, "log2FC_HD", "log2FC_PD",
                        "primary_HD", "primary_PD", "HD", "PD", COLOR_HD,
                        COLOR_PD, file.path(plot_dir,
                                          "log2FC_HD_vs_PD_common_universe.png"))

  # 4.7. Assinatura AD com pelo menos uma das outras doencas.
  AD_shared_any_genes <- intersect(AD_primary, union(HD_primary, PD_primary))
  AD_shared_any <- copy(common_effects[GeneID %in% AD_shared_any_genes])
  AD_shared_any[, `:=`(
    shared_with_HD = primary_HD,
    shared_with_PD = primary_PD,
    shared_with_both_HD_PD = primary_HD & primary_PD,
    AD_HD_same_direction = fifelse(primary_HD,
                                   direction_AD == direction_HD, NA),
    AD_PD_same_direction = fifelse(primary_PD,
                                   direction_AD == direction_PD, NA)
  )]
  write_tsv(AD_shared_any,
            file.path(tables_dir, "final_AD_primary_shared_with_HD_or_PD.tsv"))

  # 4.8. Validacao independente reciproca. Nao confundir ausencia com
  #      direcao oposta: valores nao testaveis ficam NA, nunca FALSE.
  v0 <- VAL[, .(GeneID, symbol_validation_main = symbol,
                log2FC_validation_main = log2FoldChange,
                pvalue_validation_main = pvalue,
                padj_validation_main = padj,
                primary_validation_main = DEG_primary,
                strict_validation_main = DEG_strict,
                direction_validation_main = direction)]
  validation_genes <- merge(AD_shared_any, v0, by = "GeneID", all.x = TRUE,
                            sort = FALSE)
  if (!is.null(VAL_ADJ)) {
    v1 <- VAL_ADJ[, .(GeneID, symbol_validation_adjusted = symbol,
                       log2FC_validation_adjusted = log2FoldChange,
                       pvalue_validation_adjusted = pvalue,
                       padj_validation_adjusted = padj,
                       primary_validation_adjusted = DEG_primary,
                       strict_validation_adjusted = DEG_strict,
                       direction_validation_adjusted = direction)]
    validation_genes <- merge(validation_genes, v1, by = "GeneID", all.x = TRUE,
                              sort = FALSE)
  } else {
    for (nm in c("log2FC_validation_adjusted", "pvalue_validation_adjusted",
                 "padj_validation_adjusted")) {
      validation_genes[, (nm) := NA_real_]
    }
    for (nm in c("primary_validation_adjusted", "strict_validation_adjusted")) {
      validation_genes[, (nm) := NA]
    }
    validation_genes[, direction_validation_adjusted := NA_character_]
  }
  for (model in c("main", "adjusted")) {
    lfc_nm <- paste0("log2FC_validation_", model)
    main_nm <- paste0("primary_validation_", model)
    present_nm <- paste0("present_validation_", model)
    direction_nm <- paste0("same_direction_validation_", model)
    validation_genes[, (present_nm) := is.finite(get(lfc_nm)) &
                       is.finite(log2FC_AD)]
    validation_genes[, (direction_nm) := fifelse(
      get(present_nm), sign(log2FC_AD) == sign(get(lfc_nm)), NA)]
    validation_genes[, (paste0("validated_primary_same_direction_", model)) :=
                       get(direction_nm) == TRUE & get(main_nm) == TRUE]
  }
  write_tsv(validation_genes,
            file.path(validation_dir,
                      paste0(validation_id, "_validation_transdiagnostic_AD_genes.tsv")))

  summarize_validation <- function(x, set_name) {
    out <- data.table(signature = set_name, N = nrow(x))
    for (model in c("main", "adjusted")) {
      lfc_nm <- paste0("log2FC_validation_", model)
      present_nm <- paste0("present_validation_", model)
      dir_nm <- paste0("same_direction_validation_", model)
      primary_nm <- paste0("primary_validation_", model)
      valid_n <- sum(x[[present_nm]], na.rm = TRUE)
      n_concordant <- sum(x[[dir_nm]] == TRUE, na.rm = TRUE)
      out[, (paste0("N_present_", model)) := valid_n]
      out[, (paste0("Spearman_", model)) :=
            safe_cor(x$log2FC_AD, x[[lfc_nm]])]
      out[, (paste0("N_concordant_", model)) := n_concordant]
      out[, (paste0("N_opposite_", model)) :=
            sum(x[[dir_nm]] == FALSE, na.rm = TRUE)]
      out[, (paste0("N_primary_", model)) :=
            sum(x[[primary_nm]] == TRUE & x[[present_nm]], na.rm = TRUE)]
      out[, (paste0("N_primary_concordant_", model)) :=
            sum(x[[primary_nm]] == TRUE & x[[dir_nm]] == TRUE, na.rm = TRUE)]
    }
    out
  }
  val_subsets <- list(
    AD_shared_with_HD_or_PD = validation_genes,
    AD_HD_shared_primary = validation_genes[shared_with_HD == TRUE],
    AD_PD_shared_primary = validation_genes[shared_with_PD == TRUE],
    AD_HD_PD_shared_primary = validation_genes[shared_with_both_HD_PD == TRUE],
    AD_HD_same_direction_primary = validation_genes[AD_HD_same_direction == TRUE],
    AD_PD_same_direction_primary = validation_genes[AD_PD_same_direction == TRUE]
  )
  validation_summary <- rbindlist(lapply(names(val_subsets), function(k) {
    summarize_validation(val_subsets[[k]], k)
  }), use.names = TRUE, fill = TRUE)
  write_tsv(validation_summary,
            file.path(validation_dir, paste0(validation_id, "_validation_summary.tsv")))

  # 4.9. Preservar as duas figuras de validacao (principal e ajustado).
  make_validation_scatter <- function(model) {
    lfc_nm <- paste0("log2FC_validation_", model)
    dir_nm <- paste0("same_direction_validation_", model)
    primary_nm <- paste0("primary_validation_", model)
    w <- copy(validation_genes[
      is.finite(log2FC_AD) & is.finite(get(lfc_nm))])
    w[, class := fifelse(get(dir_nm) == TRUE & get(primary_nm) == TRUE,
                         "Mesma direcao + DEG validacao",
                         fifelse(get(dir_nm) == TRUE, "Mesma direcao",
                                 "Direcao oposta"))]
    p <- ggplot(w, aes(x = log2FC_AD, y = .data[[lfc_nm]], color = class)) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      geom_point(alpha = 0.75, size = 2) +
      scale_color_manual(values = c(
        "Direcao oposta" = "#E76F51",
        "Mesma direcao" = "#2A9D8F",
        "Mesma direcao + DEG validacao" = "#457B9D"), drop = FALSE) +
      labs(title = paste0("Validacao externa AD: ", candidate, " -> ", validation_id),
           subtitle = paste0("Modelo ", model, " | Spearman rho = ",
                             fmt_num(safe_cor(w$log2FC_AD, w[[lfc_nm]]), 3),
                             " | n = ", nrow(w)),
           x = paste0("log2FC AD - ", candidate),
           y = paste0("log2FC AD - ", validation_id, " (", model, ")"),
           color = "Concordancia") +
      theme_bw(base_size = 12) +
      theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(file.path(plot_dir, paste0(validation_id,
                                     "_validation_transdiagnostic_AD_", model, ".png")),
           p, width = 8, height = 7, dpi = 300, bg = "white")
  }
  make_validation_scatter("main")
  if (!is.null(VAL_ADJ)) make_validation_scatter("adjusted")

  # 4.10. Resumo numerico e relatorio textual de cada cenario.
  integration_summary <- data.table(
    metric = c("common_testable_universe_geneids",
               "common_finite_effects_universe",
               "AD_primary", "HD_primary", "PD_primary",
               "AD_strict", "HD_strict", "PD_strict",
               "AD_nominal_full", "HD_nominal_full", "PD_nominal_full",
               "AD_HD_primary_shared", "AD_PD_primary_shared", "HD_PD_primary_shared",
               "AD_HD_PD_primary_shared", "AD_HD_primary_same_direction",
               "AD_HD_primary_opposite_direction", "AD_PD_primary_same_direction",
               "AD_PD_primary_opposite_direction", "HD_PD_primary_same_direction",
               "HD_PD_primary_opposite_direction",
               "AD_HD_PD_primary_all_same_direction",
               "AD_primary_shared_with_HD_or_PD"),
    value = c(length(common_universe), length(common_effect_testable),
              length(AD_primary), length(HD_primary), length(PD_primary),
              length(AD_strict), length(HD_strict), length(PD_strict),
              length(AD_nominal_full), length(HD_nominal_full), length(PD_nominal_full),
              length(AD_HD_shared_genes), length(AD_PD_shared_genes),
              length(HD_PD_shared_genes), length(triple_shared_genes),
              AD_HD_direction$same, AD_HD_direction$opposite,
              AD_PD_direction$same, AD_PD_direction$opposite,
              HD_PD_direction$same, HD_PD_direction$opposite,
              sum(triple_shared$all_three_same_direction, na.rm = TRUE),
              length(AD_shared_any_genes)))
  integration_summary[, `:=`(AD_candidate = candidate,
                              AD_validation = validation_id,
                              analysis_status = "EXPLORATORIA - DECISAO AD PENDENTE")]
  write_tsv(integration_summary,
            file.path(tables_dir, "final_AD_HD_PD_integration_summary.tsv"))

  report <- c(
    "FASE 01c - INTEGRACAO TRANS DIAGNOSTICA EXPLORATORIA",
    paste0("AD candidata: ", candidate),
    paste0("AD validacao externa: ", validation_id),
    "HD: GSE64810 | PD: GSE68719",
    "STATUS: nao ha coorte AD definitiva por este script.", "",
    "CRITERIOS:",
    "Principal: FDR BH < 0,05 e |log2FC| > 0,58 (universo triplo)",
    "Estrito: FDR BH < 0,05 e |log2FC| > 1 (universo triplo)",
    "Nominal: p nao ajustado < 0,05 e |log2FC| > 1 (listas completas)",
    "Universo historico = GeneIDs presentes nas 3 tabelas DESeq2.",
    "Universo finito = GeneIDs com pvalue e log2FC validos nas 3 tabelas.",
    "", "UNIVERSOS:",
    paste0("Historico triplo: ", length(common_universe)),
    paste0("Com efeitos finitos triplos: ", length(common_effect_testable)),
    "", "DEGs E INTERSECOES:",
    paste(paste(integration_summary$metric, integration_summary$value, sep = " = "),
          collapse = "\n"),
    "", "CONCORDANCIAS PAIRWISE:",
    paste(capture.output(print(pairwise_primary_summary)), collapse = "\n"),
    "", "VALIDACAO AD EXTERNA:",
    paste(capture.output(print(validation_summary)), collapse = "\n"),
    "", "LIMITACAO:",
    "Contagens de genes compartilhados nao implicam mesma direcao de regulacao.",
    "Na validacao, ausencia de log2FC valido nao foi classificada como discordancia.",
    "Comparacao por deconvolucao (MuSiC2, dtangle, Bisque) e separada desta integracao.",
    paste0("Saidas: ", output_dir)
  )
  writeLines(report, file.path(output_dir,
                              "08_integration_AD_HD_PD_candidate_report.txt"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()),
             file.path(output_dir, "08_sessionInfo.txt"), useBytes = TRUE)

  data.table(dataset = candidate, validation = validation_id,
             n_common_genes = length(common_universe),
             n_common_finite_effects = length(common_effect_testable),
             n_AD_primary = length(AD_primary),
             n_HD_primary = length(HD_primary),
             n_PD_primary = length(PD_primary),
             n_triple_primary = length(triple_shared_genes),
             n_AD_shared_HD_or_PD = length(AD_shared_any_genes),
             output_dir = output_dir)
}

# Executar AS DUAS candidatas antes da avaliacao de deconvolucao.
integration_runs <- rbindlist(list(
  run_candidate("ROSMAP", "GSE125583"),
  run_candidate("GSE125583", "ROSMAP")
), use.names = TRUE, fill = TRUE)
write_tsv(integration_runs,
          file.path(COMPARISON_DIR, "08_two_candidate_integration_summary.tsv"))

# ==============================================================================
# 5. AVALIACAO MULTICRITERIO DA ESCOLHA DE AD
# ==============================================================================
# A Fase 02 produz as proporcoes: aqui NAO fazemos deconvolucao.
# Para cada algoritmo, o arquivo LONG precisa ter colunas:
# dataset, sample_id, group, cell_type, proportion, reference_version, marker_set
# Valores de group para AD = "AD" ou "Control"; nao misturar rotulos.
# Metodo sem arquivo => PENDENTE (nao significa estimativa zero).

safe_rho <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

safe_ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  a <- x[ok]
  b <- y[ok]
  var_a <- mean((a - mean(a))^2)
  var_b <- mean((b - mean(b))^2)
  cov_ab <- mean((a - mean(a)) * (b - mean(b)))
  denominator <- var_a + var_b + (mean(a) - mean(b))^2
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  2 * cov_ab / denominator
}

read_proportions <- function(method, path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    message(method, ": estimativas indisponiveis; criterio marcado PENDENTE.")
    return(NULL)
  }
  dt <- fread(path, na.strings = c("NA", ""))
  require_columns(dt,
    c("dataset", "sample_id", "group", "cell_type", "proportion",
      "reference_version", "marker_set"), paste0("Proporcoes ", method))
  dt[, `:=`(dataset = as.character(dataset),
             sample_id = trimws(as.character(sample_id)),
             group = as.character(group),
             cell_type = as.character(cell_type),
             reference_version = as.character(reference_version),
             marker_set = as.character(marker_set),
             proportion = suppressWarnings(as.numeric(proportion)))]
  if ("method" %in% names(dt)) {
    if (any(!is.na(dt$method) & as.character(dt$method) != method)) {
      stop(method, ": coluna method inesperada no arquivo ", path)
    }
  }
  dt <- dt[reference_version == REFERENCE_VERSION &
             marker_set == MARKER_SET & dataset %in% CANDIDATES]
  if (!nrow(dt)) {
    warning(method, ": arquivo encontrado, mas sem linhas main24/top300 para AD.")
    return(NULL)
  }
  if (anyNA(dt$sample_id) || any(!nzchar(dt$sample_id)) ||
      anyNA(dt$group) || any(!dt$group %in% c("AD", "Control")) ||
      anyNA(dt$cell_type) || any(!dt$cell_type %in% CELL_TYPES) ||
      any(!is.finite(dt$proportion)) ||
      any(dt$proportion < -1e-8 | dt$proportion > 1 + 1e-8)) {
    stop(method, ": amostras, grupos, tipos celulares ou proporcoes invalidos.")
  }
  duplicates <- dt[, .N, by = .(dataset, sample_id, cell_type)]
  if (any(duplicates$N != 1L)) {
    stop(method, ": mais de uma estimativa por dataset/amostra/tipo.")
  }
  sample_check <- dt[, .(n_cells = uniqueN(cell_type),
                         n_groups = uniqueN(group),
                         sum_proportion = sum(proportion)),
                     by = .(dataset, sample_id)]
  if (any(sample_check$n_cells != length(CELL_TYPES)) ||
      any(sample_check$n_groups != 1L) ||
      any(abs(sample_check$sum_proportion - 1) > 0.02)) {
    stop(method, ": amostra sem 8 tipos, grupo ambiguo ou proporcoes nao somam 1.")
  }
  dt[, method := method]
  dt[, .(dataset, sample_id, group, cell_type, proportion, method)]
}

proportions <- lapply(METHODS, function(m) {
  read_proportions(m, PROPORTION_FILES[[m]])
})
names(proportions) <- METHODS
available_methods <- METHODS[!vapply(proportions, is.null, logical(1))]
if (length(available_methods)) {
  PROPS <- rbindlist(proportions[available_methods], use.names = TRUE)
} else {
  PROPS <- data.table(dataset = character(), sample_id = character(),
                      group = character(), cell_type = character(),
                      proportion = numeric(), method = character())
}

# Disponibilidade por coorte, mesmo que um arquivo contenha apenas ROSMAP.
method_status <- CJ(dataset = CANDIDATES, method = METHODS)
method_status[, `:=`(
  file = PROPORTION_FILES[method],
  available = FALSE,
  n_samples = 0L
)]
for (i in seq_len(nrow(method_status))) {
  ds <- method_status$dataset[i]
  m <- method_status$method[i]
  n <- uniqueN(PROPS[dataset == ds & method == m, sample_id])
  method_status[i, `:=`(available = n > 0L, n_samples = n)]
}
write_tsv(method_status,
          file.path(COMPARISON_DIR, "08_method_availability.tsv"))

if (file.exists(MUSIC2_RUN_SUMMARY)) {
  run_dt <- fread(MUSIC2_RUN_SUMMARY)
  needed <- c("dataset", "reference_version", "marker_set")
  require_columns(run_dt, needed, MUSIC2_RUN_SUMMARY)
  run_dt <- run_dt[dataset %in% CANDIDATES &
                     reference_version == REFERENCE_VERSION &
                     marker_set == MARKER_SET]
  write_tsv(run_dt, file.path(COMPARISON_DIR, "08_MuSiC2_run_audit.tsv"))
}

paired_coverage <- list()
pairwise_metrics <- list()
case_control_effects <- list()

for (ds in CANDIDATES) {
  z <- PROPS[dataset == ds]
  present <- METHODS[METHODS %in% unique(z$method)]
  if (!length(present)) {
    paired_coverage[[ds]] <- data.table(
      dataset = ds, n_methods = 0L, three_methods_available = FALSE,
      n_samples_common = 0L, n_AD_common = 0L, n_Control_common = 0L,
      sample_sets_identical = FALSE, comparability = "PENDENTE")
    next
  }

  sample_sets <- lapply(present, function(m) unique(z[method == m, sample_id]))
  common_samples <- Reduce(intersect, sample_sets)
  identical_sets <- all(vapply(sample_sets,
    function(v) setequal(v, sample_sets[[1L]]), logical(1)))
  z <- z[sample_id %in% common_samples]
  labels <- unique(z[, .(sample_id, group)])
  if (any(labels[, .N, by = sample_id]$N != 1L)) {
    stop(ds, ": grupos discordantes entre algoritmos para a mesma amostra.")
  }
  n_AD <- uniqueN(labels[group == "AD", sample_id])
  n_Control <- uniqueN(labels[group == "Control", sample_id])
  all_three <- length(present) == length(METHODS)
  paired_coverage[[ds]] <- data.table(
    dataset = ds, n_methods = length(present),
    three_methods_available = all_three,
    n_samples_common = length(common_samples),
    n_AD_common = n_AD, n_Control_common = n_Control,
    sample_sets_identical = identical_sets,
    comparability = if (all_three && length(common_samples) &&
                        n_AD > 0L && n_Control > 0L) "TRES_METODOS" else "PENDENTE")
  if (!length(common_samples) || n_AD == 0L || n_Control == 0L) next

  # Comparacao biologica exploratoria: diferenca das medianas AD - controle
  # em cada tipo celular, para as MESMAS amostras entre os metodos presentes.
  group_medians <- z[, .(median_prop = median(proportion)),
                     by = .(dataset, method, cell_type, group)]
  effects_wide <- dcast(group_medians, dataset + method + cell_type ~ group,
                        value.var = "median_prop")
  if (all(c("AD", "Control") %in% names(effects_wide))) {
    effects_wide[, `:=`(delta_median = AD - Control,
                         direction = sign(AD - Control))]
    case_control_effects[[ds]] <- effects_wide
  }
  if (length(present) < 2L) next

  for (pair in combn(present, 2L, simplify = FALSE)) {
    a <- z[method == pair[[1L]],
           .(sample_id, group, cell_type, p1 = proportion)]
    b <- z[method == pair[[2L]],
           .(sample_id, group, cell_type, p2 = proportion)]
    paired <- merge(a, b, by = c("sample_id", "group", "cell_type"),
                    all = FALSE)
    if (nrow(paired) != nrow(a) || nrow(paired) != nrow(b)) {
      stop(ds, ": grade de amostra x tipo celular diferente entre ",
           pair[[1L]], " e ", pair[[2L]], ".")
    }
    scores <- paired[, .(n_samples = .N,
                         spearman = safe_rho(p1, p2),
                         ccc_Lin = safe_ccc(p1, p2),
                         mean_absolute_diff_pp = 100 * mean(abs(p1 - p2)),
                         near_zero_fraction_1 = mean(p1 < 1e-6),
                         near_zero_fraction_2 = mean(p2 < 1e-6)),
                     by = .(cell_type)]
    scores[, `:=`(dataset = ds, method_1 = pair[[1L]],
                  method_2 = pair[[2L]])]
    pairwise_metrics[[paste(ds, paste(pair, collapse = "_"), sep = "_")]] <- scores
  }
}

coverage <- rbindlist(paired_coverage, use.names = TRUE, fill = TRUE)
metrics <- if (length(pairwise_metrics)) {
  rbindlist(pairwise_metrics, fill = TRUE)
} else {
  data.table(dataset = character(), cell_type = character(),
             method_1 = character(), method_2 = character(),
             n_samples = integer(), spearman = numeric(), ccc_Lin = numeric(),
             mean_absolute_diff_pp = numeric(),
             near_zero_fraction_1 = numeric(), near_zero_fraction_2 = numeric())
}
effects <- if (length(case_control_effects)) {
  rbindlist(case_control_effects, fill = TRUE)
} else {
  data.table(dataset = character(), method = character(),
             cell_type = character(), AD = numeric(), Control = numeric(),
             delta_median = numeric(), direction = numeric())
}
write_tsv(coverage, file.path(COMPARISON_DIR, "08_sample_coverage.tsv"))
write_tsv(metrics, file.path(COMPARISON_DIR,
                             "08_method_pairwise_concordance_by_cell_type.tsv"))
write_tsv(effects, file.path(COMPARISON_DIR,
                             "08_AD_control_composition_effects.tsv"))

# Consenso de direcao: exige os TRES estimadores disponiveis para o tipo.
if (nrow(effects)) {
  consensus <- effects[, .(
    n_methods = uniqueN(method),
    all_three_same_direction = .N == length(METHODS) &&
      all(is.finite(direction)) && length(unique(direction)) == 1L,
    signs = paste(paste(method, direction, sep = ":"), collapse = ";")),
    by = .(dataset, cell_type)]
} else {
  consensus <- data.table(dataset = character(), cell_type = character(),
                          n_methods = integer(),
                          all_three_same_direction = logical(), signs = character())
}
write_tsv(consensus, file.path(COMPARISON_DIR,
                               "08_case_control_direction_consensus.tsv"))

# Resumo DESCRITIVO. Medianas de rho ou CCC NAO definem qualidade por si so.
if (nrow(metrics)) {
  summary_metrics <- metrics[, .(
    n_cell_method_pairs = .N,
    n_valid_rho = sum(is.finite(spearman)),
    median_rho = if (any(is.finite(spearman))) median(spearman[is.finite(spearman)])
                 else NA_real_,
    median_ccc = if (any(is.finite(ccc_Lin))) median(ccc_Lin[is.finite(ccc_Lin)])
                 else NA_real_,
    median_absolute_diff_pp = median(mean_absolute_diff_pp)), by = dataset]
} else {
  summary_metrics <- data.table(dataset = character(), n_cell_method_pairs = integer(),
                                n_valid_rho = integer(), median_rho = numeric(),
                                median_ccc = numeric(), median_absolute_diff_pp = numeric())
}
summary_metrics <- merge(data.table(dataset = CANDIDATES), summary_metrics,
                         by = "dataset", all.x = TRUE, sort = FALSE)
write_tsv(summary_metrics,
          file.path(COMPARISON_DIR, "08_summary_method_concordance.tsv"))

consensus_summary <- if (nrow(consensus)) {
  consensus[n_methods == length(METHODS),
            .(n_celltypes_all_three = .N,
              n_direction_consensus = sum(all_three_same_direction)),
            by = dataset]
} else data.table(dataset = character(), n_celltypes_all_three = integer(),
                  n_direction_consensus = integer())

# 5.1. Quadro de decisao que REUNE todas as dimensoes pactuadas.
criteria <- merge(COHORT_INFO, DE_STABILITY, by = "dataset",
                  all.x = TRUE, sort = FALSE)
criteria <- merge(criteria, coverage, by = "dataset", all.x = TRUE, sort = FALSE)
criteria <- merge(criteria, summary_metrics, by = "dataset", all.x = TRUE,
                  sort = FALSE)
criteria <- merge(criteria, consensus_summary, by = "dataset", all.x = TRUE,
                  sort = FALSE)
criteria[, decision_status := fifelse(
  three_methods_available == TRUE &
    !is.na(n_celltypes_all_three) &
    n_celltypes_all_three == length(CELL_TYPES) &
    n_AD_common > 0 & n_Control_common > 0,
  "PRONTA_PARA_REVISAO_MULTICRITERIO", "PENDENTE_DOS_TRES_METODOS")]
write_tsv(criteria,
          file.path(COMPARISON_DIR, "08_criterios_multicriterio_AD.tsv"))

# 5.2. NAO ha decisao automatica nem 'melhor coorte' por rho, CCC ou MuSiC2.
if (!is.na(AD_PRIMARY)) {
  if (!AD_PRIMARY %in% CANDIDATES) stop("AD_PRIMARY invalido: ", AD_PRIMARY)
  if (!nzchar(trimws(AD_DECISION_JUSTIFICATION))) {
    stop("Para declarar AD_PRIMARY, informe AD_DECISION_JUSTIFICATION.")
  }
  if (any(criteria$decision_status != "PRONTA_PARA_REVISAO_MULTICRITERIO")) {
    stop("Nao declarar coorte principal: ainda faltam estimativas comparaveis dos ",
         "3 algoritmos nas duas coortes AD.")
  }
  write_tsv(data.table(
    selected_AD = AD_PRIMARY,
    independent_validation = setdiff(CANDIDATES, AD_PRIMARY),
    justification = AD_DECISION_JUSTIFICATION,
    requires_manual_review = TRUE),
    file.path(COMPARISON_DIR, "08_decisao_documentada_AD.tsv"))
} else {
  writeLines(c(
    "DECISAO AD AINDA NAO DEFINIDA",
    "Ambas as integracoes AD-HD-PD sao EXPLORATORIAS.",
    "A decisao deve contemplar: regiao, metadados, robustez DE e ",
    "concordancia MuSiC2/dtangle/Bisque em amostras equivalentes.",
    "Nao usar apenas convergencia do MuSiC2 ou uma correlacao isolada.",
    "Antes de declarar AD_PRIMARY, inspecionar as proporcoes por tipo celular,",
    "a plausibilidade biologica, zeros, efeito de lote e estabilidade das inferencias.",
    paste0("Tabela principal: ",
           file.path(COMPARISON_DIR, "08_criterios_multicriterio_AD.tsv"))
  ), file.path(COMPARISON_DIR, "08_DECISAO_AD_PENDENTE.txt"), useBytes = TRUE)
}

# 5.3. Auditoria de execucao e de caminhos de saida.
inputs <- data.table(
  kind = c(rep("DE", length(DE_FILES)), rep("DE_ajustado", length(DE_ADJUSTED_FILES)),
           rep("deconvolucao", length(PROPORTION_FILES))),
  id = c(names(DE_FILES), names(DE_ADJUSTED_FILES), names(PROPORTION_FILES)),
  path = c(unname(DE_FILES), unname(DE_ADJUSTED_FILES),
           unname(PROPORTION_FILES)))
inputs[, exists := !is.na(path) & nzchar(path) & file.exists(path)]
write_tsv(inputs, file.path(COMPARISON_DIR, "08_input_manifest.tsv"))
writeLines(capture.output(sessionInfo()),
           file.path(COMPARISON_DIR, "08_sessionInfo.txt"), useBytes = TRUE)
cat("\n============================================================\n")
cat("FASE 01c: duas integracoes exploratorias geradas.\n")
cat("ROSMAP e GSE125583 NAO foram eleitas automaticamente.\n")
cat("PD conferido: GSE68719; HD conferido: GSE64810.\n")
cat("Tabela de criterios: ",
    file.path(COMPARISON_DIR, "08_criterios_multicriterio_AD.tsv"), "\n", sep = "")
cat("Status dos metodos:\n")
print(method_status[, .(dataset, method, available, n_samples)])
cat("============================================================\n")
