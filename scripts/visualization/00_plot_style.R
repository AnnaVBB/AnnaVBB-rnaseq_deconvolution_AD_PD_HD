# 00_plot_style.R | Identidade visual centralizada para figuras do TCC
# Usar: source('scripts/visualization/00_plot_style.R')
# Nunca reutilizar uma cor com outro significado entre familias de graficos.

if (!requireNamespace('ggplot2', quietly = TRUE)) {
  stop('Instale ggplot2 para utilizar o modulo de visualizacao.')
}

PLOT_PALETTES <- list(
  cell = c(
    Ex = '#0072B2', In = '#D55E00', Ast = '#009E73', Oli = '#CC79A7',
    Opc = '#E69F00', Mic = '#56B4E9', End = '#7F3C8D', Per = '#8C564B'
  ),
  cell_group = c(
    Neuronios = '#284B63', Oligodendroglia = '#8F5E99',
    Vascular = '#BA5A31'
  ),
  method = c(
    dtangle = '#264653', BisqueRNA = '#B5179E', MuSiC2 = '#F77F00'
  ),
  disease = c(
    AD = '#1F5673', HD = '#9C6644', PD = '#3A7D44'
  ),
  sex = c(
    M = '#334155', F = '#E11D48'
  )
)

# Reutilizar uma cor entre familias e' proibido. Verifica hex exato.
all_colors <- unlist(PLOT_PALETTES, use.names = FALSE)
if (anyDuplicated(toupper(all_colors))) {
  stop('As paletas utilizam cores hexadecimais repetidas entre categorias.')
}

TCC_STYLE <- list(
  dpi = 300L,
  width_mm = 180,
  height_mm = 120,
  base_size = 11,
  axis_line_width = 0.40,
  histogram_border = '#202020',
  histogram_border_width = 0.28,
  background = 'white'
)

tcc_theme <- function(base_size = TCC_STYLE$base_size) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(color = '#252525', linewidth = TCC_STYLE$axis_line_width),
      axis.ticks = ggplot2::element_line(color = '#252525', linewidth = 0.35),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      plot.title = ggplot2::element_text(face = 'bold', size = ggplot2::rel(1.10)),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(.92)),
      legend.position = 'bottom',
      strip.background = ggplot2::element_rect(fill = '#F4F4F4', color = '#D6D6D6'),
      strip.text = ggplot2::element_text(face = 'bold'),
      panel.grid.major.y = ggplot2::element_line(color = '#EAEAEA', linewidth = 0.25),
      plot.background = ggplot2::element_rect(fill = 'white', color = NA),
      panel.background = ggplot2::element_rect(fill = 'white', color = NA)
    )
}

tcc_scale_color <- function(family, ...) {
  values <- PLOT_PALETTES[[family]]
  if (is.null(values)) stop('Paleta desconhecida: ', family)
  ggplot2::scale_color_manual(values = values, drop = FALSE, ...)
}

tcc_scale_fill <- function(family, ...) {
  values <- PLOT_PALETTES[[family]]
  if (is.null(values)) stop('Paleta desconhecida: ', family)
  ggplot2::scale_fill_manual(values = values, drop = FALSE, ...)
}

tcc_histogram <- function(..., bins = 30L, fill = '#ACBDCE') {
  # Contorno sempre visivel, mesmo com cores claras.
  ggplot2::geom_histogram(
    ..., bins = bins, fill = fill,
    color = TCC_STYLE$histogram_border,
    linewidth = TCC_STYLE$histogram_border_width
  )
}

tcc_save_plot <- function(plot, path, width_mm = TCC_STYLE$width_mm,
                          height_mm = TCC_STYLE$height_mm,
                          dpi = TCC_STYLE$dpi) {
  if (!inherits(plot, 'ggplot')) stop('Objeto nao e ggplot')
  if (!is.finite(width_mm) || width_mm <= 0 ||
      !is.finite(height_mm) || height_mm <= 0) stop('Dimensoes invalidas')
  if (!is.finite(dpi) || dpi < 300) stop('Resolucao inferior a 300 dpi')
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width_mm, height = height_mm,
                  units = 'mm', dpi = dpi, bg = TCC_STYLE$background, limitsize = FALSE)
  invisible(path)
}

tcc_assert_axes <- function(p) {
  # Verificacao de rotulos explicitamente declarados; revisao visual ainda e' necessaria.
  labs <- ggplot2::get_labs(p)
  if (is.null(labs$x) || is.null(labs$y) ||
      !nzchar(as.character(labs$x)) || !nzchar(as.character(labs$y))) {
    stop('Grafico sem rotulos explicitos de eixo X e/ou eixo Y')
  }
  invisible(TRUE)
}
