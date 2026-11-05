#!/usr/bin/env Rscript
# ==============================================================================
# 09 | COMPARACAO INTEGRADA DE METODOS DE DECONVOLUCAO
# Projeto TCC: Alzheimer (AD), Parkinson (PD) e Huntington (HD)
# Versao: 1.0 | Entrada: resultados previamente obtidos nas etapas 05-08f
# ==============================================================================
# PERGUNTAS DESTA ETAPA
#   1. Os tres algoritmos atribuem proporcoes semelhantes AS MESMAS amostras?
#   2. Quais estimativas estao excessivamente zeradas ou sem variabilidade?
#   3. Que contrastes casos-controles sao observados em cada metodo/coorte?
#   4. As direcoes desses contrastes sao descritivamente consistentes?
#   5. Quanto a escolha dos 8 perfis do dtangle altera as conclusoes?
#   6. A regra de agrupamento PREDEFINIDA e acionada em alguma populacao?
#
# NAO FAZ
#   - Nao reroda MuSiC2 / BisqueRNA / dtangle / 08f.
#   - Nao decide sozinho se ROSMAP ou GSE125583 sera a coorte principal de AD.
#   - Nao elege o metodo 'mais correto' pelo valor-p ou por parecer plausivel.
#   - Nao considera HD e PD independentes (44 controles compartilhados).
#   - Nao usa as estimativas como contagens absolutas de celulas.
#   - Nao confunde FDR por metodo com replicacao metodologica independente.
#
# DESENHO PRINCIPAL (todas as 4 coortes, 8 tipos):
#   MuSiC2    : main24/top300
#   BisqueRNA : main24/top300 (07 v5; use.overlap=FALSE)
#   dtangle   : 08e/pooled_counts, todos os genes, 100 marcadores automaticos
# SENSIBILIDADE: dtangle 08e/mean_nucleus_cpm, mesmos demais parametros.
# NOTA: entradas entre algoritmos nao sao perfeitamente iguais em numero de
# genes; comparacoes sao entre ESTIMATIVAS e nao testes de acuracia emparelhados.
#
# EXECUCAO, na raiz do repositorio:
#   Rscript scripts/phase02_deconvolution/09_comparar_metodos_deconvolucao.R
#
# OPCIONAL — ajustar covariaveis PREDEFINIDAS (nomes EXATOS do metadata 05):
#   DECONV09_COVARS='age,sex,rin' Rscript scripts/.../09_comparar_...R
# Em PowerShell:
#   $env:DECONV09_COVARS = 'age,sex,rin'; Rscript scripts/.../09_...R
# Se nao houver covariaveis configuradas, analises de grupo sao EXPLORATORIAS.
#   DECONV09_BOOT_N=500  # replicacoes do IC bootstrap nao ajustado
#   DECONV09_ROOT=/caminho/para/repositorio
#
# SAIDAS: results/deconvolution/Integrated_Comparison/09/
#    Tabelas TSV, figuras PNG, 09_status.tsv, manifesto e sessionInfo.
# Dependencias: data.table, ggplot2. Nao instala ou atualiza pacotes/renv.
# ==============================================================================

options(stringsAsFactors = FALSE)
REQUIRED_PACKAGES <- c("data.table", "ggplot2")
missing_pkgs <- REQUIRED_PACKAGES[!vapply(REQUIRED_PACKAGES, requireNamespace,
                                           logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Pacotes ausentes: ", paste(missing_pkgs, collapse = ", "),
       ". Instale no ambiente do projeto; o script NAO modifica renv.lock.")
}
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- 0. Parametros e convencoes ---------------------------------------------
ROOT <- Sys.getenv("DECONV09_ROOT", unset = ".")
if (!dir.exists(ROOT)) stop("DECONV09_ROOT nao existe: ", ROOT)
ROOT <- normalizePath(ROOT, mustWork = TRUE)
P <- function(...) file.path(ROOT, ...)
OUTPUT <- P("results", "deconvolution", "Integrated_Comparison", "09")
PLOTS <- file.path(OUTPUT, "plots")
dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)

DATASETS <- c("ROSMAP", "GSE125583", "GSE64810", "GSE68719")
DISEASE <- c(ROSMAP = "AD", GSE125583 = "AD", GSE64810 = "HD", GSE68719 = "PD")
EXPECTED_N <- c(ROSMAP = 144L, GSE125583 = 289L,
                GSE64810 = 69L, GSE68719 = 73L)
METHODS <- c("MuSiC2", "BisqueRNA", "dtangle")
CELL_TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
MAIN_DANGLE_MODE <- trimws(Sys.getenv("DECONV09_DTANGLE_MAIN_MODE", "pooled_counts"))
if (!MAIN_DANGLE_MODE %in% c("pooled_counts", "mean_nucleus_cpm")) {
  stop("DECONV09_DTANGLE_MAIN_MODE deve ser pooled_counts ou mean_nucleus_cpm")
}
ALT_DTANGLE_MODE <- setdiff(c("pooled_counts", "mean_nucleus_cpm"), MAIN_DANGLE_MODE)

SUM_TOL <- 1e-6
ZERO_TOL <- 1e-12
NEAR_ZERO_FRAC <- 0.001       # proporcao <= 0.1%, numero entre 0 e 1
ZERO_WARN_FRAC <- 0.90        # 90% de zeros = alerta tecnico
SD_WARN_PP <- 0.20           # SD abaixo de 0.2 ponto percentual = alerta
LOW_MEDIAN_FRAC <- 0.01      # limiar predefinido para possivel agregacao
N_METHODS_FOR_COLLAPSE <- 2L # medianas abaixo de 1% em >=2 metodos
N_DATASETS_FOR_COLLAPSE <- 3L # ... em >=3 das 4 coortes
SEED <- 20261023L
set.seed(SEED)
BOOT_N <- suppressWarnings(as.integer(Sys.getenv("DECONV09_BOOT_N", "500")))
if (is.na(BOOT_N) || BOOT_N < 100L) stop("DECONV09_BOOT_N deve ser >= 100")

cv_str <- trimws(Sys.getenv("DECONV09_COVARS", ""))
COVARS <- if (nzchar(cv_str)) unique(trimws(strsplit(cv_str, ",", fixed = TRUE)[[1L]])) else character()
COVARS <- COVARS[nzchar(COVARS)]
if (any(COVARS %in% c("sample_id", "group", "is_case", "proportion",
                      "dataset", "method", "cell_type"))) {
  stop("COVARS contem campo reservado. Selecione colunas clinicas nao reservadas.")
}

# Cores constantes nos 4 conjuntos de dados (evita legendas ambiguas).
COL_METHOD <- c(MuSiC2 = "#9B59B6", BisqueRNA = "#0F9D8A", dtangle = "#E87951")
COL_CELL <- c(Ex = "#365DCC", In = "#8251DD", Ast = "#149E79",
              Oli = "#E88B13", Opc = "#EAC028", Mic = "#BE295E",
              End = "#16A6B6", Per = "#82909B")

STATUS <- list()
MANIFEST <- list()
log_status <- function(section, check, status, message, n = NA_integer_) {
  STATUS[[length(STATUS) + 1L]] <<- data.table(
    section = as.character(section), check = as.character(check),
    status = as.character(status), n = as.integer(n), message = as.character(message))
  cat(sprintf("[%s] %s / %s: %s\n", status, section, check, message))
}
write_tsv <- function(d, name) {
  path <- file.path(OUTPUT, name)
  data.table::fwrite(d, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}
register_input <- function(label, file, required = TRUE) {
  exists <- file.exists(file)
  MANIFEST[[length(MANIFEST) + 1L]] <<- data.table(
    input = label, path = file, required = required, exists = exists,
    size_bytes = if (exists) as.numeric(file.info(file)$size) else NA_real_)
  if (!exists && required) stop("Arquivo obrigatorio ausente [", label, "]: ", file)
  if (!exists) log_status("inputs", label, "WARN", "Arquivo opcional ausente")
  invisible(exists)
}
read_tab <- function(file, label, required = TRUE) {
  exists <- register_input(label, file, required)
  if (!exists) return(NULL)
  tryCatch(data.table::fread(file, na.strings = c("NA", "", "NaN")),
           error = function(e) stop("Falha ao ler ", label, " em ", file, ": ", conditionMessage(e)))
}
assert_columns <- function(x, fields, label) {
  absent <- setdiff(fields, names(x))
  if (length(absent)) stop(label, ": faltam colunas: ", paste(absent, collapse = ", "))
}
assert_unique <- function(x, fields, label) {
  if (anyNA(x[, ..fields]) || any(vapply(x[, ..fields],
      function(z) any(trimws(as.character(z)) == "", na.rm = TRUE), logical(1)))) {
    stop(label, ": chaves com NA ou vazias")
  }
  if (anyDuplicated(x, by = fields)) stop(label, ": chaves duplicadas: ", paste(fields, collapse = ","))
}
format_method <- function(x) factor(x, levels = METHODS)
base_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 11),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 10),
          legend.position = "bottom",
          plot.margin = margin(20, 30, 25, 25))
}
save_plot <- function(plot_obj, filename, width = 15, height = 9) {
  target <- file.path(PLOTS, filename)
  problem <- tryCatch({
    ggplot2::ggsave(target, plot = plot_obj, width = width, height = height,
                    dpi = 180, bg = "white", limitsize = FALSE)
    if (!file.exists(target) || is.na(file.info(target)$size) || file.info(target)$size <= 0) {
      stop("Arquivo grafico vazio")
    }
    NA_character_
  }, error = function(e) conditionMessage(e))
  log_status("plot", filename, if (is.na(problem)) "PASS" else "WARN",
             if (is.na(problem)) "Figura salva" else problem)
  invisible(is.na(problem))
}

# ---- 1. Ler metadata canonico da etapa 05 ----------------------------------
cat("\n==== 09 / 1. Metadata da etapa 05 ====\n")
META <- rbindlist(lapply(DATASETS, function(ds) {
  file <- P("data", "processed", "deconvolution", ds, "metadata_deconvolution.tsv")
  m <- read_tab(file, paste0("metadata_", ds))
  assert_columns(m, c("sample_id", "group"), file)
  m[, sample_id := trimws(as.character(sample_id))]
  m[, group := trimws(as.character(group))]
  assert_unique(m, "sample_id", file)
  expected_groups <- c("Control", DISEASE[[ds]])
  if (!setequal(unique(m$group), expected_groups)) {
    stop(ds, ": grupos invalidos: ", paste(unique(m$group), collapse = ", "))
  }
  if (nrow(m) != EXPECTED_N[[ds]]) {
    stop(ds, ": n metadata = ", nrow(m), " mas esperado = ", EXPECTED_N[[ds]])
  }
  if ("dataset" %in% names(m)) m[, dataset := NULL]
  if ("disease" %in% names(m)) m[, disease := NULL]
  m[, `:=`(dataset = ds, disease = DISEASE[[ds]])]
  # Preservar campos clinicos, mesmo que difiram entre coortes.
  setcolorder(m, c("dataset", "disease", "sample_id", "group",
                   setdiff(names(m), c("dataset", "disease", "sample_id", "group"))))
  m
}), use.names = TRUE, fill = TRUE)
assert_unique(META, c("dataset", "sample_id"), "metadata consolidado")

if (length(COVARS)) {
  # Nao selecionar covariaveis oportunisticamente: falhar se nao existem.
  for (ds in DATASETS) {
    md <- META[dataset == ds]
    assert_columns(md, COVARS, paste0("DECONV09_COVARS para ", ds))
    for (cov in COVARS) {
      if (all(is.na(md[[cov]]))) stop(ds, ": covariavel ", cov, " inteiramente ausente")
    }
  }
  log_status("design", "covariaveis", "PASS", paste(COVARS, collapse = ", "))
} else {
  log_status("design", "covariaveis", "WARN",
             "Sem covariaveis: contrastes caso-controle somente exploratorios; nao inferir causalidade")
}

# Controle compartilhado entre as duas coortes — jamais somar n sem deduplicar.
shared_controls <- intersect(
  META[dataset == "GSE64810" & group == "Control", sample_id],
  META[dataset == "GSE68719" & group == "Control", sample_id])
write_tsv(data.table(sample_id = sort(shared_controls)), "09_controles_compartilhados_HD_PD.tsv")
log_status("metadata", "shared_HD_PD", if (length(shared_controls) == 44L) "PASS" else "WARN",
           paste(length(shared_controls), "controles compartilhados (esperado 44)"), length(shared_controls))

# ---- 2. Importar resultados primarios: tres metodos, SEM duplicar variantes -
cat("\n==== 09 / 2. Importando estimativas primarias ====\n")
F06 <- P("results", "deconvolution", "MuSiC2", "06_proportions_long_all.tsv")
F07 <- P("results", "deconvolution", "BisqueRNA", "07_proportions_long_all.tsv")
F08E <- P("results", "deconvolution", "dtangle", "08e_replicacao_orientador",
          "08e_proportions_long_all.tsv")

read_props <- function(path, label, method, reference_mode = NULL) {
  x <- read_tab(path, label)
  assert_columns(x, c("dataset", "sample_id", "group", "cell_type", "proportion"), label)
  if (method == "MuSiC2") {
    assert_columns(x, c("reference_version", "marker_set"), label)
    x <- x[reference_version == "main24" & marker_set == "top300"]
  }
  if (method == "BisqueRNA") {
    if ("reference_version" %in% names(x)) x <- x[reference_version == "main24"]
    if ("marker_set" %in% names(x)) x <- x[marker_set == "top300"]
  }
  if (method == "dtangle") {
    assert_columns(x, "reference_mode", label)
    requested_mode <- as.character(reference_mode)
    if (length(requested_mode) != 1L || !requested_mode %in%
        c(MAIN_DANGLE_MODE, ALT_DTANGLE_MODE)) stop("Modo dtangle nao reconhecido")
    x <- x[reference_mode == requested_mode]
  }
  if (!nrow(x)) stop(label, ": filtragem nao retornou nenhuma amostra")
  x <- copy(x[, .(dataset, sample_id, group, cell_type, proportion)])
  for (col in c("dataset", "sample_id", "group", "cell_type")) {
    set(x, j = col, value = trimws(as.character(x[[col]])))
  }
  x[, proportion := as.numeric(proportion)]
  x[, method := method]
  x
}

# Selecionar referencia principal explicitamente; NUNCA comparar duas variantes
# do mesmo algoritmo como se fossem metodos independentes.
P06 <- read_props(F06, "MuSiC2_main24_top300", "MuSiC2")
P07 <- read_props(F07, "BisqueRNA_main24_top300", "BisqueRNA")
P08 <- read_props(F08E, "dtangle_08e_pooled_counts", "dtangle", MAIN_DANGLE_MODE)
D <- rbindlist(list(P06, P07, P08), use.names = TRUE, fill = TRUE)
rm(P06, P07, P08)

# Checagem de alinhamento SEM normalizar de novo (nunca consertar silenciosamente).
if (any(!D$dataset %in% DATASETS) || any(!D$cell_type %in% CELL_TYPES)) {
  stop("Dataset ou tipo celular inesperado nos arquivos: nao remover linhas silenciosamente")
}
if (anyNA(D$proportion) || any(!is.finite(D$proportion))) stop("Proporcoes faltantes ou nao finitas")
if (any(D$proportion < -SUM_TOL | D$proportion > 1 + SUM_TOL)) {
  stop("Proporcoes fora de [0, 1]")
}
assert_unique(D, c("dataset", "method", "sample_id", "cell_type"), "proporcoes")
if (!setequal(unique(D$dataset), DATASETS) || !setequal(unique(D$method), METHODS)) {
  stop("Coortes ou metodos incompletos nos arquivos long consolidados")
}

# Uma linha por amostra e metodo; todas devem ter exatamente os oito tipos.
map_sum <- D[, .(n_cell_types = .N, sum_fraction = sum(proportion),
                 min_fraction = min(proportion), max_fraction = max(proportion)),
             by = .(dataset, method, sample_id)]
if (any(map_sum$n_cell_types != length(CELL_TYPES))) stop("Algumas amostras nao tem os oito tipos")
if (any(abs(map_sum$sum_fraction - 1) > SUM_TOL)) {
  bad <- map_sum[abs(sum_fraction - 1) > SUM_TOL]
  write_tsv(bad, "09_AUDIT_FALHA_somas.tsv")
  stop("Uma ou mais amostras nao somam 1. Consultar 09_AUDIT_FALHA_somas.tsv")
}

# Os conjuntos de IDs precisam coincidir exatamente com o metadata em cada metodo.
coverage <- D[, .(n_rows = .N, n_samples = uniqueN(sample_id)), by = .(dataset, method)]
for (ds in DATASETS) {
  expected <- sort(META[dataset == ds, sample_id])
  for (mt in METHODS) {
    observed <- sort(D[dataset == ds & method == mt, unique(sample_id)])
    if (!identical(observed, expected)) {
      stop(ds, "/", mt, ": conjunto de amostras difere do metadata 05. ",
           "Faltam: ", paste(head(setdiff(expected, observed), 10), collapse = ","),
           "; extras: ", paste(head(setdiff(observed, expected), 10), collapse = ","))
    }
  }
}

# Conferir group dos arquivos contra metadata: nao permitir rotulo obsoleto.
D <- merge(D, META[, .(dataset, sample_id, meta_group = group, disease)],
           by = c("dataset", "sample_id"), all.x = TRUE, sort = FALSE)
if (anyNA(D$meta_group) || any(D$group != D$meta_group)) {
  stop("Grupos dos arquivos 06/07/08e nao correspondem ao metadata da etapa 05")
}
D[, meta_group := NULL]
setcolorder(D, c("dataset", "disease", "method", "sample_id", "group", "cell_type", "proportion"))
write_tsv(coverage, "09_cobertura_amostras.tsv")
write_tsv(map_sum, "09_auditoria_somas_amostras.tsv")
log_status("input", "3methods_4cohorts", "PASS",
           paste(nrow(D), "linhas; 3 metodos x 4 coortes, oito tipos e ID correto"), nrow(D))

# Gravar LONG harmonizado, para ser lido posteriormente pelo 02_deconvolution.Rmd.
write_tsv(D, "09_proportions_long_harmonized.tsv")

# ---- 3. Diagnostico de estimativas zeradas, concentracao e variabilidade ----
cat("\n==== 09 / 3. Diagnosticos de qualidade por tipo ====\n")
variability <- D[, .(
  n_samples = .N,
  mean_fraction = mean(proportion),
  median_fraction = median(proportion),
  sd_pp = 100 * sd(proportion),
  iqr_pp = 100 * IQR(proportion),
  min_pct = 100 * min(proportion),
  max_pct = 100 * max(proportion),
  zero_fraction = mean(proportion <= ZERO_TOL),
  near_zero_fraction = mean(proportion <= NEAR_ZERO_FRAC),
  distinct_vectors_component_6dec = uniqueN(round(proportion, 6))
), by = .(dataset, method, cell_type)]
variability[, `:=`(
  flag_zero = zero_fraction >= ZERO_WARN_FRAC,
  flag_near_zero = near_zero_fraction >= ZERO_WARN_FRAC,
  flag_low_sd = sd_pp < SD_WARN_PP
)]
variability[, any_technical_flag := flag_zero | flag_near_zero | flag_low_sd]
variability[, technical_interpretation := fifelse(
  any_technical_flag,
  "LIMITADO: zero/variabilidade insuficiente para inferencia deste tipo",
  "Sem alerta simples de zero/SD (NAO e validacao biologica)")]
# Restricoes adicionais fundamentadas nas auditorias 06/07/08, sem inventar
# indice quantitativo de acuracia para as coortes reais.
variability[, method_context := fifelse(method == "MuSiC2",
   "Iteracao nao convergente em 06; conferir 06_run_summary e zeros de In/Opc",
   fifelse(method == "BisqueRNA",
     "07d: forte influencia da referencia; estimativas absolutas nao validadas",
     "08f: vieses por tipo; perfil principal 08e e 100 marcadores"))]
write_tsv(variability, "09_celltype_variability_quality.tsv")

# ---- 4. Concordancia pareada entre metodos -------------------------------
cat("\n==== 09 / 4. Comparacao pareada das mesmas amostras ====\n")
METHOD_PAIRS <- combn(METHODS, 2L, simplify = FALSE)
paired_list <- lapply(METHOD_PAIRS, function(pair) {
  a <- D[method == pair[[1]], .(dataset, sample_id, group, cell_type,
                                 proportion_a = proportion)]
  b <- D[method == pair[[2]], .(dataset, sample_id, group, cell_type,
                                 proportion_b = proportion)]
  z <- merge(a, b, by = c("dataset", "sample_id", "group", "cell_type"),
             all = TRUE, sort = FALSE)
  if (anyNA(z$proportion_a) || anyNA(z$proportion_b)) {
    stop("Comparacao pareada incompleta: ", paste(pair, collapse = " vs "))
  }
  z[, `:=`(method_a = pair[[1]], method_b = pair[[2]])]
  z[, delta_pp := 100 * (proportion_a - proportion_b)]
  z
})
PAIRED <- rbindlist(paired_list, use.names = TRUE)
write_tsv(PAIRED, "09_comparacoes_pareadas_amostras.tsv")

safe_cor <- function(a, b, method = "spearman") {
  if (length(a) < 3 || sd(a) < 1e-12 || sd(b) < 1e-12) return(NA_real_)
  suppressWarnings(stats::cor(a, b, method = method, use = "complete.obs"))
}
PAIR_METRICS <- PAIRED[, .(
  n_samples = .N,
  spearman = safe_cor(proportion_a, proportion_b, "spearman"),
  pearson = safe_cor(proportion_a, proportion_b, "pearson"),
  mae_pp = mean(abs(delta_pp)),
  mean_signed_diff_pp = mean(delta_pp),
  rmse_pp = sqrt(mean(delta_pp^2)),
  max_abs_diff_pp = max(abs(delta_pp))
), by = .(dataset, cell_type, method_a, method_b)]
PAIR_METRICS[, pair := paste(method_a, "vs", method_b)]
PAIR_METRICS[, correlation_interpretable := is.finite(spearman)]
write_tsv(PAIR_METRICS, "09_concordancia_metodos.tsv")

# ---- 5. Contrastes caso-controle EXPLORATORIOS -----------------------------
cat("\n==== 09 / 5. Contrastes caso-controle ====\n")
# ATENCAO: os resultados sao composicionais (cada linha soma 1), os metodos
# compartilham as MESMAS amostras, e as coortes PD/HD compartilham controles.
# Valores-p abaixo nao equivalem a replicacoes independentes e NAO demonstram
# alteracoes causais no numero absoluto de celulas.
bootstrap_delta_ci <- function(case, control, n_boot = BOOT_N) {
  bc <- replicate(n_boot, mean(sample(case, length(case), replace = TRUE)) -
                             mean(sample(control, length(control), replace = TRUE)))
  100 * as.numeric(stats::quantile(bc, probs = c(0.025, 0.975), na.rm = TRUE,
                                    names = FALSE, type = 7))
}
contrasts <- D[, {
  dname <- as.character(disease[[1]])
  control <- proportion[group == "Control"]
  cases <- proportion[group == dname]
  if (length(control) < 5 || length(cases) < 5) {
    stop("Grupo com menos de 5 amostras em ", dataset[[1]], "/", method[[1]], "/", cell_type[[1]])
  }
  ci <- bootstrap_delta_ci(cases, control)
  pv <- tryCatch(suppressWarnings(stats::wilcox.test(cases, control,
                     exact = FALSE, correct = FALSE)$p.value), error = function(e) NA_real_)
  list(n_control = length(control), n_cases = length(cases),
       mean_control_pct = 100 * mean(control),
       mean_case_pct = 100 * mean(cases),
       median_control_pct = 100 * median(control),
       median_case_pct = 100 * median(cases),
       delta_mean_pp = 100 * (mean(cases) - mean(control)),
       delta_median_pp = 100 * (median(cases) - median(control)),
       boot_ci_low_pp = ci[[1]], boot_ci_high_pp = ci[[2]],
       p_wilcoxon_exploratory = pv)
}, by = .(dataset, disease, method, cell_type)]
contrasts[, fdr_bh_wilcoxon_8types := stats::p.adjust(p_wilcoxon_exploratory,
                                                    method = "BH"),
          by = .(dataset, method)]
contrasts[, status_inference := "EXPLORATORIO_NAO_AJUSTADO_PARA_CONFOUNDERS"]
contrasts <- merge(contrasts, variability[, .(dataset, method, cell_type,
    technical_flag = any_technical_flag, zero_fraction, near_zero_fraction, sd_pp)],
  by = c("dataset", "method", "cell_type"), all.x = TRUE, sort = FALSE)
contrasts[, direction := fifelse(delta_mean_pp > 0, "increase",
                          fifelse(delta_mean_pp < 0, "decrease", "zero"))]
contrasts[, interpretation_warning := fifelse(technical_flag,
  "Alerta de zero/baixa variabilidade: nao inferir mudanca deste tipo",
  "Contraste exploratorio; efeitos relativos e possivel confundimento")]
write_tsv(contrasts, "09_case_control_exploratory.tsv")

# Opcional: regressao linear ajustada com covariaveis previamente declaradas.
# A analise sempre registra o N usado e NAO substitui o contraste exploratorio.
# Quando ausentes, produz tabela com status, sem simular um ajuste.
adjusted_list <- list()
if (length(COVARS)) {
  for (ds in DATASETS) {
    md <- copy(META[dataset == ds, c("sample_id", "group", COVARS), with = FALSE])
    md[, is_case := as.integer(as.character(group) == DISEASE[[ds]])]
    for (cov in COVARS) {
      z <- md[[cov]]
      if (is.character(z)) set(md, j = cov, value = factor(z))
    }
    for (mt in METHODS) for (ct in CELL_TYPES) {
      y <- D[dataset == ds & method == mt & cell_type == ct,
             .(sample_id, proportion)]
      dat <- merge(y, md, by = "sample_id", all.x = TRUE, sort = FALSE)
      dat <- dat[stats::complete.cases(dat)]
      n_ctrl <- sum(dat$is_case == 0)
      n_case <- sum(dat$is_case == 1)
      common <- data.table(dataset = ds, method = mt, cell_type = ct,
                           n_used = nrow(dat), n_control = n_ctrl, n_cases = n_case)
      if (n_ctrl < 10L || n_case < 10L || nrow(dat) <= length(COVARS) + 5L) {
        common[, `:=`(status = "SKIP_INSUFFICIENT_COMPLETE_CASES",
                       effect_pp = NA_real_, ci_low_pp = NA_real_,
                       ci_high_pp = NA_real_, p_value = NA_real_)]
      } else {
        formula <- stats::reformulate(c("is_case", COVARS), response = "proportion")
        fit <- tryCatch(stats::lm(formula, data = dat), error = function(e) e)
        if (inherits(fit, "error") || is.na(stats::coef(fit)["is_case"]) ||
            fit$rank < length(stats::coef(fit))) {
          common[, `:=`(status = "SKIP_RANK_OR_MODEL_ERROR",
                         effect_pp = NA_real_, ci_low_pp = NA_real_,
                         ci_high_pp = NA_real_, p_value = NA_real_)]
        } else {
          sm <- summary(fit)$coefficients
          ci <- stats::confint(fit, parm = "is_case", level = 0.95)
          common[, `:=`(status = "PASS_OLS_ADJUSTED_EXPLORATORY",
                         effect_pp = 100 * as.numeric(sm["is_case", "Estimate"]),
                         ci_low_pp = 100 * as.numeric(ci[[1]]),
                         ci_high_pp = 100 * as.numeric(ci[[2]]),
                         p_value = as.numeric(sm["is_case", "Pr(>|t|)"]))]
        }
      }
      adjusted_list[[length(adjusted_list) + 1L]] <- common
    }
  }
  adj <- rbindlist(adjusted_list, fill = TRUE)
  adj[, fdr_bh_8types := stats::p.adjust(p_value, "BH"), by = .(dataset, method)]
  adj[, covariates := paste(COVARS, collapse = ",")]
  write_tsv(adj, "09_case_control_covariate_adjusted.tsv")
  log_status("models", "adjusted", "PASS", "Modelo OLS por tipo; verificar rank/N/modelo")
} else {
  write_tsv(data.table(status = "NOT_RUN_NO_PRESPECIFIED_COVARIATES",
                       note = "Use DECONV09_COVARS=col1,col2; nunca selecionar por valor-p"),
            "09_case_control_covariate_adjusted.tsv")
}

# ---- 6. Concordancia das DIRECOES, sem transformar metodos em 'votos' ----
cat("\n==== 09 / 6. Consistencia direcional ====\n")
robustness <- contrasts[, .(
  n_methods = .N,
  n_increase = sum(direction == "increase"),
  n_decrease = sum(direction == "decrease"),
  n_zero = sum(direction == "zero"),
  n_technical_flags = sum(technical_flag),
  n_fdr05_exploratory = sum(fdr_bh_wilcoxon_8types < 0.05, na.rm = TRUE),
  methods_increase = paste(method[direction == "increase"], collapse = ";"),
  methods_decrease = paste(method[direction == "decrease"], collapse = ";")
), by = .(dataset, disease, cell_type)]
robustness[, direction_status := fifelse(
  n_increase == n_methods, "ALL_INCREASE",
  fifelse(n_decrease == n_methods, "ALL_DECREASE",
          "DISCORDANT_OR_ZERO"))]
robustness[, `:=`(
  biological_validation = "NAO_AVALIADA",
  conclusion = "Concordancia descritiva; metodos correlacionados e auditorias com vieses"
)]
write_tsv(robustness, "09_directional_robustness_descriptive.tsv")

# ---- 7. Regra a priori para colapso (nunca aplica por conveniencia) -------
cat("\n==== 09 / 7. Regra de agregacao registrada a priori ====\n")
# Regra: um tipo so e candidato se sua mediana for <1% em pelo menos DOIS
# metodos distintos dentro da MESMA coorte, e isso acontecer em pelo menos
# TRES coortes. Nao usar a mesma analise de sensibilidade como quarto metodo.
collapse_per_cohort <- variability[, .(
  n_methods_below_1pct = sum(median_fraction < LOW_MEDIAN_FRAC),
  methods_below_1pct = paste(method[median_fraction < LOW_MEDIAN_FRAC], collapse = ";")),
  by = .(dataset, cell_type)]
collapse_per_cohort[, cohort_trigger := n_methods_below_1pct >= N_METHODS_FOR_COLLAPSE]
collapse_global <- collapse_per_cohort[, .(
  n_cohorts_triggered = sum(cohort_trigger),
  coortes = paste(dataset[cohort_trigger], collapse = ";")), by = cell_type]
collapse_global[, rule_triggered := n_cohorts_triggered >= N_DATASETS_FOR_COLLAPSE]
collapse_global[, scope := "SOMENTE_SENSIBILIDADE_SE_APLICAVEL"]
write_tsv(collapse_per_cohort, "09_collapse_rule_by_cohort.tsv")
write_tsv(collapse_global, "09_collapse_rule_global.tsv")

# Agregacao depois do ajuste NAO recupera um tipo zerado pelo solver original.
# Apenas disponibilizar a soma como analise secundaria se a regra ocorreu.
aggregated <- list()
if (collapse_global[cell_type == "In", rule_triggered]) {
  aggregated[[length(aggregated) + 1L]] <- D[cell_type %in% c("Ex", "In"),
    .(proportion = sum(proportion)), by = .(dataset, disease, method, sample_id, group)][
      , aggregate_type := "Neurons_Ex_plus_In"]
}
if (collapse_global[cell_type == "Opc", rule_triggered]) {
  aggregated[[length(aggregated) + 1L]] <- D[cell_type %in% c("Oli", "Opc"),
    .(proportion = sum(proportion)), by = .(dataset, disease, method, sample_id, group)][
      , aggregate_type := "Oligodendroglial_Oli_plus_Opc"]
}
if (length(aggregated)) {
  agg <- rbindlist(aggregated, use.names = TRUE)
  write_tsv(agg, "09_secondary_aggregates_long.tsv")
  log_status("collapse", "secondary", "WARN",
             "Agregacao pos-estimativa acionada; nao corrige fracoes colapsadas no ajuste original")
} else {
  write_tsv(data.table(status = "NOT_TRIGGERED", note = "Regra predefinida nao acionada"),
            "09_secondary_aggregates_long.tsv")
  log_status("collapse", "secondary", "PASS", "Regra nao acionada; oito tipos preservados")
}

# ---- 8. Sensibilidade: referencia dtangle 08e, sem 'quarto voto' ----------
cat("\n==== 09 / 8. Sensibilidade dtangle 08e ====\n")
DALT <- read_props(F08E, "dtangle_08e_mean_nucleus_cpm", "dtangle", ALT_DTANGLE_MODE)
assert_unique(DALT, c("dataset", "sample_id", "cell_type"), "dtangle 08e sensibilidade")
DMAIN <- D[method == "dtangle", .(dataset, sample_id, group, cell_type, pooled = proportion)]
DALT <- DALT[, .(dataset, sample_id, group, cell_type, alternative = proportion)]
SA <- merge(DMAIN, DALT, by = c("dataset", "sample_id", "group", "cell_type"),
            all = TRUE, sort = FALSE)
if (anyNA(SA$pooled) || anyNA(SA$alternative) || nrow(SA) != nrow(DMAIN)) {
  stop("Pareamento incompleto das duas configuracoes de referencia 08e")
}
SA[, delta_pp := 100 * (pooled - alternative)]
ref_mode_summary <- SA[, .(
  n_samples = .N,
  mean_abs_delta_pp = mean(abs(delta_pp)),
  signed_mean_delta_pp = mean(delta_pp),
  max_abs_delta_pp = max(abs(delta_pp)),
  spearman = safe_cor(pooled, alternative),
  mean_pooled_pct = 100 * mean(pooled),
  mean_alt_pct = 100 * mean(alternative)
), by = .(dataset, cell_type)]
write_tsv(ref_mode_summary, "09_dtangle_reference_sensitivity.tsv")
write_tsv(SA, "09_dtangle_reference_sensitivity_paired.tsv")
# Delta dos contrastes caso-controle sob as duas referencias.
SA_contrasts <- SA[, .(
  delta_case_control_pooled_pp = 100 * (mean(pooled[group != "Control"]) -
                                      mean(pooled[group == "Control"])),
  delta_case_control_alternative_pp = 100 * (mean(alternative[group != "Control"]) -
                                           mean(alternative[group == "Control"]))
), by = .(dataset, cell_type)]
SA_contrasts[, direction_stable := (
  sign(delta_case_control_pooled_pp) == sign(delta_case_control_alternative_pp))]
write_tsv(SA_contrasts, "09_dtangle_sensitivity_case_control_direction.tsv")

# ---- 9. Auditorias 06, 07d, 08f: REFERENCIAS, nao combinacao de pontuacoes -
cat("\n==== 09 / 9. Inventario de validacoes preexistentes ====\n")
OPTIONAL <- list(
  MuSiC2_run = P("results", "deconvolution", "MuSiC2", "06_run_summary.tsv"),
  Bisque_run = P("results", "deconvolution", "BisqueRNA", "07_run_summary.tsv"),
  Bisque_C = P("results", "deconvolution", "BisqueRNA_Audit", "07d", "C", "07d_C_metrics_by_type.tsv"),
  dtangle_08f = P("results", "deconvolution", "dtangle", "08f_synthetic_validation", "08f_metrics_by_type.tsv"),
  dtangle_08d = P("results", "deconvolution", "dtangle", "08d_synthetic_validation", "08d_metrics_by_type.tsv")
)
for (nm in names(OPTIONAL)) register_input(nm, OPTIONAL[[nm]], required = FALSE)

# 9.1 Resumos dos ajustes anteriores: registrar seus avisos reais, nao apenas
# uma nota narrativa. Sao OPTIONAL porque a comparacao primaria nao depende
# de execucoes historicas; se existirem, contudo, nao ignorar seus problemas.
if (file.exists(OPTIONAL$MuSiC2_run)) {
  run06 <- data.table::fread(OPTIONAL$MuSiC2_run)
  assert_columns(run06, c("dataset", "reference_version", "marker_set",
                          "converged", "cases", "n_id_not_converged"), "06_run_summary")
  run06 <- run06[reference_version == "main24" & marker_set == "top300"]
  if (setequal(run06$dataset, DATASETS) && nrow(run06) == 4L) {
    run06 <- copy(run06[, .(dataset, converged, cases, n_id_not_converged)])
    run06[, fraction_cases_nonconverged := n_id_not_converged / cases]
    write_tsv(run06, "09_MuSiC2_convergence_context.tsv")
    log_status("historical", "MuSiC2_run", if (all(run06$converged)) "PASS" else "WARN",
               paste(sum(!run06$converged), "coortes com status nao convergente"))
  } else {
    log_status("historical", "MuSiC2_run", "WARN",
               "Resumo 06 nao possui exatamente quatro runs main24/top300")
  }
}
if (file.exists(OPTIONAL$Bisque_run)) {
  run07 <- data.table::fread(OPTIONAL$Bisque_run)
  assert_columns(run07, c("dataset", "degenerate_uniform", "degenerate_constant",
                          "bulk_samples"), "07_run_summary")
  if ("reference_version" %in% names(run07)) run07 <- run07[reference_version == "main24"]
  if ("marker_set" %in% names(run07)) run07 <- run07[marker_set == "top300"]
  if (setequal(run07$dataset, DATASETS) && nrow(run07) == 4L) {
    run07 <- copy(run07[, .(dataset, bulk_samples, degenerate_uniform, degenerate_constant)])
    write_tsv(run07, "09_BisqueRNA_run_degeneracy_context.tsv")
    log_status("historical", "BisqueRNA_run", "PASS",
               "Flags locais registrados (ausencia de flag nao valida calibracao)")
  } else {
    log_status("historical", "BisqueRNA_run", "WARN",
               "Resumo 07 nao possui exatamente quatro runs main24/top300")
  }
}
# 9.2 Resultados sinteticos por TIPO; comparar apenas como auditoria documental,
# pois os procedimentos de treino e estimandos podem ser diferentes.
val_tables <- list()
if (file.exists(OPTIONAL$Bisque_C)) {
  v <- data.table::fread(OPTIONAL$Bisque_C)
  assert_columns(v, c("cell_type", "mae_nucleus_pp", "mae_rna_pp"), "07d_C")
  v <- copy(v[, .(cell_type, mae_nucleus_pp, mae_rna_pp)])
  v[, `:=`(method = "BisqueRNA", validation = "07d_C")]
  val_tables[[length(val_tables)+1L]] <- v
}
if (file.exists(OPTIONAL$dtangle_08f)) {
  v <- data.table::fread(OPTIONAL$dtangle_08f)
  assert_columns(v, c("reference_mode", "cell_type", "mae_nucleus_pp", "mae_rna_pp"), "08f")
  v <- copy(v[reference_mode == MAIN_DANGLE_MODE,
               .(cell_type, mae_nucleus_pp, mae_rna_pp)])
  v[, `:=`(method = "dtangle", validation = paste0("08f_", MAIN_DANGLE_MODE))]
  val_tables[[length(val_tables)+1L]] <- v
}
if (length(val_tables)) {
  validation_table <- rbindlist(val_tables, use.names = TRUE, fill = TRUE)
  setcolorder(validation_table, c("method", "validation", "cell_type", "mae_nucleus_pp", "mae_rna_pp"))
  write_tsv(validation_table, "09_validation_synthetic_context_by_type.tsv")
}

VALIDATION_CONTEXT <- data.table(
  method = METHODS,
  existing_evidence = c(
    "06: nao convergencia/zeros In-Opc; 06f simulacoes MuSiC vs NNLS",
    "07d: ancoragem forte a referencia, efeito do contexto bulk, 36 misturas",
    "08d/08f: validacao LODO; diferenciar fracao nuclear/contribuicao RNA"
  ),
  evidence_strength = c("LIMITACOES_GRAVES", "CALIBRACAO_LIMITADA", "TIPO_DEPENDENTE"),
  use_as_absolute_cell_fraction = c(FALSE, FALSE, FALSE),
  reason = c(
    "Nao convergiu globalmente na configuracao testada; In/Opc zerados",
    "Estimativas influenciadas por composicao referencia/contexto",
    "Erro nuclear por tipo ainda elevado, apesar de sinal de RNA mais bem recuperado"
  )
)
write_tsv(VALIDATION_CONTEXT, "09_validity_context_methods.tsv")

# ---- 10. Tabela de decisao ROSMAP versus GSE125583: SEM decisao automatica -
cat("\n==== 09 / 10. Comparacao entre coortes AD ====\n")
ad_evidence <- rbindlist(lapply(c("ROSMAP", "GSE125583"), function(ds) {
  n_control <- META[dataset == ds & group == "Control", .N]
  n_case <- META[dataset == ds & group == "AD", .N]
  all_cons <- robustness[dataset == ds, sum(direction_status != "DISCORDANT_OR_ZERO")]
  n_flags <- variability[dataset == ds, sum(any_technical_flag)]
  m <- META[dataset == ds]
  tissue_cols <- intersect(names(m), c("tissue_subtype", "tissue_region", "brain_region",
                                      "region", "tissue"))
  region <- if (length(tissue_cols)) {
    paste(vapply(tissue_cols, function(k) {
      paste0(k, "=", paste(unique(na.omit(as.character(m[[k]]))), collapse = ";"))
    }, character(1)), collapse = " | ")
  } else "NAO_DOCUMENTADA_EM_METADATA_05"
  data.table(dataset = ds, disease = "AD", n_controls = n_control,
             n_cases = n_case, n_total = nrow(m), region_metadata = region,
             n_celltypes_all_methods_same_direction = all_cons,
             n_dataset_method_celltype_technical_warnings = n_flags,
             DE_robustness = "PENDENTE_AVALIACAO_FASE01",
             covariate_completeness = if (length(COVARS))
               paste0(round(100 * mean(stats::complete.cases(m[, ..COVARS])), 1), "%")
               else "COVARIAVEIS_NAO_CONFIGURADAS",
             decision = "INDEFINIDA_EXIGE_CONFRONTO_CLINICO_DE_E_RESULTADOS")
}), use.names = TRUE)
write_tsv(ad_evidence, "09_AD_ROSMAP_vs_GSE125583_decision_matrix.tsv")

# ---- 11. Graficos reprodutiveis -------------------------------------------
cat("\n==== 09 / 11. Gerando figuras ====\n")
fmt_cell <- function(x) factor(x, levels = CELL_TYPES)
fmt_ds <- function(x) factor(x, levels = DATASETS)

# 11.1 Calor de frequencia de zeros.
q <- copy(variability)
q[, `:=`(dataset = fmt_ds(dataset), method = format_method(method),
         cell_type = fmt_cell(cell_type), value = 100 * zero_fraction)]
p <- ggplot(q, aes(x = cell_type, y = method, fill = value)) +
  geom_tile(color = "white", linewidth = .45) +
  geom_text(aes(label = sprintf("%.0f", value)), size = 3) +
  facet_wrap(~dataset, ncol = 2) +
  scale_fill_gradient(low = "#F8FAFC", high = "#9B2C2C", limits = c(0,100)) +
  labs(title = "Proporcoes exatamente zero por metodo e coorte",
       subtitle = "Valores em % de amostras; zeros nao provam ausencia biologica",
       x = "Tipo celular", y = "Metodo", fill = "% zeros") + base_theme()
save_plot(p, "09_01_heatmap_zeros.png", 15, 9)

# 11.2 Variabilidade (SD, pp) por metodo/coorte/tipo.
p <- ggplot(q, aes(x = cell_type, y = method, fill = sd_pp)) +
  geom_tile(color = "white", linewidth = .4) +
  geom_text(aes(label = sprintf("%.1f", sd_pp)), size = 3) +
  facet_wrap(~dataset, ncol = 2) +
  scale_fill_gradient(low = "#F8FAFC", high = "#176B73") +
  labs(title = "Variabilidade entre individuos",
       subtitle = "Desvio-padrao das proporcoes em pontos percentuais",
       x = "Tipo celular", y = "Metodo", fill = "SD (pp)") + base_theme()
save_plot(p, "09_02_heatmap_sd.png", 15, 9)

# 11.3 Barras de composicao (media por grupo); nunca apaga ou renormaliza.
comp <- D[, .(mean_pct = 100 * mean(proportion)),
          by = .(dataset, method, group, cell_type)]
check_comp <- comp[, .(sum_means = sum(mean_pct), n_types = .N),
                   by = .(dataset, method, group)]
if (any(abs(check_comp$sum_means - 100) > 1e-4) || any(check_comp$n_types != 8L)) {
  stop("Barras de composicao: medias nao fecham 100%; nunca corrigir silenciosamente")
}
comp[, `:=`(sample_label = paste(dataset, group, sep = " | "),
            method = format_method(method), cell_type = fmt_cell(cell_type))]
comp[, sample_label := factor(sample_label,
             levels = unique(comp[order(match(dataset, DATASETS), group), sample_label]))]
p <- ggplot(comp, aes(x = sample_label, y = mean_pct, fill = cell_type)) +
  geom_col(width = .78, color = "white", linewidth = .18) +
  facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = COL_CELL, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Composicao relativa media por grupo",
       subtitle = "Descritivo: metodos diferentes nao sao escalas biologicas equivalentes",
       x = NULL, y = "Contribuicao estimada (%)", fill = "Tipo") +
  base_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p, "09_03_composicao_media.png", 17, 12)

# 11.4 Boxplots: mesma escala Y por TIPO entre metodos e grupos.
# O painel de cada tipo possui escala propria ("free_y"); nunca comparar
# a altura de facetas de tipos diferentes como se tivessem escala igual.
for (ds in DATASETS) {
  d <- copy(D[dataset == ds]); d[, pct := 100 * proportion]
  d[, `:=`(cell_type = fmt_cell(cell_type), method = format_method(method))]
  disease_label <- DISEASE[[ds]]   # fora da avaliacao data.table (evita factor duplicado)
  d[, group := factor(as.character(group), levels = c("Control", disease_label))]
  p <- ggplot(d, aes(x = group, y = pct, fill = method)) +
    geom_boxplot(position = position_dodge(width = .75), width = .65,
                 outlier.shape = NA, alpha = .77, linewidth = .45) +
    geom_point(aes(color = method),
               position = position_jitterdodge(jitter.width = .12,
                         jitter.height = 0, dodge.width = .75, seed = SEED),
               size = .55, alpha = .22) +
    facet_wrap(~cell_type, ncol = 4, scales = "free_y") +
    scale_fill_manual(values = COL_METHOD, drop = FALSE) +
    scale_color_manual(values = COL_METHOD, drop = FALSE, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(.12, .18))) +
    coord_cartesian(clip = "off") +
    labs(title = paste("Estimativas por metodo e grupo —", ds),
         subtitle = "Cada tipo tem sua escala Y, compartilhada entre grupos e metodos daquele painel",
         x = "Grupo", y = "Proporcao (%)", fill = "Metodo") + base_theme()
  save_plot(p, paste0("09_04_boxplots_", ds, ".png"), 18, 10.5)
}

# 11.5 Concordancia por tipo, incluindo NA na Spearman quando os valores sao constantes.
h <- copy(PAIR_METRICS)
h[, `:=`(dataset = fmt_ds(dataset), cell_type = fmt_cell(cell_type))]
p <- ggplot(h, aes(x = cell_type, y = pair, fill = spearman)) +
  geom_tile(color = "white", linewidth = .35) +
  geom_text(aes(label = ifelse(is.finite(spearman),
                               sprintf("%.2f", spearman), "NA")), size = 3) +
  facet_wrap(~dataset, ncol = 2) +
  scale_fill_gradient2(low = "#A93A46", mid = "#F4F5F7", high = "#19796C",
                       midpoint = 0, limits = c(-1,1), na.value = "#D5D5D5") +
  labs(title = "Correlacao de Spearman entre metodos",
       subtitle = "NA = variabilidade insuficiente; correlacao alta nao valida exatidao",
       x = "Tipo celular", y = NULL, fill = "Spearman") + base_theme()
save_plot(p, "09_05_spearman_metodos.png", 15, 9)

# 11.6 Diferencas pareadas globais das estimativas em pontos percentuais.
p <- ggplot(h, aes(x = cell_type, y = mae_pp, fill = pair)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  facet_wrap(~dataset, ncol = 2) +
  labs(title = "Distancia media absoluta entre estimativas",
       subtitle = "Mesmo paciente e tipo celular; MAE entre metodos (nao versus verdade)",
       x = "Tipo", y = "Diferenca absoluta media (pp)", fill = "Par") +
  base_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p, "09_06_distancia_entre_metodos.png", 16, 9)

# 11.6b Dispersoes pareadas AMOSTRA A AMOSTRA. Sao 12 figuras separadas,
# uma por coorte/par de metodos; nao confundir com medias agregadas.
for (ds in DATASETS) {
  for (pair in METHOD_PAIRS) {
    selected <- copy(PAIRED[dataset == ds & method_a == pair[[1]] & method_b == pair[[2]]])
    selected[, `:=`(pct_a = 100 * proportion_a, pct_b = 100 * proportion_b,
                     cell_type = fmt_cell(cell_type))]
    pp <- ggplot(selected, aes(x = pct_a, y = pct_b)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#7C8794") +
      geom_point(color = "#176B73", alpha = .35, size = .7) +
      facet_wrap(~cell_type, scales = "free", ncol = 4) +
      scale_x_continuous(expand = expansion(mult = c(.07, .13))) +
      scale_y_continuous(expand = expansion(mult = c(.07, .13))) +
      coord_cartesian(clip = "off") +
      labs(title = paste0(ds, ": ", pair[[1]], " vs ", pair[[2]]),
           subtitle = "Cada ponto: uma amostra; escalas diferentes entre tipos; linha = identidade",
           x = paste0(pair[[1]], " (%)"), y = paste0(pair[[2]], " (%)")) + base_theme()
    file_n <- paste0("09_06b_paired_", ds, "_", pair[[1]], "_", pair[[2]], ".png")
    save_plot(pp, file_n, width = 17, height = 10.5)
  }
}

# 11.7 Contraste caso-controle nao ajustado, delta das medias, IC bootstrap.
e <- copy(contrasts)
e[, `:=`(dataset = fmt_ds(dataset), cell_type = fmt_cell(cell_type),
          method = format_method(method))]
p <- ggplot(e, aes(x = cell_type, y = delta_mean_pp, color = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7C8794") +
  geom_linerange(aes(ymin = boot_ci_low_pp, ymax = boot_ci_high_pp),
                 position = position_dodge(width = .58), linewidth = .6) +
  geom_point(aes(shape = technical_flag),
             position = position_dodge(width = .58), size = 2.1) +
  facet_wrap(~dataset, ncol = 2, scales = "free_y") +
  scale_color_manual(values = COL_METHOD) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                     labels = c(`FALSE` = "Sem alerta simples", `TRUE` = "Alerta tecnico")) +
  labs(title = "Mudanca media casos - controles",
       subtitle = "IC bootstrap 95% nao ajustado; nao inferir quantidade absoluta de celulas",
       x = "Tipo celular", y = "Diferenca relativa (pp)", color = "Metodo", shape = "QC") +
  base_theme()
save_plot(p, "09_07_case_control_delta_bootstrap.png", 17, 10)

# 11.8 Sensibilidade da construcao dos oito perfis dtangle.
s <- copy(ref_mode_summary)
s[, `:=`(dataset = fmt_ds(dataset), cell_type = fmt_cell(cell_type))]
p <- ggplot(s, aes(x = cell_type, y = mean_abs_delta_pp, fill = cell_type)) +
  geom_col(width = .72) + facet_wrap(~dataset, ncol = 2) +
  scale_fill_manual(values = COL_CELL, guide = "none") +
  labs(title = "Sensibilidade a forma de construir os oito perfis do dtangle",
       subtitle = "MAE entre pooled_counts e mean_nucleus_cpm para a MESMA amostra",
       x = "Tipo", y = "Diferenca media (pp)") + base_theme()
save_plot(p, "09_08_dtangle_modos_referencia.png", 16, 9)

# 11.9 Matriz de acordo direcional (apenas descritiva).
r <- copy(robustness)
r[, `:=`(dataset = fmt_ds(dataset), cell_type = fmt_cell(cell_type))]
p <- ggplot(r, aes(x = cell_type, y = dataset, fill = direction_status)) +
  geom_tile(color = "white", linewidth = .6) +
  geom_text(aes(label = paste0(n_increase, "+ / ", n_decrease, "-")), size = 3) +
  scale_fill_manual(values = c(ALL_INCREASE = "#56A38D", ALL_DECREASE = "#D7787A",
                               DISCORDANT_OR_ZERO = "#C8CED7")) +
  labs(title = "Concordancia descritiva das direcoes caso-controle",
       subtitle = "Tres algoritmos analisam as mesmas amostras; isso NAO e replicacao independente",
       x = "Tipo celular", y = "Coorte", fill = "Direcao") + base_theme()
save_plot(p, "09_09_direcoes_caso_controle.png", 15, 6.5)

# ---- 12. Resumo de interpretacao e rastreabilidade ------------------------
cat("\n==== 09 / 12. Finalizacao e status ====\n")
notes <- c(
  "09 COMPARACAO INTEGRADA DE DECONVOLUCAO - NOTAS DE INTERPRETACAO",
  "=================================================================",
  "Entrada principal: 06 MuSiC2/main24/top300; 07 BisqueRNA/main24/top300;",
  paste0("08e dtangle/", MAIN_DANGLE_MODE, "/todos os genes/100 marcadores por tipo."),
  paste0("Sensibilidade: dtangle 08e/", ALT_DTANGLE_MODE, "; NAO constitui quarto metodo."),
  "Nenhum algoritmo foi rerodado e nenhuma proporcao foi renormalizada.",
  "Contrastes nao ajustados: variacao relativa entre grupos, nao causalidade.",
  "FDR BH calculado separadamente para oito tipos por coorte e metodo.",
  "Mesmo conjunto de amostras em 3 metodos: NAO sao replicacoes independentes.",
  sprintf("HD/PD compartilham %d controles; nao sao duas replicacoes independentes.",
          length(shared_controls)),
  "MuSiC2 nao convergiu nas rodadas 06; zeros In/Opc nao sao ausencias biologicas.",
  "Bisque 07d tem ancoragem a referencia e influencia do contexto bulk.",
  "dtangle 08f tem erros elevados para alguns tipos na verdade nuclear.",
  "As metricas aqui calculadas NAO sao uma medida externa de acuracia.",
  "Nenhuma decisao ROSMAP x GSE125583 sera tomada automaticamente.",
  "A regra de agregacao usa mediana <1% em >=2 metodos, >=3 coortes, apenas sensibilidade.",
  "Somar fracoes apos o ajuste NAO soluciona falta de identificabilidade do solver.",
  "Analise inferencial final exige auditar idade/sexo/PMI/RIN/regiao/batch etc.",
  "Relatorio definitivo ainda sera reports/02_deconvolution.Rmd, depois do bMIND."
)
writeLines(notes, file.path(OUTPUT, "09_LEIA_ANTES_DE_INTERPRETAR.txt"), useBytes = TRUE)
write_tsv(rbindlist(MANIFEST, use.names = TRUE, fill = TRUE), "09_input_manifest.tsv")
write_tsv(rbindlist(STATUS, use.names = TRUE, fill = TRUE), "09_status.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUTPUT, "09_sessionInfo.txt"))

log_status("final", "tables", "PASS", "Tabelas geradas; consultar manifesto, QC e notas")
write_tsv(rbindlist(STATUS, use.names = TRUE, fill = TRUE), "09_status.tsv")
cat("\n===========================================================\n")
cat("09 CONCLUIDO. Saidas em: ", OUTPUT, "\n", sep = "")
cat("Cobertura: ", nrow(coverage), " ajustes coorte/metodo.\n", sep = "")
cat("Resultados: 09_case_control_exploratory.tsv, 09_concordancia_metodos.tsv, ",
    "09_directional_robustness_descriptive.tsv\n", sep = "")
cat("ATENCAO: efeitos exploratorios; escolha metodologica NAO automatizada.\n")
