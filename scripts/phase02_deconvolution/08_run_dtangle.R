#!/usr/bin/env Rscript
# ==============================================================================
# 08_run_dtangle.R | TCC AD–PD–HD | Fase 02: deconvolução celular
# CORREÇÃO V2 (2026-10-21): consulta de REGISTRY por índice explícito;
# evita conflito data.table entre variável local e coluna 'dataset'.
# ==============================================================================
# OBJETIVO
#   Estimar a contribuição transcriptômica relativa dos oito tipos celulares
#   definidos na referência snRNA-seq Mathys et al. (2019), pathology-free,
#   usando o método clássico dtangle (Hunt et al., Bioinformatics, 2019).
#
# COORTES (padrão: as quatro)
#   ROSMAP    | AD, principal    | 144 amostras
#   GSE125583 | AD, validação    | 289 amostras
#   GSE64810  | HD, principal    |  69 amostras
#   GSE68719  | PD, principal    |  73 amostras
#
# REFERÊNCIA E ASSINATURA
#   - 24 doadores sem patologia AD, 35.140 núcleos, 8 tipos.
#   - Os núcleos são agregados em perfis DONOR x CELL_TYPE (181 observados).
#   - Cada perfil é um pseudobulk de counts agregados do próprio doador.
#   - Normalize DEPOIS de agregar, com library size daquele pseudobulk.
#   - Os marcadores top300 vêm do script 04, NÃO são redescobertos no bulk.
#   - Importante: alguns GeneIDs aparecem em mais de um conjunto top300.
#     O artigo dtangle supõe conjuntos de marcadores disjuntos. Por isso,
#     uma atribuição por GeneID é escolhida por regra determinística:
#       menor marker_rank; empate -> maior score (quando presente);
#       novo empate -> ordem fixa Ex, In, Ast, Oli, Opc, Mic, End, Per.
#     As atribuições e as sobreposições são auditadas, sem modificar o 04.
#
# ENTRADAS E ESCALAS
#   - BULK: counts brutos top300, com library_sizes.tsv da matriz COMPLETA.
#     CPM = 1e6 * counts(top300) / library_size_full_input.
#   - REFERÊNCIA: counts por núcleo (17394 GeneIDs) são agregados por
#     doador/tipo; o CPM usa a library size de TODOS os genes mapeados.
#   - Para ambas as fontes: log2(CPM + 1). Dimensões entregues ao dtangle:
#       Y:          amostras bulk x genes
#       references: pseudobulks ref x genes
#       pure_samples: índices DAS LINHAS de references por tipo celular
#       markers:    ÍNDICES DAS COLUNAS de Y / references por tipo celular
#   - data_type="rna-seq" (gamma escolhido pelo pacote; versionado).
#   - Sem quantile normalization por padrão. Seria uma sensibilidade nova,
#     não equivalente à análise principal definida aqui.
#
# QUALIDADE / INTERPRETAÇÃO
#   - Preserva resultados brutos; não força proporções ao intervalo [0, 1].
#   - Se houver valores não finitos, negativos relevantes ou soma != 1,
#     o ajuste falha: NÃO há correção numérica silenciosa.
#   - Detecta zeros, quase zeros, vetores repetidos, variância muito baixa e
#     soluções uniformes (12,5% x oito tipos). Tais sinais geram WARN.
#   - Comparação com MuSiC2/Bisque é DESCRITIVA, não seleção automática.
#   - Sem agrupar tipos; sem inferência casos x controles; sem mudar coortes.
#   - HD e PD compartilham 44 controles: NÃO são validações independentes.
#
# EXECUÇÃO (a partir da RAIZ do repositório):
#   Rscript scripts/phase02_deconvolution/08_run_dtangle.R
#   DTANGLE_DATASET=GSE64810 Rscript scripts/phase02_deconvolution/08_run_dtangle.R
#   DTANGLE_FORCE=1 Rscript scripts/phase02_deconvolution/08_run_dtangle.R
#
# Não instala pacotes, não modifica renv.lock, não sobrescreve etapas 06/07.
# Doc. API: https://cran.r-project.org/package=dtangle
# ==============================================================================

# ──────────────────────────────────────────────────────────────────────────────
# 0) CONFIGURAÇÃO DO ESTUDO (NÃO escolher por significância dos resultados)
# ──────────────────────────────────────────────────────────────────────────────
EXPECTED_TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
DATASETS <- c("ROSMAP", "GSE125583", "GSE64810", "GSE68719")
REGISTRY <- data.table::data.table(
  dataset = c("ROSMAP", "GSE125583", "GSE64810", "GSE68719"),
  disease = c("AD", "AD", "HD", "PD"),
  cohort_role = c("principal", "validation", "principal", "principal"),
  expected_n = c(144L, 289L, 69L, 73L)
)
REF_VERSION <- "main24"
MARKER_SET <- "top300"
EXPECTED_NUCLEI <- 35140L
EXPECTED_DONORS <- 24L
EXPECTED_REFERENCE_GENES <- 17394L
EXPECTED_MARKER_UNION <- 2242L
EXPECTED_PSEUDOBULKS <- 181L
BASE_SEED <- 20261021L

# Critérios de alerta técnico, NÃO critérios de exclusão por resultado de doença.
ZERO_TOL <- 1e-10
NEAR_ZERO_TOL <- 1e-4         # fração abaixo de 0,01%
SUM_TOL <- 1e-6
NEG_TOL <- 1e-8
LOW_SD_PP <- 0.20             # alerta de SD < 0,2 ponto percentual
LOW_RANGE_PP <- 0.50          # alerta de amplitude < 0,5 pp
UNIFORM_RMSE_PP <- 2.0       # alerta de proximidade com 12,5% por tipo
MIN_MARKERS_PER_TYPE <- 15L # falha apenas se marcadores praticamente ausentes

# Variáveis de ambiente: filtros de execução; padrão processa as quatro coortes.
selection <- trimws(Sys.getenv("DTANGLE_DATASET", unset = ""))
if (nzchar(selection)) {
  selected <- trimws(strsplit(selection, ",", fixed = TRUE)[[1]])
  if (!length(selected) || any(!selected %in% DATASETS)) {
    stop("DTANGLE_DATASET inválido. Use ROSMAP,GSE125583,GSE64810,GSE68719.")
  }
  DATASETS <- unique(selected)
}
FORCE_RECOMPUTE <- identical(Sys.getenv("DTANGLE_FORCE", "0"), "1")

BULK_ROOT <- file.path("data", "processed", "deconvolution")
REF_ROOT <- file.path("data", "reference", "mathys2019", "phase02_prepared", "reference")
MARKER_ROOT <- file.path("data", "reference", "mathys2019", "phase02_prepared", "markers")
OUT_ROOT <- file.path("results", "deconvolution", "dtangle")
REF_COUNTS_FILE <- file.path(REF_ROOT, "mathys2019_pathology_free_geneid_counts_sparse.rds")
REF_META_FILE <- file.path(REF_ROOT, "mathys2019_pathology_free_cell_metadata_aligned.tsv")
MARKER_FILE <- file.path(MARKER_ROOT, "mathys2019_top300_markers.tsv")
STEP05_STATUS <- file.path("results", "deconvolution", "bulk_preparation", "05_status.tsv")
MUSIC2_LONG <- file.path("results", "deconvolution", "MuSiC2", "06_proportions_long_all.tsv")
BISQUE_LONG <- file.path("results", "deconvolution", "BisqueRNA", "07_proportions_long_all.tsv")

CELL_COLORS <- c(Ex="#3B5BDB", In="#7950F2", Ast="#0CA678",
                 Oli="#F08C00", Opc="#FCC419", Mic="#C2255C",
                 End="#15AABF", Per="#868E96")
GROUP_COLORS <- c(Control="#94A3B8", AD="#7B2CBF",
                  HD="#277DA1", PD="#E75480")

# ──────────────────────────────────────────────────────────────────────────────
# 1) DEPENDÊNCIAS E FUNÇÕES UTILITÁRIAS
# ──────────────────────────────────────────────────────────────────────────────
REQUIRED <- c("data.table", "Matrix", "ggplot2", "dtangle")
missing <- REQUIRED[!vapply(REQUIRED, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop(
  "Pacotes ausentes: ", paste(missing, collapse = ", "),
  ". Instale no renv do projeto (sem snapshot automático)."
)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(sprintf(...), "\n", sep = "")
fail <- function(...) stop(paste0(...), call. = FALSE)
clean <- function(x) trimws(as.character(x))
write_tab <- function(dt, filename) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(dt, filename, sep = "\t", quote = FALSE, na = "NA")
  invisible(filename)
}

check_ids <- function(x, label) {
  x <- clean(x)
  if (!length(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) {
    fail(label, ": IDs vazios, duplicados ou ausentes.")
  }
  invisible(TRUE)
}
check_counts <- function(x, label) {
  if (is.null(dim(x)) || is.null(rownames(x)) || is.null(colnames(x))) {
    fail(label, ": matriz sem nomes de genes ou amostras.")
  }
  check_ids(rownames(x), paste0(label, " genes"))
  check_ids(colnames(x), paste0(label, " amostras"))
  v <- if (inherits(x, "sparseMatrix")) x@x else as.numeric(x)
  if (any(!is.finite(v)) || any(v < 0) || any(abs(v - round(v)) > 1e-7)) {
    fail(label, ": counts devem ser inteiros, finitos e não negativos.")
  }
  invisible(TRUE)
}
file_stamp <- function(path) {
  i <- file.info(path)
  if (nrow(i) != 1 || is.na(i$size)) fail("Arquivo inacessível: ", path)
  paste(normalizePath(path, winslash = "/", mustWork = TRUE),
        as.numeric(i$size), as.numeric(i$mtime), sep = "|")
}
run_status <- function(ds, item, status, observed, note = "") {
  data.table(dataset = ds, reference_version = REF_VERSION,
             marker_set = MARKER_SET, check = item, status = status,
             observed = as.character(observed), note = as.character(note))
}
# Coeficiente não interpretável quando uma entrada é constante.
safe_spearman <- function(x, y) {
  if (length(x) < 3L || sd(x) < 1e-12 || sd(y) < 1e-12) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

# Conferência automática dos argumentos públicos: evita usar uma API divergente.
api_names <- names(formals(dtangle::dtangle))
api_need <- c("Y", "references", "pure_samples", "markers", "data_type")
if (!all(api_need %in% api_names)) {
  fail("A API da versão instalada de dtangle difere da esperada: ",
       paste(setdiff(api_need, api_names), collapse = ", "))
}
api_tab <- data.table(package = "dtangle",
                      version = as.character(utils::packageVersion("dtangle")),
                      public_arguments = paste(api_names, collapse = ","),
                      expected_arguments_present = TRUE,
                      input_scale = "log2(CPM+1)", data_type = "rna-seq",
                      markers = "top300 disjuntos por rank/score",
                      reference = "181 pseudobulks donor x cell_type")
write_tab(api_tab, file.path(OUT_ROOT, "08_api_audit.tsv"))

# ──────────────────────────────────────────────────────────────────────────────
# 2) VALIDAÇÃO DO ESTUDO E REFERÊNCIA Mathys
# ──────────────────────────────────────────────────────────────────────────────
needed <- c(REF_COUNTS_FILE, REF_META_FILE, MARKER_FILE, STEP05_STATUS)
if (any(!file.exists(needed))) {
  fail("Entradas ausentes. Execute etapas 03-05: ",
       paste(needed[!file.exists(needed)], collapse = "; "))
}
status05 <- fread(STEP05_STATUS)
if (!"status" %in% names(status05)) fail("05_status.tsv sem coluna status.")
if (any(status05$status == "FAIL", na.rm = TRUE)) {
  fail("Etapa 05 possui falhas (05_status.tsv). Corrija antes do 08.")
}

say("========================================================")
say("08 | dtangle %s | Mathys pathology-free / top300", api_tab$version[[1]])
say("Coortes nesta chamada: %s", paste(DATASETS, collapse = ", "))
say("========================================================")

sc <- readRDS(REF_COUNTS_FILE)
if (!inherits(sc, "dgCMatrix")) sc <- as(sc, "dgCMatrix")
check_counts(sc, "Referência Mathys")
if (nrow(sc) != EXPECTED_REFERENCE_GENES || ncol(sc) != EXPECTED_NUCLEI) {
  fail("Dimensões inesperadas da referência: ", nrow(sc), " x ", ncol(sc))
}
ref_meta <- fread(REF_META_FILE)
for (f in c("TAG", "donor_id", "cell_type")) {
  if (!f %in% names(ref_meta)) fail("Metadata Mathys sem coluna ", f)
}
ref_meta[, `:=`(TAG = clean(TAG), donor_id = clean(donor_id),
                cell_type = clean(cell_type))]
check_ids(ref_meta$TAG, "Metadata TAG")
idx <- match(colnames(sc), ref_meta$TAG)
if (anyNA(idx)) fail("Núcleos da matriz não correspondem aos metadados.")
ref_meta <- ref_meta[idx]
if (!identical(ref_meta$TAG, colnames(sc))) fail("Ordenação incorreta da referência.")
if (uniqueN(ref_meta$donor_id) != EXPECTED_DONORS ||
    !setequal(unique(ref_meta$cell_type), EXPECTED_TYPES)) {
  fail("A referência deve ter 24 doadores e os oito tipos previstos.")
}

# 2A) Um pseudobulk por combinação OBSERVADA doador-tipo, sem imputar ausências.
ref_groups <- unique(ref_meta[, .(donor_id, cell_type)])
ref_groups[, ct_order := match(cell_type, EXPECTED_TYPES)]
setorder(ref_groups, ct_order, donor_id)
ref_groups[, ct_order := NULL]
if (nrow(ref_groups) != EXPECTED_PSEUDOBULKS) {
  fail("Esperados 181 pseudobulks observados; encontrados ", nrow(ref_groups))
}
ref_groups[, ref_id := sprintf("REF_%03d", .I)]
ref_key <- paste(ref_groups$donor_id, ref_groups$cell_type, sep = "___")
assigned <- match(paste(ref_meta$donor_id, ref_meta$cell_type, sep = "___"), ref_key)
if (anyNA(assigned)) fail("Falha na atribuição núcleo -> pseudobulk.")
ref_groups[, n_nuclei := tabulate(assigned, nbins = .N)]

# Produto de matrizes esparsas: genes x núcleos MULTIPLICADO por núcleos x grupos.
# Não densificar 35.140 núcleos. Usa os 17.394 genes para biblioteca completa.
assignment_matrix <- Matrix::sparseMatrix(
  i = seq_along(assigned), j = assigned, x = rep(1, length(assigned)),
  dims = c(ncol(sc), nrow(ref_groups))
)
ref_pb <- sc %*% assignment_matrix
rownames(ref_pb) <- rownames(sc)
colnames(ref_pb) <- ref_groups$ref_id
ref_lib_full <- as.numeric(Matrix::colSums(ref_pb))
if (any(ref_lib_full <= 0) || any(!is.finite(ref_lib_full))) {
  fail("Há pseudobulks de referência com biblioteca inválida.")
}
ref_groups[, library_size_mapped_allgenes := ref_lib_full]
ref_groups[, library_size_per_nucleus := ref_lib_full / n_nuclei]
write_tab(ref_groups, file.path(OUT_ROOT, "08_reference_pseudobulk_audit.tsv"))
ref_support <- ref_groups[, .(n_donors = .N,
                              n_nuclei = sum(n_nuclei),
                              min_nuclei_per_donor = min(n_nuclei),
                              median_nuclei_per_donor = median(n_nuclei)),
                          by = cell_type]
ref_support[, cell_type := factor(cell_type, levels = EXPECTED_TYPES)]
setorder(ref_support, cell_type)
ref_support[, cell_type := as.character(cell_type)]
write_tab(ref_support, file.path(OUT_ROOT, "08_reference_support.tsv"))

# A ordem de references e pure_samples é explícita e fixa.
ref_pure_samples <- setNames(lapply(EXPECTED_TYPES, function(ct) {
  which(ref_groups$cell_type == ct)
}), EXPECTED_TYPES)
if (any(lengths(ref_pure_samples) < 2L)) fail("Tipo sem doadores suficientes.")
if (!identical(sort(unlist(ref_pure_samples, use.names = FALSE)),
               seq_len(nrow(ref_groups)))) {
  fail("Os índices pure_samples não cobrem exatamente os pseudobulks.")
}

# ──────────────────────────────────────────────────────────────────────────────
# 3) MARCADORES TOP300: RESOLVER SOBREPOSIÇÕES ENTRE TIPOS
# ──────────────────────────────────────────────────────────────────────────────
markers_raw <- fread(MARKER_FILE)
for (f in c("gene_id", "cell_type", "marker_rank")) {
  if (!f %in% names(markers_raw)) fail("Arquivo top300 sem coluna ", f)
}
markers_raw[, `:=`(gene_id = clean(gene_id), cell_type = clean(cell_type),
                   marker_rank = as.integer(marker_rank))]
if (anyNA(markers_raw$marker_rank) || any(markers_raw$marker_rank < 1) ||
    any(markers_raw$marker_rank > 300) ||
    !setequal(unique(markers_raw$cell_type), EXPECTED_TYPES)) {
  fail("Ranks ou tipos celulares inválidos nos marcadores do 04.")
}
if (anyDuplicated(markers_raw[, .(gene_id, cell_type)])) {
  fail("Marcadores duplicados no mesmo tipo celular.")
}
marker_sizes <- markers_raw[, .N, by = cell_type]
if (any(marker_sizes$N != 300L) || uniqueN(markers_raw$gene_id) != EXPECTED_MARKER_UNION) {
  fail("Top300 difere do conjunto do script 04 (300/tipo, união 2242).")
}

# O escore auxiliar só desempata ranks iguais, quando disponível. Se ausente,
# a decisão é feita por rank e pela ordem fixa dos oito tipos celulares.
possible_scores <- c("score", "marker_score", "ranking_score", "combined_score")
score_column <- intersect(possible_scores, names(markers_raw))
score_column <- if (length(score_column)) score_column[[1]] else NA_character_
markers_raw[, score_for_tie := if (!is.na(score_column))
              suppressWarnings(as.numeric(get(score_column))) else NA_real_]
markers_raw[!is.finite(score_for_tie), score_for_tie := -Inf]
markers_raw[, ct_order := match(cell_type, EXPECTED_TYPES)]
setorderv(markers_raw,
         cols = c("gene_id", "marker_rank", "score_for_tie", "ct_order"),
         order = c(1L, 1L, -1L, 1L))
markers_raw[, choice := !duplicated(gene_id)]
markers_raw[, assignments_for_gene := .N, by = gene_id]
markers_raw[, overlap := assignments_for_gene > 1L]
markers_primary <- copy(markers_raw[choice == TRUE])
if (anyDuplicated(markers_primary$gene_id)) fail("Falha ao remover conflitos top300.")
write_tab(markers_raw[, .(gene_id, cell_type, marker_rank, score_for_tie,
                          overlap, assignments_for_gene, choice)],
          file.path(OUT_ROOT, "08_marker_overlap_resolution.tsv"))
write_tab(markers_primary[, .(gene_id, cell_type, marker_rank, score_for_tie)],
          file.path(OUT_ROOT, "08_markers_unique_assignment.tsv"))
say("Top300: %d atribuições, %d GeneIDs únicos, %d genes com conflito.",
    nrow(markers_raw), nrow(markers_primary),
    uniqueN(markers_raw[overlap == TRUE, gene_id]))

# Evitar manter a matriz de TODOS os núcleos na memória após construir pseudobulk.
rm(sc, assignment_matrix)
invisible(gc(verbose = FALSE))

# ──────────────────────────────────────────────────────────────────────────────
# 4) PREPARAÇÃO POR COORTE: GENES, METADADOS, CPM, LOG2
# ──────────────────────────────────────────────────────────────────────────────
read_inputs <- function(dataset) {
  # NÃO usar REGISTRY[REGISTRY$dataset == dataset]: dentro de data.table,
  # o nome 'dataset' pode resolver para a coluna da própria tabela.
  # match() evita avaliação não padrão e seleciona UMA única coorte.
  dataset_key <- as.character(dataset)
  if (length(dataset_key) != 1L || is.na(dataset_key) ||
      !nzchar(dataset_key)) fail("Identificador de coorte inválido.")
  registry_index <- match(dataset_key, REGISTRY[["dataset"]])
  if (is.na(registry_index)) fail("Coorte desconhecida: ", dataset_key)
  reg <- REGISTRY[registry_index]
  if (nrow(reg) != 1L ||
      !identical(as.character(reg[["dataset"]][[1L]]), dataset_key)) {
    fail("Falha de integridade ao consultar registro da coorte: ", dataset_key)
  }
  disease_label <- reg$disease[[1]]
  bulk_file <- file.path(BULK_ROOT, dataset, "bulk_counts_top300.rds")
  meta_file <- file.path(BULK_ROOT, dataset, "metadata_deconvolution.tsv")
  library_file <- file.path(BULK_ROOT, dataset, "library_sizes.tsv")
  paths <- c(bulk_file, meta_file, library_file)
  if (any(!file.exists(paths))) {
    fail(dataset, ": arquivo(s) ausente(s): ",
         paste(paths[!file.exists(paths)], collapse = ", "))
  }
  bulk <- as.matrix(readRDS(bulk_file))
  storage.mode(bulk) <- "double"
  check_counts(bulk, paste0(dataset, " bulk"))
  if (ncol(bulk) != reg$expected_n[[1]]) {
    fail(dataset, ": n amostras inesperado ", ncol(bulk),
         "; esperado ", reg$expected_n[[1]])
  }
  md <- fread(meta_file)
  req_md <- c("sample_id", "group", "disease", "cohort_role")
  if (length(setdiff(req_md, names(md)))) fail(dataset, ": metadados incompletos.")
  md[, sample_id := clean(sample_id)]
  check_ids(md$sample_id, paste0(dataset, " metadata"))
  if (!setequal(md$sample_id, colnames(bulk))) {
    fail(dataset, ": bulk/metadata não possuem as mesmas amostras.")
  }
  md <- md[match(colnames(bulk), md$sample_id)]
  if (!identical(as.character(md$sample_id), as.character(colnames(bulk)))) {
    fail(dataset, ": ordem bulk/metadata incompatível.")
  }
  if (!setequal(unique(md$group), c("Control", disease_label)) ||
      any(md$disease != disease_label) ||
      any(md$cohort_role != reg$cohort_role[[1]])) {
    fail(dataset, ": grupos/atributos diferentes do registro.")
  }
  libs <- fread(library_file)
  if (!all(c("sample_id", "library_size_full_input") %in% names(libs))) {
    fail(dataset, ": library_sizes.tsv inválido.")
  }
  libs[, sample_id := clean(sample_id)]
  check_ids(libs$sample_id, paste0(dataset, " library sample_ids"))
  lib <- as.numeric(libs$library_size_full_input[match(md$sample_id, libs$sample_id)])
  if (any(!is.finite(lib)) || any(lib <= 0)) fail(dataset, ": library sizes inválidas.")
  observed_lib <- colSums(bulk)
  if (any(lib + 1e-6 < observed_lib)) {
    fail(dataset, ": library_size_full_input é inferior à soma top300.")
  }

  # Alinhar GeneIDs SEM imputar ausências; manter mesma ordem nas duas fontes.
  raw_genes <- rownames(bulk)
  common <- raw_genes[raw_genes %in% rownames(ref_pb)]
  if (length(common) < 500L) fail(dataset, ": menos de 500 genes em comum.")
  bulk <- bulk[common, , drop = FALSE]
  ref_sub <- ref_pb[common, , drop = FALSE]
  keep <- rowSums(bulk) > 0 & Matrix::rowSums(ref_sub) > 0
  genes <- common[keep]
  if (length(genes) < 500L) fail(dataset, ": poucos genes após excluir linhas sem expressão.")
  bulk <- bulk[genes, , drop = FALSE]
  ref_sub <- ref_sub[genes, , drop = FALSE]
  if (!identical(rownames(bulk), rownames(ref_sub))) fail("GeneIDs fora de ordem.")

  # CPM usa denominadores COMPLETOS, jamais apenas a soma top300.
  bulk_cpm <- sweep(bulk, 2L, lib, "/") * 1e6
  ref_cpm <- sweep(as.matrix(ref_sub), 2L, ref_lib_full, "/") * 1e6
  y_bulk <- t(log2(bulk_cpm + 1))
  y_ref <- t(log2(ref_cpm + 1))
  rownames(y_ref) <- ref_groups$ref_id
  if (!identical(colnames(y_bulk), colnames(y_ref))) {
    fail("GeneIDs entre bulk e referência não estão alinhados nas matrizes finais.")
  }
  if (any(!is.finite(y_bulk)) || any(!is.finite(y_ref))) {
    fail(dataset, ": log2CPM contém valores não finitos.")
  }

  # Marcadores do 04, únicos em um só tipo celular, presentes na coorte.
  selected <- copy(markers_primary[gene_id %in% genes])
  selected[, marker_column_index := match(gene_id, genes)]
  if (anyNA(selected$marker_column_index)) fail("Não foi possível indexar os marcadores.")
  marker_list <- setNames(lapply(EXPECTED_TYPES, function(ct) {
    as.integer(selected[cell_type == ct, marker_column_index])
  }), EXPECTED_TYPES)
  if (any(lengths(marker_list) < MIN_MARKERS_PER_TYPE)) {
    fail(dataset, ": tipo com menos de ", MIN_MARKERS_PER_TYPE,
         " marcadores exclusivos disponíveis: ",
         paste(names(marker_list), lengths(marker_list), sep = "=", collapse = "; "))
  }
  if (anyDuplicated(unlist(marker_list, use.names = FALSE))) {
    fail(dataset, ": mesma coluna foi atribuída a mais de um tipo.")
  }
  list(dataset = dataset, disease = disease_label,
       cohort_role = reg$cohort_role[[1]],
       bulk_file = bulk_file, meta_file = meta_file,
       library_file = library_file, metadata = md,
       lib = lib, genes = genes, common = common,
       removed = common[!keep], bulk_log = y_bulk,
       ref_log = y_ref, markers = marker_list, marker_table = selected,
       bulk_count_total = observed_lib)
}

# ──────────────────────────────────────────────────────────────────────────────
# 5) ESTIMAÇÃO DO DTANGLE E INTEGRIDADE (SEM RENORMALIZAR RESULTADOS)
# ──────────────────────────────────────────────────────────────────────────────
fit_signature <- function(x) {
  list(
    dataset = x$dataset,
    sample_ids = as.character(x$metadata$sample_id),
    genes = x$genes,
    marker_assignment = paste(x$marker_table$gene_id, x$marker_table$cell_type,
                              x$marker_table$marker_column_index, sep = ":"),
    ref_donor_n = EXPECTED_DONORS,
    ref_pb_n = nrow(ref_groups),
    input_scale = "log2(CPM+1)",
    dtangle_version = as.character(utils::packageVersion("dtangle")),
    source_stamp = vapply(c(REF_COUNTS_FILE, REF_META_FILE, MARKER_FILE,
                             x$bulk_file, x$meta_file, x$library_file),
                          file_stamp, character(1))
  )
}
run_dtangle <- function(x, run_dir) {
  cp_file <- file.path(run_dir, "08_fit_checkpoint.rds")
  signature <- fit_signature(x)
  if (!FORCE_RECOMPUTE && file.exists(cp_file)) {
    old <- tryCatch(readRDS(cp_file), error = function(e) NULL)
    if (is.list(old) && identical(old$signature, signature) &&
        is.list(old$result) && !is.null(old$result$estimates)) {
      say("  Checkpoint compatível reutilizado.")
      return(list(result = old$result, checkpoint_reused = TRUE,
                  elapsed_seconds = old$elapsed_seconds))
    }
    say("  Checkpoint incompatível/incompleto: recalculando.")
  }
  set.seed(BASE_SEED)
  tic <- proc.time()[["elapsed"]]
  result <- dtangle::dtangle(
    Y = x$bulk_log,
    references = x$ref_log,
    pure_samples = ref_pure_samples,
    markers = x$markers,
    data_type = "rna-seq"
  )
  seconds <- unname(proc.time()[["elapsed"]] - tic)
  # Checkpoint SEMPRE antes do pós-processamento/gráficos.
  saveRDS(list(signature = signature, result = result,
               elapsed_seconds = seconds), cp_file, compress = TRUE)
  list(result = result, checkpoint_reused = FALSE,
       elapsed_seconds = seconds)
}

validate_result <- function(x, result) {
  e <- as.matrix(result$estimates)
  if (!is.numeric(e) || nrow(e) != nrow(x$metadata) ||
      ncol(e) != length(EXPECTED_TYPES)) {
    fail(x$dataset, ": dimensões do dtangle inesperadas.")
  }
  # A documentação garante ordem pelos índices de pure_samples; ainda assim,
  # conferir quaisquer nomes antes de atribuir nomes fixos às colunas.
  if (!is.null(colnames(e)) &&
      !identical(as.character(colnames(e)), EXPECTED_TYPES)) {
    if (setequal(colnames(e), EXPECTED_TYPES)) {
      e <- e[, EXPECTED_TYPES, drop = FALSE]
    } else {
      fail(x$dataset, ": rótulos dos oito tipos divergentes: ",
           paste(colnames(e), collapse = ","))
    }
  }
  if (is.null(colnames(e))) colnames(e) <- EXPECTED_TYPES
  if (!is.null(rownames(e)) &&
      !identical(as.character(rownames(e)), as.character(x$metadata$sample_id))) {
    if (setequal(rownames(e), x$metadata$sample_id)) {
      e <- e[x$metadata$sample_id, , drop = FALSE]
    } else {
      fail(x$dataset, ": amostras devolvidas pelo dtangle não correspondem ao bulk.")
    }
  }
  rownames(e) <- x$metadata$sample_id
  if (any(!is.finite(e))) fail(x$dataset, ": dtangle retornou NA/Inf.")
  if (any(e < -NEG_TOL)) {
    fail(x$dataset, ": dtangle retornou proporções negativas relevantes: ", min(e))
  }
  if (any(e > 1 + NEG_TOL)) {
    fail(x$dataset, ": dtangle retornou proporções maiores que 1.")
  }
  if (max(abs(rowSums(e) - 1)) > SUM_TOL) {
    fail(x$dataset, ": dtangle não retornou linhas somando 1. Maior erro: ",
         max(abs(rowSums(e) - 1)))
  }
  e
}

# ──────────────────────────────────────────────────────────────────────────────
# 6) TABELAS DE RESULTADOS E ALERTAS POR COORTE
# ──────────────────────────────────────────────────────────────────────────────
make_outputs <- function(x, e) {
  sample_ids <- x$metadata$sample_id
  wide <- data.table(sample_id = sample_ids)
  for (ct in EXPECTED_TYPES) wide[[ct]] <- as.numeric(e[, ct])
  wide[, `:=`(dataset = x$dataset, disease = x$disease,
             cohort_role = x$cohort_role, reference_version = REF_VERSION,
             marker_set = MARKER_SET)]
  setcolorder(wide, c("dataset", "disease", "cohort_role", "reference_version",
                      "marker_set", "sample_id", EXPECTED_TYPES))

  long <- melt(copy(wide), id.vars = setdiff(names(wide), EXPECTED_TYPES),
               measure.vars = EXPECTED_TYPES,
               variable.name = "cell_type", value.name = "proportion")
  long[, cell_type := as.character(cell_type)]
  gr <- x$metadata[, .(sample_id, group)]
  long <- merge(long, gr, by = "sample_id", all.x = TRUE, sort = FALSE)
  setcolorder(long, c("dataset", "disease", "cohort_role", "reference_version",
                      "marker_set", "sample_id", "group", "cell_type", "proportion"))
  if (anyNA(long$group)) fail("Metadados de grupo incompletos após merge.")
  if (nrow(long) != nrow(wide) * length(EXPECTED_TYPES)) {
    fail("long output: número de linhas inesperado.")
  }

  grp <- long[, .(n_samples = .N,
                  mean_proportion = mean(proportion),
                  median_proportion = median(proportion),
                  sd_pp = 100 * sd(proportion),
                  min_proportion = min(proportion),
                  q1_proportion = as.numeric(quantile(proportion, 0.25)),
                  q3_proportion = as.numeric(quantile(proportion, 0.75)),
                  max_proportion = max(proportion),
                  fraction_zero = mean(proportion <= ZERO_TOL),
                  fraction_near_zero = mean(proportion < NEAR_ZERO_TOL)),
              by = .(dataset, disease, cohort_role, reference_version,
                     marker_set, group, cell_type)]
  variation <- long[, .(n_samples = .N,
                        mean_proportion = mean(proportion),
                        median_proportion = median(proportion),
                        sd_pp = 100 * sd(proportion),
                        range_pp = 100 * diff(range(proportion)),
                        fraction_zero = mean(proportion <= ZERO_TOL),
                        fraction_near_zero = mean(proportion < NEAR_ZERO_TOL),
                        n_distinct_6dec = uniqueN(round(proportion, 6))),
                    by = .(dataset, reference_version, marker_set, cell_type)]
  variation[, `:=`(low_sd = sd_pp < LOW_SD_PP,
                   low_range = range_pp < LOW_RANGE_PP)]

  uniform <- rep(1 / length(EXPECTED_TYPES), length(EXPECTED_TYPES))
  diag <- data.table(dataset = x$dataset,
                     reference_version = REF_VERSION, marker_set = MARKER_SET,
                     sample_id = sample_ids,
                     group = x$metadata$group,
                     sum_proportion = rowSums(e),
                     n_exact_zero = rowSums(e <= ZERO_TOL),
                     n_near_zero = rowSums(e < NEAR_ZERO_TOL),
                     rmse_to_uniform_pp = 100 * sqrt(rowMeans(
                       (e - matrix(uniform, nrow = nrow(e), ncol = ncol(e),
                                   byrow = TRUE))^2)))
  median_rmse <- median(diag$rmse_to_uniform_pp)
  vector_string <- apply(round(e, 6), 1L, paste, collapse = "|")
  unique_vectors <- uniqueN(vector_string)
  degenerate_uniform <- median_rmse < UNIFORM_RMSE_PP &&
    all(variation$low_sd)
  degenerate_constant <- unique_vectors <= 1L ||
    all(variation$low_sd & variation$low_range)
  flags <- data.table(dataset = x$dataset,
                      reference_version = REF_VERSION, marker_set = MARKER_SET,
                      uniform_rmse_median_pp = median_rmse,
                      unique_vectors_6dec = unique_vectors,
                      degenerate_uniform = degenerate_uniform,
                      degenerate_constant = degenerate_constant)
  list(wide = wide, long = long, group = grp, variation = variation,
       sample_diag = diag, flags = flags)
}

# ──────────────────────────────────────────────────────────────────────────────
# 7) VISUALIZAÇÕES: DIMENSÕES AMPLAS, EIXOS NÃO CORTADOS
# ──────────────────────────────────────────────────────────────────────────────
plot_theme <- function(size = 12) {
  theme_minimal(base_size = size) +
    theme(plot.title = element_text(face = "bold", size = 17),
          plot.subtitle = element_text(size = 10.5),
          plot.margin = margin(18, 26, 36, 30),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          axis.text.x = element_text(color = "#3F444C"),
          legend.position = "bottom")
}
save_plot <- function(p, path, width, height, dpi = 185) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = p, width = width, height = height,
                  dpi = dpi, units = "in", limitsize = FALSE, bg = "white")
  data.table(file = path, status = "PASS", error = NA_character_)
}
make_plots <- function(out, x, plot_dir) {
  disease_label <- x$disease   # evita conflito com coluna 'disease' em data.table
  d <- copy(out$long)
  d[, group := factor(group, levels = c("Control", disease_label))]
  d[, cell_type := factor(cell_type, levels = EXPECTED_TYPES)]
  d[, proportion_pp := 100 * proportion]
  attempts <- list()
  safe_plot <- function(key, expr) {
    # Somente a figura pode falhar; resultados numéricos continuam disponíveis.
    attempts[[key]] <<- tryCatch(expr, error = function(e)
      data.table(file = key, status = "WARN", error = conditionMessage(e)))
  }

  # Facetas com eixo Y independente; expansão impede clipping de pontos extremos.
  p1 <- ggplot(d, aes(x = group, y = proportion_pp, fill = group)) +
    geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.77,
                 linewidth = 0.40, color = "#313842") +
    geom_jitter(aes(color = group), width = 0.14, height = 0,
                size = 0.84, alpha = 0.42, show.legend = FALSE) +
    facet_wrap(~ cell_type, ncol = 4, scales = "free_y") +
    scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
    scale_color_manual(values = GROUP_COLORS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.18))) +
    coord_cartesian(clip = "off") +
    labs(title = paste0("dtangle • ", x$dataset, " • proporções individuais"),
         subtitle = "Diagnóstico exploratório; escala Y independente por tipo celular",
         x = "Grupo", y = "Contribuição estimada (%)", fill = "Grupo") +
    plot_theme(12) +
    theme(axis.text.x = element_text(angle = 32, hjust = 1, vjust = 1))
  safe_plot("boxplots", save_plot(p1, file.path(plot_dir,
                                  "08_boxplots_by_celltype.png"), 17.2, 10.9))

  mean_dt <- d[, .(mean_pp = mean(proportion_pp)), by = .(group, cell_type)]
  # Checar antes do empilhamento para não gerar gráfico de soma truncada.
  totals <- mean_dt[, .(sum_of_type_means_pp = sum(mean_pp)), by = group]
  if (any(abs(totals$sum_of_type_means_pp - 100) > 1e-4)) {
    fail(x$dataset, ": médias por grupo não somam 100%.")
  }
  write_tab(totals, file.path(dirname(plot_dir), "08_plot_group_totals.tsv"))
  p2 <- ggplot(mean_dt, aes(x = group, y = mean_pp, fill = cell_type)) +
    geom_col(width = 0.64, color = "white", linewidth = 0.20) +
    scale_fill_manual(values = CELL_COLORS, breaks = EXPECTED_TYPES,
                      drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0.025))) +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    labs(title = paste0("dtangle • ", x$dataset, " • composição média"),
         subtitle = "Médias de contribuições relativas; não equivalem a abundância absoluta",
         x = "Grupo", y = "Contribuição média estimada (%)", fill = "Tipo celular") +
    plot_theme(12)
  safe_plot("composicao", save_plot(p2, file.path(plot_dir,
                                      "08_mean_composition_stacked.png"), 10.5, 6.4))

  v <- copy(out$variation)
  v[, cell_type := factor(cell_type, levels = EXPECTED_TYPES)]
  p3 <- ggplot(v, aes(x = cell_type, y = sd_pp, fill = cell_type)) +
    geom_col(width = 0.73, show.legend = FALSE) +
    geom_hline(yintercept = LOW_SD_PP, color = "#B42318", linetype = "dashed") +
    scale_fill_manual(values = CELL_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    labs(title = paste0("dtangle • ", x$dataset, " • variabilidade entre amostras"),
         subtitle = "Linha tracejada: alerta operacional de SD < 0,20 pp",
         x = "Tipo celular", y = "Desvio-padrão (pontos percentuais)") +
    plot_theme(12)
  safe_plot("variabilidade", save_plot(p3, file.path(plot_dir,
                                      "08_variability_sd.png"), 10.5, 5.8))

  p4 <- ggplot(v, aes(x = cell_type, y = 100 * fraction_zero,
                      fill = cell_type)) +
    geom_col(width = 0.73, show.legend = FALSE) +
    scale_fill_manual(values = CELL_COLORS) +
    scale_y_continuous(breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0.01, 0.04))) +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    labs(title = paste0("dtangle • ", x$dataset, " • estimativas iguais a zero"),
         subtitle = "Zero estimado não demonstra ausência celular biológica",
         x = "Tipo celular", y = "Amostras com estimativa zero (%)") +
    plot_theme(12)
  safe_plot("zeros", save_plot(p4, file.path(plot_dir,
                                 "08_zero_fraction.png"), 10.5, 5.8))
  rbindlist(attempts, fill = TRUE)
}

# ──────────────────────────────────────────────────────────────────────────────
# 8) EXECUTAR POR COORTE; PRESERVAR AS SAÍDAS MESMO SE FIGURA FALHAR
# ──────────────────────────────────────────────────────────────────────────────
# Preflight explícito: detectar IDs ausentes ou repetidos ANTES de executar ajustes.
if (anyDuplicated(REGISTRY[["dataset"]])) {
  fail("REGISTRY possui nomes de coorte duplicados.")
}
selected_registry_indices <- match(DATASETS, REGISTRY[["dataset"]])
if (anyNA(selected_registry_indices)) {
  fail("Coortes não encontradas no REGISTRY: ",
       paste(DATASETS[is.na(selected_registry_indices)], collapse = ", "))
}
plan <- copy(REGISTRY[selected_registry_indices])
if (!identical(as.character(plan[["dataset"]]), as.character(DATASETS))) {
  fail("Ordem das coortes no registro não corresponde à seleção.")
}
plan[, `:=`(reference_version = REF_VERSION, marker_set = MARKER_SET)]
write_tab(plan, file.path(OUT_ROOT, "08_run_plan.tsv"))

all_results <- list()
all_status <- list()
all_plot_status <- list()
run_errors <- character()
for (ds in DATASETS) {
  cat("\n------------------------------------------------------------\n")
  say("dtangle: %s | %s | %s", ds, MARKER_SET, REF_VERSION)
  run_dir <- file.path(OUT_ROOT, ds, REF_VERSION, MARKER_SET)
  plot_dir <- file.path(run_dir, "plots")
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  one <- tryCatch({
    input <- read_inputs(ds)
    input_audit <- data.table(
      dataset = ds, reference_version = REF_VERSION, marker_set = MARKER_SET,
      bulk_samples = nrow(input$bulk_log),
      reference_nuclei = EXPECTED_NUCLEI,
      reference_donors = EXPECTED_DONORS,
      pseudobulks = nrow(input$ref_log),
      genes_common_before_filter = length(input$common),
      genes_used = length(input$genes),
      genes_removed_all_zero = length(input$removed),
      marker_assignments_used = sum(lengths(input$markers)),
      shared_bulk_full_libraries = TRUE,
      bulk_lib_min = min(input$lib),
      bulk_lib_median = median(input$lib),
      bulk_lib_max = max(input$lib)
    )
    write_tab(input_audit, file.path(run_dir, "input_audit.tsv"))
    write_tab(input$marker_table[, .(gene_id, cell_type, marker_rank,
                                     score_for_tie, marker_column_index)],
              file.path(run_dir, "markers_used.tsv"))
    marker_counts <- data.table(cell_type = EXPECTED_TYPES,
                                n_markers = as.integer(lengths(input$markers)))
    write_tab(marker_counts, file.path(run_dir, "marker_count_by_type.tsv"))
    write_tab(data.table(gene_id = input$genes),
              file.path(run_dir, "genes_used.tsv"))
    write_tab(data.table(gene_id = input$removed),
              file.path(run_dir, "genes_removed_allzero.tsv"))
    pars <- data.table(parameter = c("package", "version", "dataset",
                                     "reference_version", "marker_set",
                                     "expression_input", "normalization",
                                     "library_size_bulk", "library_size_ref",
                                     "reference_unit", "marker_overlap_rule",
                                     "data_type", "force_recompute", "seed"),
                       value = c("dtangle", api_tab$version[[1]], ds,
                                 REF_VERSION, MARKER_SET, "raw counts",
                                 "log2(CPM+1)", "full input library sizes",
                                 "all mapped reference genes (17394)",
                                 "one observed donor x cell_type pseudobulk",
                                 "min rank; tie max score; tie fixed type order",
                                 "rna-seq", as.character(FORCE_RECOMPUTE),
                                 as.character(BASE_SEED)))
    write_tab(pars, file.path(run_dir, "parameters.tsv"))
    say("Genes: %d | bulk: %d | pseudobulks: %d | markers únicos: %d",
        length(input$genes), nrow(input$bulk_log), nrow(input$ref_log),
        sum(lengths(input$markers)))

    fit <- run_dtangle(input, run_dir)
    # Preserva saída bruta separadamente do checkpoint para inspeção.
    saveRDS(fit$result, file.path(run_dir, "dtangle_result.rds"), compress = TRUE)
    est <- validate_result(input, fit$result)
    out <- make_outputs(input, est)
    write_tab(out$wide, file.path(run_dir, "proportions_wide.tsv"))
    write_tab(out$long, file.path(run_dir, "proportions_long.tsv"))
    write_tab(out$group, file.path(run_dir, "proportion_summary_by_group.tsv"))
    write_tab(out$variation, file.path(run_dir, "celltype_variability.tsv"))
    write_tab(out$sample_diag, file.path(run_dir, "sample_diagnostics.tsv"))
    write_tab(out$flags, file.path(run_dir, "degeneracy_flags.tsv"))
    summary <- data.table(dataset = ds, disease = input$disease,
                          cohort_role = input$cohort_role,
                          reference_version = REF_VERSION,
                          marker_set = MARKER_SET,
                          n_samples = nrow(est),
                          genes_used = length(input$genes),
                          n_pseudobulks = nrow(input$ref_log),
                          n_markers = sum(lengths(input$markers)),
                          gamma = as.numeric(fit$result$gamma),
                          elapsed_seconds = fit$elapsed_seconds,
                          checkpoint_reused = fit$checkpoint_reused,
                          degenerate_uniform = out$flags$degenerate_uniform,
                          degenerate_constant = out$flags$degenerate_constant,
                          unique_vectors_6dec = out$flags$unique_vectors_6dec,
                          median_uniform_rmse_pp = out$flags$uniform_rmse_median_pp)
    write_tab(summary, file.path(run_dir, "run_summary.tsv"))

    checks <- rbindlist(list(
      run_status(ds, "sample_mapping", "PASS", nrow(est)),
      run_status(ds, "gene_alignment", "PASS", length(input$genes)),
      run_status(ds, "reference_support", "PASS", nrow(ref_groups)),
      run_status(ds, "markers_unique", "PASS", sum(lengths(input$markers))),
      run_status(ds, "finite_estimates", "PASS", TRUE),
      run_status(ds, "proportions_in_unit_interval", "PASS", TRUE),
      run_status(ds, "row_sums_one", "PASS", max(abs(rowSums(est) - 1))),
      run_status(ds, "nonconstant_estimates",
                 if (out$flags$degenerate_constant) "WARN" else "PASS",
                 out$flags$unique_vectors_6dec),
      run_status(ds, "nonuniform_estimates",
                 if (out$flags$degenerate_uniform) "WARN" else "PASS",
                 out$flags$uniform_rmse_median_pp)
    ), fill = TRUE)
    write_tab(checks, file.path(run_dir, "status.tsv"))

    plot_status <- make_plots(out, input, plot_dir)
    plot_status[, `:=`(dataset = ds, reference_version = REF_VERSION,
                      marker_set = MARKER_SET)]
    setcolorder(plot_status, c("dataset", "reference_version", "marker_set",
                               "file", "status", "error"))
    write_tab(plot_status, file.path(run_dir, "08_plot_status.tsv"))
    if (any(plot_status$status == "WARN")) {
      warning(ds, ": um ou mais gráficos não foram gerados. Saídas numéricas preservadas.")
    }
    if (out$flags$degenerate_uniform || out$flags$degenerate_constant) {
      warning(ds, ": alerta de degenerescência; não usar como evidência sem auditoria.")
    }
    if (file.exists(file.path(run_dir, "ERROR.tsv")))
      unlink(file.path(run_dir, "ERROR.tsv"))
    say("CONCLUÍDO: %s | gamma=%s | degenerescência uniforme=%s",
        ds, format(summary$gamma, digits = 3), out$flags$degenerate_uniform)
    list(summary = summary, out = out, status = checks, plots = plot_status)
  }, error = function(e) {
    message("FALHA em ", ds, ": ", conditionMessage(e))
    dt <- data.table(dataset = ds, reference_version = REF_VERSION,
                     marker_set = MARKER_SET, error = conditionMessage(e))
    write_tab(dt, file.path(run_dir, "ERROR.tsv"))
    run_errors <<- c(run_errors, ds)
    NULL
  })
  if (!is.null(one)) {
    all_results[[ds]] <- one
    all_status[[ds]] <- one$status
    all_plot_status[[ds]] <- one$plots
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# 9) CONSOLIDAÇÃO E COMPARAÇÃO DESCRITIVA COM AS ETAPAS 06/07
# ──────────────────────────────────────────────────────────────────────────────
if (length(all_results)) {
  summary_all <- rbindlist(lapply(all_results, `[[`, "summary"), fill = TRUE)
  long_all <- rbindlist(lapply(all_results, function(a) a$out$long), fill = TRUE)
  grp_all <- rbindlist(lapply(all_results, function(a) a$out$group), fill = TRUE)
  var_all <- rbindlist(lapply(all_results, function(a) a$out$variation), fill = TRUE)
  diag_all <- rbindlist(lapply(all_results, function(a) a$out$sample_diag), fill = TRUE)
  status_all <- rbindlist(all_status, fill = TRUE)
  plots_all <- rbindlist(all_plot_status, fill = TRUE)
  write_tab(summary_all, file.path(OUT_ROOT, "08_run_summary.tsv"))
  write_tab(long_all, file.path(OUT_ROOT, "08_proportions_long_all.tsv"))
  write_tab(grp_all, file.path(OUT_ROOT, "08_proportion_summary_all.tsv"))
  write_tab(var_all, file.path(OUT_ROOT, "08_celltype_variability_all.tsv"))
  write_tab(diag_all, file.path(OUT_ROOT, "08_sample_diagnostics_all.tsv"))
  write_tab(status_all, file.path(OUT_ROOT, "08_status.tsv"))
  write_tab(plots_all, file.path(OUT_ROOT, "08_plot_status_all.tsv"))

  # Pares por MESMO dataset + sample_id + cell_type, nunca entre indivíduos
  # não pareados. Não calculamos p-valores aqui e não escolhemos método vencedor.
  pairs <- list()
  for (sp in list(list(path = MUSIC2_LONG, name = "MuSiC2"),
                  list(path = BISQUE_LONG, name = "BisqueRNA"))) {
    if (!file.exists(sp$path)) next
    other <- fread(sp$path)
    req <- c("dataset", "sample_id", "cell_type", "proportion",
             "reference_version", "marker_set")
    if (!all(req %in% names(other))) next
    other <- other[reference_version == REF_VERSION & marker_set == MARKER_SET,
                   .(dataset, sample_id, cell_type,
                     proportion_other = proportion)]
    if (anyDuplicated(other[, .(dataset, sample_id, cell_type)])) {
      warning("Comparação com ", sp$name, ": IDs duplicados; ignorando.")
      next
    }
    matched <- merge(long_all[, .(dataset, sample_id, cell_type,
                                   proportion_dtangle = proportion)],
                     other, by = c("dataset", "sample_id", "cell_type"),
                     all = FALSE, sort = FALSE)
    if (!nrow(matched)) next
    cmp <- matched[, .(n_samples = .N,
                       mean_absolute_difference_pp =
                         mean(abs(proportion_dtangle - proportion_other)) * 100,
                       mean_signed_difference_pp =
                         mean(proportion_dtangle - proportion_other) * 100,
                       spearman = safe_spearman(proportion_dtangle,
                                                proportion_other)),
                   by = .(dataset, cell_type)]
    cmp[, method_other := sp$name]
    pairs[[sp$name]] <- cmp
  }
  if (length(pairs)) write_tab(rbindlist(pairs, fill = TRUE),
                              file.path(OUT_ROOT, "08_vs_MuSiC2_BisqueRNA.tsv"))

  # Gráficos globais: zeros e variabilidade. Uma mesma escala por painel.
  global_plots <- list()
  try_plot <- function(key, expr) {
    global_plots[[key]] <<- tryCatch(expr, error = function(e)
      data.table(file = key, status = "WARN", error = conditionMessage(e)))
  }
  global_dir <- file.path(OUT_ROOT, "plots")
  z <- copy(var_all)
  z[, cell_type := factor(cell_type, levels = EXPECTED_TYPES)]
  pz <- ggplot(z, aes(x = cell_type, y = 100 * fraction_zero,
                      fill = cell_type)) +
    geom_col(width = 0.70, show.legend = FALSE) +
    facet_wrap(~ dataset, ncol = 2) +
    scale_fill_manual(values = CELL_COLORS) +
    scale_y_continuous(breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0.01, 0.05))) +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    labs(title = "dtangle • proporções exatamente iguais a zero",
         subtitle = "Comparação descritiva por coorte e tipo celular",
         x = "Tipo celular", y = "Amostras com estimativa zero (%)") + plot_theme(11)
  try_plot("zeros_all", save_plot(pz, file.path(global_dir,
                                 "08_zero_fraction_all.png"), 13, 8.5))
  ps <- ggplot(z, aes(x = cell_type, y = sd_pp, fill = cell_type)) +
    geom_col(width = 0.70, show.legend = FALSE) +
    geom_hline(yintercept = LOW_SD_PP, color = "#B42318", linetype = "dashed") +
    facet_wrap(~ dataset, ncol = 2) +
    scale_fill_manual(values = CELL_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(title = "dtangle • desvio-padrão das proporções",
         subtitle = "Linha tracejada: alerta de SD < 0,20 ponto percentual",
         x = "Tipo celular", y = "SD entre amostras (pp)") + plot_theme(11)
  try_plot("variability_all", save_plot(ps, file.path(global_dir,
                                   "08_variability_all.png"), 13, 8.5))
  write_tab(rbindlist(global_plots, fill = TRUE),
            file.path(OUT_ROOT, "08_global_plot_status.tsv"))
}

# Auditoria adicional dos 44 controles comuns. Não agregá-los como independentes.
if (all(c("GSE64810", "GSE68719") %in% names(all_results))) {
  h <- all_results[["GSE64810"]]$out$long[group == "Control",
                  .(sample_id, cell_type, hd = proportion)]
  p <- all_results[["GSE68719"]]$out$long[group == "Control",
                  .(sample_id, cell_type, pd = proportion)]
  shared <- merge(h, p, by = c("sample_id", "cell_type"))
  if (nrow(shared)) {
    shared[, difference_pp := 100 * (hd - pd)]
    write_tab(shared, file.path(OUT_ROOT, "08_shared_HD_PD_controls.tsv"))
    shared_stats <- shared[, .(n_shared = .N,
                              mean_absolute_difference_pp = mean(abs(difference_pp)),
                              max_absolute_difference_pp = max(abs(difference_pp)),
                              spearman = safe_spearman(hd, pd)), by = cell_type]
    write_tab(shared_stats,
              file.path(OUT_ROOT, "08_shared_HD_PD_controls_summary.tsv"))
  }
}

manifest <- data.table(path = list.files(OUT_ROOT, recursive = TRUE,
                                        full.names = TRUE))
manifest <- manifest[!dir.exists(path)]
if (nrow(manifest)) {
  manifest[, bytes := file.info(path)$size]
  write_tab(manifest, file.path(OUT_ROOT, "08_output_manifest.tsv"))
}
sink(file.path(OUT_ROOT, "08_sessionInfo.txt"))
print(sessionInfo())
sink()

say("\n=== ETAPA 08 FINALIZADA ===")
say("Coortes processadas: %s", paste(names(all_results), collapse = ", "))
say("Resultados: %s", OUT_ROOT)
if (length(run_errors)) {
  fail("Falha em ", paste(run_errors, collapse = ", "),
       ". Consulte ERROR.tsv na respectiva pasta; os demais ajustes foram preservados.")
}
say("Verificar: 08_run_summary.tsv, 08_status.tsv, 08_celltype_variability_all.tsv")
say("Atenção: saídas são estimativas relativas; não demonstram exatidão biológica.")
