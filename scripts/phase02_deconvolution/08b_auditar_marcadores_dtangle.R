#!/usr/bin/env Rscript
# ============================================================================
# 08b | AUDITORIA DOS MARCADORES E DA REFERÊNCIA DO DTANGLE — SEM REESTIMAR
# ============================================================================
# Objetivo: investigar estimativas vascular End+Per excessivamente elevadas
# nas quatro coortes. Esta auditoria NÃO executa dtangle, NÃO altera os scripts
# 04–08 e NÃO substitui proporções. Deve ser iniciada na raiz do repositório.
#
# Verificações:
# 1) Alinhamento núcleo -> tipo celular e gene -> referência;
# 2) Especificidade real de cada marcador top300 exclusivo (08);
# 3) Suporte por doador e tipo celular, sobretudo End/Per;
# 4) Média vascular estimada em cada coorte, somente descritiva.
#
# Atenção: a especificidade aqui é calculada sobre CPM agregado por tipo,
# e NÃO equivale exatamente ao critério de seleção da etapa 04 ou ao estimador
# interno do dtangle. É uma auditoria independente, não uma reclassificação.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
REF <- "data/reference/mathys2019/phase02_prepared/reference"
OUT08 <- "results/deconvolution/dtangle"
OUT <- file.path(OUT08, "08b_auditoria_marcadores")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(REF, "mathys2019_pathology_free_geneid_counts_sparse.rds")
meta_file <- file.path(REF, "mathys2019_pathology_free_cell_metadata_aligned.tsv")
marker_file <- file.path(OUT08, "08_markers_unique_assignment.tsv")
ref_support_file <- file.path(OUT08, "08_reference_pseudobulk_audit.tsv")
var_file <- file.path(OUT08, "08_celltype_variability_all.tsv")
paths <- c(counts_file, meta_file, marker_file, ref_support_file, var_file)
missing <- paths[!file.exists(paths)]
if (length(missing)) stop("Arquivos ausentes: ", paste(missing, collapse = "; "))

cat("08b | Auditoria independente da especificidade e suporte da referência\n")
sc <- readRDS(counts_file)
if (!inherits(sc, "Matrix")) sc <- Matrix(sc, sparse = TRUE)
if (is.null(rownames(sc)) || is.null(colnames(sc))) {
  stop("A matriz de referência não possui GeneIDs ou IDs nucleares.")
}
meta <- fread(meta_file)
stopifnot(all(c("TAG", "donor_id", "cell_type") %in% names(meta)))
idx <- match(colnames(sc), as.character(meta[["TAG"]]))
if (anyNA(idx)) stop("Núcleos não encontrados em metadata.")
meta <- meta[idx]
if (!setequal(unique(meta[["cell_type"]]), TYPES)) {
  stop("Tipos celulares distintos dos oito definidos no estudo.")
}
markers <- fread(marker_file)
stopifnot(all(c("gene_id", "cell_type") %in% names(markers)))
markers[, gene_id := as.character(gene_id)]
markers[, cell_type := as.character(cell_type)]
if (anyDuplicated(markers[["gene_id"]])) stop("Genes duplicados na atribuição exclusiva de marcadores.")
if (!all(markers[["cell_type"]] %in% TYPES)) stop("Tipo celular inválido nos marcadores.")
if (!all(markers[["gene_id"]] %in% rownames(sc))) {
  stop("Alguns marcadores exclusivos não estão na referência.")
}

# CPM calculado com a biblioteca COMPLETA de cada tipo, não com top300.
# Agrupar as contagens primeiro reduz o efeito de núcleos individuais rasos.
genes <- markers[["gene_id"]]
indices <- match(genes, rownames(sc))
ref_markers <- sc[indices, , drop = FALSE]
type_cpm <- matrix(NA_real_, nrow = length(genes), ncol = length(TYPES),
                   dimnames = list(genes, TYPES))
support <- vector("list", length(TYPES))
for (j in seq_along(TYPES)) {
  typ <- TYPES[[j]]
  which_type <- which(as.character(meta[["cell_type"]]) == typ)
  denom <- sum(Matrix::colSums(sc[, which_type, drop = FALSE]))
  if (!is.finite(denom) || denom <= 0) stop("Biblioteca inválida: ", typ)
  numer <- Matrix::rowSums(ref_markers[, which_type, drop = FALSE])
  type_cpm[, j] <- as.numeric(numer) / denom * 1e6
  support[[j]] <- data.table(cell_type = typ,
                              n_nuclei = length(which_type),
                              n_donors = uniqueN(meta[["donor_id"]][which_type]),
                              library_counts = denom)
}

# Para cada gene, comparar expressão no tipo designado com o melhor concorrente.
assigned_idx <- match(markers[["cell_type"]], TYPES)
target <- type_cpm[cbind(seq_len(nrow(markers)), assigned_idx)]
tmp <- type_cpm
# Descartar a coluna do tipo-alvo da busca pelo maior concorrente.
tmp[cbind(seq_len(nrow(markers)), assigned_idx)] <- -Inf
rival_idx <- max.col(tmp, ties.method = "first")
rival <- tmp[cbind(seq_len(nrow(markers)), rival_idx)]
if (any(!is.finite(target) | !is.finite(rival))) stop("CPM inválidos.")

specificity <- copy(markers)
specificity[, `:=`(
  target_cpm = as.numeric(target),
  largest_other_cpm = as.numeric(rival),
  strongest_other_cell_type = TYPES[rival_idx],
  target_is_highest = as.logical(target >= rival),
  target_over_other_ratio_pseudocount1 = as.numeric((target + 1) / (rival + 1)),
  log2_target_over_other_pseudocount1 = as.numeric(log2((target + 1) / (rival + 1)))
)]
setorder(specificity, cell_type, log2_target_over_other_pseudocount1)
fwrite(specificity, file.path(OUT, "08b_especificidade_por_marcador.tsv"), sep = "\t")

summary_marker <- specificity[, .(
  n_markers = .N,
  fraction_highest_in_assigned_type = mean(target_is_highest),
  fraction_target_cpm_below_1 = mean(target_cpm < 1),
  median_target_cpm = median(target_cpm),
  median_strongest_other_cpm = median(largest_other_cpm),
  median_log2_target_over_other = median(log2_target_over_other_pseudocount1),
  fraction_target_expression_below_competitor = mean(target_cpm < largest_other_cpm)
), by = cell_type]
summary_marker[, cell_type := factor(cell_type, levels = TYPES)]
setorder(summary_marker, cell_type)
summary_marker[, cell_type := as.character(cell_type)]
fwrite(summary_marker, file.path(OUT, "08b_resumo_especificidade.tsv"), sep = "\t")

# Marcadores com a pior separação para End/Per são um ponto de inspeção,
# e NÃO são removidos automaticamente do estudo.
weak_vascular <- specificity[cell_type %in% c("End", "Per")][
  order(log2_target_over_other_pseudocount1)
]
weak_vascular <- head(weak_vascular, 40L)
fwrite(weak_vascular, file.path(OUT, "08b_marcadores_vasculares_inspecao.tsv"), sep = "\t")

pb <- fread(ref_support_file)
stopifnot(all(c("cell_type", "n_nuclei", "donor_id") %in% names(pb)))

# IMPORTANTE (correção 08b-v2): data.table exige que cada coluna de uma
# agregação `by=` tenha a MESMA classe de retorno em todos os grupos.
# stats::median() sobre contagens inteiras pode retornar `integer` para
# um grupo e `double` para outro, causando:
#   Column 4 of result for group ... type 'integer' but expecting 'double'
# Converter as contagens para `numeric` e definir explicitamente as classes
# das estatísticas resolve o problema SEM alterar seus valores.
pb[, n_nuclei := as.numeric(n_nuclei)]
if (anyNA(pb[["n_nuclei"]]) || any(!is.finite(pb[["n_nuclei"]])) ||
    any(pb[["n_nuclei"]] < 0)) {
  stop("08b: n_nuclei contém NA, valores não finitos ou negativos.")
}

pb_support <- pb[, .(
  n_observed_donors = as.integer(.N),
  total_nuclei = as.numeric(sum(n_nuclei)),
  min_nuclei_per_donor = as.numeric(min(n_nuclei)),
  median_nuclei_per_donor = as.numeric(stats::median(n_nuclei)),
  n_donors_with_1_or_2_nuclei = as.integer(sum(n_nuclei <= 2)),
  n_donors_with_5_or_fewer_nuclei = as.integer(sum(n_nuclei <= 5))
), by = cell_type]
fwrite(pb_support, file.path(OUT, "08b_suporte_por_tipo.tsv"), sep = "\t")

# Extrair a fração vascular estimada por coorte (média de proporções).
v <- fread(var_file)
stopifnot(all(c("dataset", "cell_type", "mean_proportion", "sd_pp", "fraction_zero") %in% names(v)))
vascular <- v[cell_type %in% c("End", "Per")]
vascular_wide <- dcast(vascular, dataset ~ cell_type, value.var = "mean_proportion")
vascular_wide[, `:=`(
  vascular_mean_fraction = End + Per,
  vascular_mean_percent = 100 * (End + Per)
)]
fwrite(vascular_wide, file.path(OUT, "08b_contribuicao_vascular_por_coorte.tsv"), sep = "\t")

# Documento de interpretação e limites (não classificar como causa demonstrada).
report <- c(
  "08b | AUDITORIA TÉCNICA DA ETAPA 08",
  "=================================",
  "Validação computacional dos resultados não implica validação biológica.",
  "Especificidade calculada como CPM de contagens agregadas por tipo.",
  "A seleção do top300 da etapa 04 usa critério distinto: não são equivalentes.",
  "End/Per são escassos e pseudobulks com pouquíssimos núcleos podem ser instáveis.",
  "Marcadores com expressão mais alta fora da categoria atribuída merecem revisão.",
  "Nenhum marcador foi excluído e nenhuma estimativa dtangle foi recalculada.",
  "Correlação entre métodos não prova acurácia. Não interpretar percentuais",
  "vasculares elevados sem validação de expressão/marcadores e sensibilidade.",
  "", "Resultados de especificidade:",
  capture.output(print(summary_marker)),
  "", "Composição vascular estimada por coorte:",
  capture.output(print(vascular_wide))
)
writeLines(report, file.path(OUT, "08b_relatorio_interpretacao.txt"), useBytes = TRUE)

cat("\n--- Especificidade por tipo celular ---\n")
print(summary_marker)
cat("\n--- Suporte por doador ---\n")
print(pb_support)
cat("\n--- End + Per por coorte (estimativa; não verdade biológica) ---\n")
print(vascular_wide)
cat("\nAUDITORIA 08b CONCLUÍDA. Arquivos em: ", OUT, "\n", sep = "")
