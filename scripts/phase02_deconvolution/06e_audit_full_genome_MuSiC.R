#!/usr/bin/env Rscript
# ==============================================================================
# 06e_audit_full_genome_MuSiC.R
# FASE 02 - AUDITORIA DIAGNÓSTICA (SEM ALTERAR A ANÁLISE PRINCIPAL)
# ==============================================================================
# OBJETIVO
#  Comparar, somente nos controles, os resultados de MuSiC convencional usando:
#  (A) top300 por tipo celular (a união de 2.242 GeneIDs, ou genes disponíveis);
#  (B) todos os genes compartilhados entre bulk e referência Mathys.
#
#  Para AMBAS as entradas, avaliar duas saídas do mesmo ajuste:
#     Est.prop.weighted -> estimativa ponderada do MuSiC;
#     Est.prop.allgene  -> estimativa NNLS de referência, sem pesos do MuSiC.
#
#  Isso permite testar SEPARADAMENTE:
#   1. efeito da seleção de genes (top300 vs todos os genes);
#   2. efeito da ponderação interna (MuSiC ponderado vs NNLS).
#
#  Também audita diferenças de contagens por núcleo (proxy de tamanho de
#  biblioteca, não medida biológica de tamanho da célula).
#
#  REQUISITOS/RESTRIÇÕES
#  - mantém 8 tipos, 24 doadores, mesmos controles e parâmetros do 06d;
#  - usa raw counts e normalização dentro do MuSiC igual a 06d;
#  - não preenche genes não observados com zero;
#  - cada dataset possui sua própria interseção com a referência;
#  - não chama MuSiC2, não substitui outputs de 06/06c/06d;
#  - não realiza testes caso-controle nem prova acurácia biológica.
#
# EXECUTAR DA RAIZ DO REPOSITÓRIO:
#   Rscript scripts/phase02_deconvolution/06e_audit_full_genome_MuSiC.R
#
# SAÍDAS:
# results/deconvolution/MuSiC_FullGene_Audit/
#   06e_proportions_long.tsv
#   06e_summary_by_dataset_celltype.tsv
#   06e_paired_contrasts.tsv
#   06e_input_audit.tsv
#   06e_reference_library_by_celltype.tsv
#   06e_top300_reproduction.tsv
#   06e_status.tsv
#   06e_sessionInfo.txt
#   plots/06e_zero_heatmap.png
#   plots/06e_mean_proportion.png
# ==============================================================================

DATASETS <- c("ROSMAP", "GSE125583")
REFERENCE_VERSION <- "main24"
CELL_TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
GENE_SETS <- c("top300", "allgenes")
NEAR_ZERO <- 0.001 # fracao: 0.001 = 0.1%
ZERO_TOL <- 0
EXPECTED_DONORS <- 24L
EXPECTED_NUCLEI <- 35140L
SUM_TOLERANCE <- 1e-6

pkgs <- c("data.table", "Matrix", "MuSiC", "SingleCellExperiment", "S4Vectors")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Pacotes não instalados: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(MuSiC)
  library(SingleCellExperiment)
  library(S4Vectors)
})

outdir <- file.path("results", "deconvolution", "MuSiC_FullGene_Audit")
plotdir <- file.path(outdir, "plots")
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)
refdir <- file.path("data", "reference", "mathys2019", "phase02_prepared", "reference")
refcountsfile <- file.path(refdir, "mathys2019_pathology_free_geneid_counts_sparse.rds")
refmetafile <- file.path(refdir, "mathys2019_pathology_free_cell_metadata_aligned.tsv")
bulkroot <- file.path("data", "processed", "deconvolution")
prior_music2 <- file.path("results", "deconvolution", "MuSiC2", "06_proportions_long_all.tsv")

stop_if_missing <- function(paths) {
  m <- paths[!file.exists(paths)]
  if (length(m)) stop("Arquivos ausentes: ", paste(m, collapse = "; "))
}
stop_if_missing(c(refcountsfile, refmetafile, prior_music2))
previous <- fread(prior_music2, na.strings = c("NA", ""))
if (!all(c("dataset", "reference_version", "marker_set", "group", "sample_id",
           "cell_type", "proportion") %in% names(previous))) {
  stop("Tabela de proporções 06 sem colunas esperadas.")
}

cat("Carregando Mathys main24...\n")
ref <- readRDS(refcountsfile)
if (!inherits(ref, "dgCMatrix")) ref <- as(ref, "dgCMatrix")
meta_ref <- fread(refmetafile, na.strings = c("NA", ""))
if (!all(c("TAG", "donor_id", "cell_type") %in% names(meta_ref))) {
  stop("Metadados Mathys não têm TAG, donor_id ou cell_type.")
}
if (is.null(rownames(ref)) || is.null(colnames(ref)) || anyDuplicated(colnames(ref)) ||
    anyDuplicated(rownames(ref)) || anyDuplicated(meta_ref$TAG)) {
  stop("Dimnames GeneID/TAG da referência ausentes ou duplicados.")
}
idx <- match(colnames(ref), as.character(meta_ref$TAG))
if (anyNA(idx)) stop("Falha no alinhamento dos TAGs na referência.")
meta_ref <- meta_ref[idx]
if (!identical(as.character(meta_ref$TAG), colnames(ref)) ||
    uniqueN(meta_ref$donor_id) != EXPECTED_DONORS ||
    ncol(ref) != EXPECTED_NUCLEI ||
    !setequal(unique(meta_ref$cell_type), CELL_TYPES)) {
  stop("A referência carregada não corresponde ao main24 auditado.")
}

# Apenas DIAGNÓSTICO das bibliotecas de núcleo (não é 'cell size' real).
ref_total_counts <- Matrix::colSums(ref)
ref_genes_detected <- Matrix::colSums(ref > 0)
ref_lib <- data.table(
  cell_type = as.character(meta_ref$cell_type),
  donor_id = as.character(meta_ref$donor_id),
  total_counts_per_nucleus = as.numeric(ref_total_counts),
  detected_genes_per_nucleus = as.numeric(ref_genes_detected)
)
ref_lib_sum <- ref_lib[, .(
  n_nuclei = .N,
  n_donors = uniqueN(donor_id),
  library_counts_median = median(total_counts_per_nucleus),
  library_counts_mean = mean(total_counts_per_nucleus),
  detected_genes_median = median(detected_genes_per_nucleus)
), by = cell_type]
ref_lib_sum[, order_i := match(cell_type, CELL_TYPES)]
setorder(ref_lib_sum, order_i)
ref_lib_sum[, order_i := NULL]
fwrite(ref_lib_sum, file.path(outdir, "06e_reference_library_by_celltype.tsv"), sep = "\t")

standardize <- function(x, sample_ids, name) {
  if (is.null(x)) stop(name, ": campo de proporções ausente.")
  p <- as.matrix(x)
  if (!is.null(rownames(p)) && all(sample_ids %in% rownames(p))) {
    p <- p[sample_ids, , drop = FALSE]
  } else if (!is.null(colnames(p)) && all(sample_ids %in% colnames(p))) {
    p <- t(p[, sample_ids, drop = FALSE])
  } else stop(name, ": não é possível alinhar as amostras.")
  if (!setequal(colnames(p), CELL_TYPES)) {
    stop(name, ": MuSiC retornou tipos celulares inesperados.")
  }
  p <- p[, CELL_TYPES, drop = FALSE]
  storage.mode(p) <- "double"
  if (anyNA(p) || any(!is.finite(p)) || any(p < -1e-8)) {
    stop(name, ": proporções inválidas.")
  }
  if (max(abs(rowSums(p) - 1)) > SUM_TOLERANCE) {
    stop(name, ": soma de proporções difere de 1. Não renormalizar silenciosamente.")
  }
  p
}

as_long <- function(p, ds, gene_set, method) {
  x <- data.table::melt(
    as.data.table(p, keep.rownames = "sample_id"),
    id.vars = "sample_id", variable.name = "cell_type", value.name = "proportion"
  )
  x[, `:=`(dataset = ds, gene_set = gene_set, method = method,
           cell_type = as.character(cell_type))]
  setcolorder(x, c("dataset", "gene_set", "method", "sample_id", "cell_type", "proportion"))
  x
}

all_results <- list()
input_logs <- list()
reproduction_checks <- list()
fit_count <- 0L
for (ds in DATASETS) {
  metafile <- file.path(bulkroot, ds, "metadata_deconvolution.tsv")
  stop_if_missing(metafile)
  bulk_meta <- fread(metafile, na.strings = c("NA", ""))
  if (!all(c("sample_id", "group") %in% names(bulk_meta)) ||
      anyDuplicated(bulk_meta$sample_id)) stop(ds, ": metadata bulk inválido.")
  ids <- as.character(bulk_meta[group == "Control", sample_id])
  if (length(ids) < 2L) stop(ds, ": menos de 2 controles.")
  cat("\nDataset: ", ds, " | controles: ", length(ids), "\n", sep = "")

  for (gs in GENE_SETS) {
    bulkfile <- file.path(
      bulkroot, ds,
      if (gs == "top300") "bulk_counts_top300.rds" else "bulk_counts_reference_intersection.rds"
    )
    stop_if_missing(bulkfile)
    bulk <- as.matrix(readRDS(bulkfile))
    storage.mode(bulk) <- "double"
    if (is.null(rownames(bulk)) || is.null(colnames(bulk)) ||
        anyDuplicated(rownames(bulk)) || anyDuplicated(colnames(bulk)) ||
        anyNA(bulk) || any(!is.finite(bulk)) || any(bulk < 0)) {
      stop(ds, " / ", gs, ": matriz bulk inválida.")
    }
    if (!all(ids %in% colnames(bulk))) stop(ds, ": controles ausentes no bulk.")
    common <- rownames(bulk)[rownames(bulk) %in% rownames(ref)]
    if (length(common) < 100L) stop(ds, " / ", gs, ": <100 genes compartilhados.")
    r <- ref[common, , drop = FALSE]
    b_full <- bulk[common, , drop = FALSE]
    # Filtrar como no 06: zero somente se zerado em TODO o bulk (controles + casos)
    # ou em toda a referência main24. NÃO filtrar por resultado de doença.
    keep <- rowSums(b_full) > 0 & Matrix::rowSums(r) > 0
    genes <- common[keep]
    if (length(genes) < 100L) stop(ds, " / ", gs, ": <100 genes após QC.")
    sc <- r[genes, , drop = FALSE]
    b <- b_full[genes, ids, drop = FALSE]
    sce <- SingleCellExperiment(
      assays = list(counts = sc),
      colData = S4Vectors::DataFrame(
        cellType = factor(meta_ref$cell_type, levels = CELL_TYPES),
        sampleID = factor(meta_ref$donor_id),
        row.names = as.character(meta_ref$TAG)
      )
    )
    rownames(sce) <- genes
    colnames(sce) <- meta_ref$TAG
    if (!identical(rownames(b), rownames(sce))) {
      stop(ds, " / ", gs, ": referência e bulk em ordens diferentes.")
    }
    cat("  ", gs, ": ", length(genes), " genes; iniciando music_prop()...\n", sep = "")
    begin <- Sys.time()
    # Parâmetros equivalentes ao 06d, para comparação controlada.
    fit <- MuSiC::music_prop(
      bulk.mtx = b, sc.sce = sce,
      clusters = "cellType", samples = "sampleID", select.ct = CELL_TYPES,
      markers = NULL, cell_size = NULL, ct.cov = FALSE,
      iter.max = 1000, nu = 0.0001, eps = 0.01,
      centered = FALSE, normalize = FALSE, verbose = FALSE
    )
    elapsed <- as.numeric(difftime(Sys.time(), begin, units = "secs"))
    p_music <- standardize(fit$Est.prop.weighted, ids, paste(ds, gs, "MuSiC"))
    p_nnls <- standardize(fit$Est.prop.allgene, ids, paste(ds, gs, "NNLS"))
    if (gs == "top300") {
      expected <- previous[
        dataset == ds & reference_version == "main24" &
          marker_set == "top300" & group == "Control",
        .(sample_id, cell_type, proportion)
      ]
      expected_mat <- data.table::dcast(expected, sample_id ~ cell_type,
                                         value.var = "proportion")
      expected_mat <- expected_mat[match(ids, expected_mat$sample_id)]
      if (nrow(expected_mat) != length(ids) || anyNA(expected_mat$sample_id) ||
          !identical(as.character(expected_mat$sample_id), ids)) {
        stop(ds, ": falha ao alinhar top300 com o script 06.")
      }
      p_expected <- as.matrix(expected_mat[, ..CELL_TYPES])
      storage.mode(p_expected) <- "double"
      maxdiff <- max(abs(p_expected - p_music))
      reproduction_checks[[length(reproduction_checks) + 1L]] <- data.table(
        dataset = ds, gene_set = gs, max_abs_diff_vs_06 = maxdiff,
        status = if (maxdiff <= 1e-8) "PASS" else "FAIL"
      )
      if (maxdiff > 1e-8) {
        stop(ds, ": top300 MuSiC não reproduziu as estimativas controle do script 06; diferença máxima = ", maxdiff)
      }
    }
    fit_count <- fit_count + 1L
    all_results[[length(all_results) + 1L]] <- as_long(p_music, ds, gs, "MuSiC_weighted")
    all_results[[length(all_results) + 1L]] <- as_long(p_nnls, ds, gs, "NNLS_unweighted")
    input_logs[[length(input_logs) + 1L]] <- data.table(
      dataset = ds, gene_set = gs, controls = length(ids),
      shared_geneids_initial = length(common), genes_used = length(genes),
      removed_allzero_bulk = sum(rowSums(b_full) == 0),
      removed_allzero_reference = sum(Matrix::rowSums(r) == 0),
      reference_donors = uniqueN(meta_ref$donor_id), reference_nuclei = ncol(sc),
      elapsed_seconds = elapsed, MuSiC_version = as.character(packageVersion("MuSiC"))
    )
    cat("  ", gs, " concluído em ", round(elapsed, 1), " s\n", sep = "")
    # Salvar resultados parciais após cada ajuste para proteger o trabalho.
    partial <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
    fwrite(partial, file.path(outdir, "06e_proportions_partial.tsv"), sep = "\t", na = "NA")
    fwrite(rbindlist(input_logs), file.path(outdir, "06e_input_audit_partial.tsv"), sep = "\t", na = "NA")
    rm(bulk, b_full, r, sc, b, sce, fit, p_music, p_nnls, partial)
    invisible(gc())
  }
}

long <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
if (fit_count != length(DATASETS) * length(GENE_SETS)) {
  stop("Nem todos os ajustes esperados foram realizados.")
}
if (anyDuplicated(long[, .(dataset, gene_set, method, sample_id, cell_type)])) {
  stop("Duplicação de (dataset, gene_set, method, sample_id, cell_type).")
}
summary <- long[, .(
  n_samples = .N,
  mean_proportion = mean(proportion),
  median_proportion = median(proportion),
  n_zero = sum(proportion == ZERO_TOL),
  fraction_zero = mean(proportion == ZERO_TOL),
  fraction_near_zero = mean(proportion < NEAR_ZERO)
), by = .(dataset, gene_set, method, cell_type)]

# Comparações pareadas: cada contraste modifica SOMENTE uma dimensão.
base <- long[gene_set == "top300" & method == "MuSiC_weighted",
             .(dataset, sample_id, cell_type, base_top300_weighted = proportion)]
allc <- merge(long, base, by = c("dataset", "sample_id", "cell_type"), all.x = TRUE)
if (anyNA(allc$base_top300_weighted)) stop("Falha no pareamento com top300 weighted.")
paired <- allc[, .(
  n_samples = .N,
  mean_absolute_difference_pp = 100 * mean(abs(proportion - base_top300_weighted)),
  median_difference_pp = 100 * median(proportion - base_top300_weighted),
  spearman = if (sd(proportion) > 0 && sd(base_top300_weighted) > 0) {
    suppressWarnings(cor(proportion, base_top300_weighted, method = "spearman"))
  } else NA_real_
), by = .(dataset, gene_set, method, cell_type)]

fwrite(long, file.path(outdir, "06e_proportions_long.tsv"), sep = "\t", na = "NA")
fwrite(summary, file.path(outdir, "06e_summary_by_dataset_celltype.tsv"), sep = "\t", na = "NA")
fwrite(paired, file.path(outdir, "06e_paired_contrasts.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(input_logs), file.path(outdir, "06e_input_audit.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(reproduction_checks), file.path(outdir, "06e_top300_reproduction.tsv"), sep = "\t", na = "NA")
status <- data.table(
  check = c("all_four_fits_completed", "all_eight_types_present", "all_proportions_nonnegative", "all_proportions_sum_one", "top300_reproduces_MuSiC2_controls"),
  status = c(
    if (fit_count == 4L) "PASS" else "FAIL",
    if (all(long[, .(types = uniqueN(cell_type)), by = .(dataset, gene_set, method)]$types == 8L)) "PASS" else "FAIL",
    if (all(long$proportion >= -1e-8)) "PASS" else "FAIL",
    if (all(abs(long[, .(total = sum(proportion)), by = .(dataset, gene_set, method, sample_id)]$total - 1) <= SUM_TOLERANCE)) "PASS" else "FAIL",
    if (length(reproduction_checks) == length(DATASETS) &&
        all(rbindlist(reproduction_checks)$status == "PASS")) "PASS" else "FAIL"
  ),
  note = rep(NA_character_, 5L)
)
fwrite(status, file.path(outdir, "06e_status.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(outdir, "06e_sessionInfo.txt"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  plotzero <- copy(summary)
  plotzero[, condition := paste(gene_set, method, sep = " / ")]
  plotzero[, condition := factor(condition, levels = c(
    "top300 / MuSiC_weighted", "allgenes / MuSiC_weighted",
    "top300 / NNLS_unweighted", "allgenes / NNLS_unweighted"
  ))]
  p0 <- ggplot2::ggplot(plotzero, ggplot2::aes(
    x = condition, y = 100 * fraction_zero, fill = condition
  )) +
    ggplot2::geom_col(width = 0.76) +
    ggplot2::facet_grid(cell_type ~ dataset) +
    ggplot2::scale_fill_manual(values = c(
      "top300 / MuSiC_weighted" = "#64748B",
      "allgenes / MuSiC_weighted" = "#00ADA2",
      "top300 / NNLS_unweighted" = "#7B2CBF",
      "allgenes / NNLS_unweighted" = "#E75480"
    ), drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::labs(
      title = "06e: zeros por seleção de genes e ponderação",
      subtitle = "Controles; oito tipos celulares; referência Mathys main24",
      x = "Genes / estimador", y = "Estimativas iguais a zero (%)",
      caption = "Diagnóstico, não validação de acurácia nem estimativa de diferenças caso-controle"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, hjust = 1),
                   legend.position = "none", panel.grid.minor = ggplot2::element_blank())
  ggplot2::ggsave(file.path(plotdir, "06e_zero_heatmap.png"), p0,
                  width = 13, height = 11, dpi = 180)

  meanplot <- copy(summary[cell_type %in% c("Ex", "In", "Oli", "Opc")])
  meanplot[, condition := paste(gene_set, method, sep = " / ")]
  p1 <- ggplot2::ggplot(meanplot, ggplot2::aes(
    x = condition, y = 100 * mean_proportion, fill = condition
  )) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_grid(cell_type ~ dataset, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(
      "top300 / MuSiC_weighted" = "#64748B",
      "allgenes / MuSiC_weighted" = "#00ADA2",
      "top300 / NNLS_unweighted" = "#7B2CBF",
      "allgenes / NNLS_unweighted" = "#E75480"
    )) +
    ggplot2::labs(
      title = "06e: proporções médias por genes e ponderação",
      subtitle = "Cada facet possui sua própria escala; observar o eixo Y de In/Opc",
      x = "Genes / estimador", y = "Proporção média (%)",
      caption = "Resultados provisórios, somente controles"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, hjust = 1),
                   legend.position = "none", panel.grid.minor = ggplot2::element_blank())
  ggplot2::ggsave(file.path(plotdir, "06e_mean_proportion.png"), p1,
                  width = 13, height = 11, dpi = 180)
} else {
  message("ggplot2 ausente: tabelas foram geradas, gráficos omitidos.")
}

if (any(status$status == "FAIL")) stop("A auditoria 06e terminou com falhas; confira 06e_status.tsv.")
cat("\n=== 06e CONCLUÍDO ===\n")
print(rbindlist(input_logs)[, .(dataset, gene_set, controls, genes_used, elapsed_seconds)])
cat("\nResumo In/Opc:\n")
print(summary[cell_type %in% c("In", "Opc"),
              .(dataset, gene_set, method, n_zero, fraction_zero, mean_proportion)])
cat("\nArquivos em: ", outdir, "\n", sep = "")
