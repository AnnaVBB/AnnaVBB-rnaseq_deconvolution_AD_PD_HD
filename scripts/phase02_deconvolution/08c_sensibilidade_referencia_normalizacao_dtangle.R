#!/usr/bin/env Rscript
# =============================================================================
# 08c | dtangle | SENSIBILIDADE AO SUPORTE VASCULAR E NORMALIZACAO
# TCC AD-PD-HD / Mathys 2019 pathology-free
# =============================================================================
# OBJETIVO
#   Investigar, SEM SUBSTITUIR o resultado 08, dois fatores que podem explicar
#   estimativas End+Per de 56-62%:
#      (A) pseudobulks de doador/tipo com 1-2 nucleos vasculares;
#      (B) ausencia de normalizacao quantilica conjunta bulk + referencia.
#
# Quatro condicoes, em cada coorte:
#   1. donor_raw       = 181 pseudobulks (reproduz o script 08 original)
#   2. donor_quantile  = 181 pseudobulks + normalizacao quantilica conjunta
#   3. pooled_raw      = 8 pseudobulks (contagens reunidas por tipo)
#   4. pooled_quantile = 8 pseudobulks + normalizacao quantilica conjunta
#
# OBSERVACOES IMPORTANTES
# - Mantem OS OITO TIPOS e as atribuicoes top300 feitas no script 08.
# - O pooling por tipo privilegia nucleos de doadores mais amostrados; nao e
#   necessariamente metodologicamente superior. E apenas uma sensibilidade.
# - Normalizacao quantilica tambem pode remover diferencas biologicas reais;
#   nao escolher a condicao com base no percentual vascular mais plausivel.
# - Este script NAO edita nem recalcula os resultados existentes do 08.
# - Comeca por reproduzir o baseline 08 com tolerancia estrita (< 1e-6).
# - Nenhuma significancia de casos/controles e calculada aqui.
# - Resultados vao EXCLUSIVAMENTE a results/deconvolution/dtangle/08c_sensitivity/.
#
# EXECUCAO (raiz do repositorio, com renv ativo)
#   Rscript scripts/phase02_deconvolution/08c_sensibilidade_referencia_normalizacao_dtangle.R
#   DTANGLE08C_DATASET=GSE64810 Rscript scripts/phase02_deconvolution/08c_sensibilidade_referencia_normalizacao_dtangle.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(dtangle)
  library(limma)
  library(ggplot2)
})

TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
ALL_COHORTS <- c("ROSMAP", "GSE125583", "GSE64810", "GSE68719")
chosen <- trimws(Sys.getenv("DTANGLE08C_DATASET", ""))
COHORTS <- if (nzchar(chosen)) trimws(strsplit(chosen, ",", fixed=TRUE)[[1]]) else ALL_COHORTS
if (any(!COHORTS %in% ALL_COHORTS)) stop("DTANGLE08C_DATASET invalido")
COHORTS <- unique(COHORTS)
BASE <- "results/deconvolution/dtangle"
OUT <- file.path(BASE, "08c_sensitivity")
REF_DIR <- "data/reference/mathys2019/phase02_prepared/reference"
BULK_DIR <- "data/processed/deconvolution"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

required_paths <- c(
  file.path(REF_DIR, "mathys2019_pathology_free_geneid_counts_sparse.rds"),
  file.path(REF_DIR, "mathys2019_pathology_free_cell_metadata_aligned.tsv"),
  file.path(BASE, "08_reference_pseudobulk_audit.tsv"),
  file.path(BASE, "08_markers_unique_assignment.tsv")
)
if (any(!file.exists(required_paths))) stop("Arquivos da referencia/08 ausentes: ",
  paste(required_paths[!file.exists(required_paths)], collapse="; "))

cat("08c | Reconstruindo referencia do 08 de modo reproduzivel\n")
sc <- readRDS(required_paths[[1]])
if (!inherits(sc, "dgCMatrix")) sc <- as(sc, "dgCMatrix")
stopifnot(nrow(sc) == 17394L, ncol(sc) == 35140L)
meta <- fread(required_paths[[2]])
stopifnot(all(c("TAG", "donor_id", "cell_type") %in% names(meta)))
m_idx <- match(colnames(sc), as.character(meta[["TAG"]]))
if (anyNA(m_idx)) stop("Metadados da referencia nao correspondem aos nucleos")
meta <- meta[m_idx]
stopifnot(identical(as.character(meta[["TAG"]]), colnames(sc)))
refgrp <- fread(required_paths[[3]])
stopifnot(all(c("donor_id", "cell_type", "ref_id", "n_nuclei") %in% names(refgrp)))
stopifnot(nrow(refgrp) == 181L)
stopifnot(setequal(unique(refgrp$cell_type), TYPES))
ids_ref <- paste(refgrp$donor_id, refgrp$cell_type, sep="\r")
ids_sc <- paste(meta$donor_id, meta$cell_type, sep="\r")
grp_index <- match(ids_sc, ids_ref)
if (anyNA(grp_index)) stop("Nucleo sem grupo de referencia")
if (!identical(as.integer(tabulate(grp_index, nrow(refgrp))), as.integer(refgrp$n_nuclei))) {
  stop("Contagem de nucleos por grupo nao reproduz a etapa 08")
}
A <- Matrix::sparseMatrix(i=seq_len(ncol(sc)), j=grp_index,
  x=rep(1, ncol(sc)), dims=c(ncol(sc), nrow(refgrp)))
ref_pb <- sc %*% A
rownames(ref_pb) <- rownames(sc)
colnames(ref_pb) <- refgrp$ref_id
ref_lib <- as.numeric(Matrix::colSums(ref_pb))
if (any(ref_lib <= 0) || any(!is.finite(ref_lib))) stop("Biblioteca de pseudobulk invalida")

# Construir um perfil por tipo COM contagens somadas antes de normalizar.
ref_pooled <- sapply(TYPES, function(ct) {
  as.numeric(Matrix::rowSums(ref_pb[, refgrp$cell_type == ct, drop=FALSE]))
})
rownames(ref_pooled) <- rownames(ref_pb)
ref_pooled_lib <- as.numeric(colSums(ref_pooled))
stopifnot(all(ref_pooled_lib > 0))
pure_donor <- setNames(lapply(TYPES, function(ct) {
  which(refgrp$cell_type == ct)
}), TYPES)
pure_pool <- setNames(lapply(seq_along(TYPES), function(i) as.integer(i)), TYPES)

marker_tbl <- fread(required_paths[[4]])
stopifnot(all(c("gene_id", "cell_type") %in% names(marker_tbl)))
marker_tbl[, `:=`(gene_id=as.character(gene_id), cell_type=as.character(cell_type))]
if (anyDuplicated(marker_tbl$gene_id) || !setequal(unique(marker_tbl$cell_type), TYPES)) {
  stop("Marcadores exclusivos do 08 invalidos")
}

results <- list()
ct_rows <- list()
check_rows <- list()
failures <- list()

qnorm_both <- function(ref, bulk) {
  # limma espera genes nas linhas e AMOSTRAS nas colunas.
  joint <- cbind(t(ref), t(bulk))
  q <- limma::normalizeBetweenArrays(joint, method="quantile")
  nr <- nrow(ref)
  list(ref=t(q[, seq_len(nr), drop=FALSE]),
       bulk=t(q[, nr + seq_len(nrow(bulk)), drop=FALSE]))
}

validate_estimates <- function(mat, sample_ids, label) {
  e <- as.matrix(mat)
  if (!is.numeric(e) || nrow(e) != length(sample_ids) || ncol(e) != 8L) {
    stop(label, ": estimativas com dimensoes incorretas")
  }
  if (!is.null(colnames(e))) {
    if (!setequal(colnames(e), TYPES)) stop(label, ": tipos celulares incorretos")
    e <- e[, TYPES, drop=FALSE]
  }
  if (!is.null(rownames(e)) && setequal(rownames(e), sample_ids)) {
    e <- e[sample_ids, , drop=FALSE]
  }
  rownames(e) <- sample_ids
  colnames(e) <- TYPES
  if (any(!is.finite(e)) || min(e) < -1e-8 || max(e) > 1+1e-8 ||
      max(abs(rowSums(e) - 1)) > 1e-6) {
    stop(label, ": estimativas invalidas ou nao somam 1")
  }
  e
}

for (cohort in COHORTS) {
  cat("\n----------------------------------------\n08c | ", cohort, "\n", sep="")
  tryCatch({
    counts_path <- file.path(BULK_DIR, cohort, "bulk_counts_top300.rds")
    lib_path <- file.path(BULK_DIR, cohort, "library_sizes.tsv")
    old_path <- file.path(BASE, cohort, "main24", "top300", "proportions_wide.tsv")
    needed <- c(counts_path, lib_path, old_path)
    if (any(!file.exists(needed))) stop("Entrada necessaria ausente: ",
      paste(needed[!file.exists(needed)], collapse=", "))
    bulk <- as.matrix(readRDS(counts_path))
    storage.mode(bulk) <- "double"
    sample_ids <- colnames(bulk)
    if (!length(sample_ids) || anyNA(sample_ids) || anyDuplicated(sample_ids)) {
      stop("IDs bulk invalidos")
    }
    libs <- fread(lib_path)
    if (!all(c("sample_id", "library_size_full_input") %in% names(libs))) {
      stop("Arquivo library_sizes.tsv inesperado")
    }
    lib <- as.numeric(libs$library_size_full_input[match(sample_ids, libs$sample_id)])
    if (any(!is.finite(lib)) || any(lib <= 0) || any(lib < colSums(bulk)-1e-6)) {
      stop("Bibliotecas bulk invalidas")
    }
    common <- rownames(bulk)[rownames(bulk) %in% rownames(ref_pb)]
    b <- bulk[common, , drop=FALSE]
    r <- ref_pb[common, , drop=FALSE]
    keep <- rowSums(b) > 0 & as.numeric(Matrix::rowSums(r)) > 0
    genes <- common[keep]
    if (length(genes) < 500L) stop("Poucos genes apos filtro")
    b <- b[genes, , drop=FALSE]
    r <- r[genes, , drop=FALSE]
    y_bulk <- t(log2(sweep(b, 2, lib, "/") * 1e6 + 1))
    y_donor <- t(log2(sweep(as.matrix(r), 2, ref_lib, "/") * 1e6 + 1))
    y_pool <- t(log2(sweep(ref_pooled[genes, , drop=FALSE], 2,
                           ref_pooled_lib, "/") * 1e6 + 1))
    markers <- setNames(lapply(TYPES, function(ct) {
      as.integer(match(marker_tbl$gene_id[marker_tbl$cell_type == ct], genes,
                       nomatch=0L))
    }), TYPES)
    markers <- lapply(markers, function(x) x[x > 0L])
    if (any(lengths(markers) < 15L) || anyDuplicated(unlist(markers, use.names=FALSE))) {
      stop("Marcadores nao reproduzem atribuicao exclusiva do 08")
    }
    if (!identical(colnames(y_bulk), colnames(y_donor)) ||
        !identical(colnames(y_bulk), colnames(y_pool))) {
      stop("Desalinhamento dos genes")
    }
    # Fazer a reproducao PRIMEIRO; apenas depois executar as sensibilidades.
    run_one <- function(y, ref, pure) {
      z <- dtangle::dtangle(Y=y, references=ref, pure_samples=pure,
                            markers=markers, data_type="rna-seq")
      validate_estimates(z$estimates, sample_ids, cohort)
    }
    cat("Amostras:", length(sample_ids), " genes:", length(genes), "\n")
    est_base <- run_one(y_bulk, y_donor, pure_donor)
    baseline_file <- fread(old_path)
    if (!all(c("sample_id", TYPES) %in% names(baseline_file))) {
      stop("Resultado anterior 08 com colunas inesperadas")
    }
    ii <- match(sample_ids, baseline_file$sample_id)
    if (anyNA(ii) || length(ii) != nrow(baseline_file)) {
      stop("Amostras atuais diferem do baseline do 08")
    }
    old_est <- as.matrix(baseline_file[ii, ..TYPES])
    storage.mode(old_est) <- "double"
    err <- max(abs(est_base - old_est))
    cat("Diferenca maxima para o 08 original:", signif(err, 5), "\n")
    if (!is.finite(err) || err > 1e-6) {
      stop("NAO reproduziu a etapa 08; verificar versao, IDs e matrizes")
    }
    check_rows[[length(check_rows)+1L]] <- data.table(
      dataset=cohort, check="baseline_reproduced", status="PASS",
      observed=err, note="tolerancia 1e-6")

    # QN aplica a MESMA transformacao a referencias e bulks.
    qd <- qnorm_both(y_donor, y_bulk)
    qp <- qnorm_both(y_pool, y_bulk)
    all_est <- list(
      donor_raw=est_base,
      donor_quantile=run_one(qd$bulk, qd$ref, pure_donor),
      pooled_raw=run_one(y_bulk, y_pool, pure_pool),
      pooled_quantile=run_one(qp$bulk, qp$ref, pure_pool)
    )
    for (condition in names(all_est)) {
      e <- all_est[[condition]]
      long <- as.data.table(e, keep.rownames="sample_id")
      long <- melt(long, id.vars="sample_id", variable.name="cell_type",
                   value.name="proportion")
      long[, `:=`(dataset=cohort, condition=condition)]
      results[[length(results)+1L]] <- long[, .(dataset, condition, sample_id,
                                                cell_type, proportion)]
      bytype <- long[, .(
        n_samples=.N, mean_proportion=mean(proportion),
        sd_pp=100*stats::sd(proportion),
        fraction_zero=mean(proportion == 0)
      ), by=cell_type]
      # MAE individual calculado diretamente sobre as matrizes pareadas.
      mae_vec <- colMeans(abs(e - est_base))*100
      bytype[, mean_abs_change_from_baseline_pp := as.numeric(mae_vec[match(cell_type, TYPES)])]
      bytype[, `:=`(dataset=cohort, condition=condition)]
      ct_rows[[length(ct_rows)+1L]] <- bytype
      vascular <- mean(rowSums(e[, c("End", "Per"), drop=FALSE]))*100
      cat(sprintf("  %-18s End+Per=%.2f%% In=%.2f%% Opc=%.2f%%\n",
                  condition, vascular, 100*mean(e[,"In"]), 100*mean(e[,"Opc"])))
    }
    status <- data.table(dataset=cohort, check="all_4_conditions", status="PASS",
                         observed=4L, note=paste(length(genes), "genes"))
    check_rows[[length(check_rows)+1L]] <- status
  }, error=function(e) {
    cat("FALHA em ", cohort, ": ", conditionMessage(e), "\n", sep="")
    failures[[length(failures)+1L]] <<- data.table(dataset=cohort, error=conditionMessage(e))
    check_rows[[length(check_rows)+1L]] <<- data.table(
      dataset=cohort, check="run", status="FAIL", observed=NA_real_,
      note=conditionMessage(e))
  })
}

if (length(results)) {
  long_all <- rbindlist(results)
  stats_all <- rbindlist(ct_rows)
  fwrite(long_all, file.path(OUT, "08c_proportions_long_all.tsv"), sep="\t")
  fwrite(stats_all, file.path(OUT, "08c_celltype_summary.tsv"), sep="\t")
  vas <- stats_all[cell_type %in% c("End", "Per"), .(
    vascular_mean_percent=100*sum(mean_proportion)),
    by=.(dataset, condition)]
  neuron <- dcast(stats_all[cell_type %in% c("In", "Opc")],
                  dataset+condition~cell_type, value.var="mean_proportion")
  vas <- merge(vas[, .(dataset, condition, vascular_mean_percent)],
               neuron, by=c("dataset", "condition"), all.x=TRUE)
  setnames(vas, c("In", "Opc"), c("In_mean_fraction", "Opc_mean_fraction"))
  vas[, `:=`(In_mean_percent=100*In_mean_fraction,
             Opc_mean_percent=100*Opc_mean_fraction)]
  vas[, c("In_mean_fraction", "Opc_mean_fraction") := NULL]
  fwrite(vas, file.path(OUT, "08c_sensitivity_summary.tsv"), sep="\t")
  p <- ggplot(vas, aes(x=condition, y=vascular_mean_percent, fill=condition)) +
    geom_col(width=.65) + facet_wrap(~dataset, nrow=2) +
    scale_y_continuous(limits=c(0,100), expand=expansion(mult=c(0,.06))) +
    labs(title="dtangle: sensibilidade da estimativa End + Per",
         subtitle="Quatro condicoes metodologicas; nenhuma constitui validacao externa",
         x="Configuracao", y="End + Per estimado (%)") +
    theme_minimal(base_size=12) +
    theme(axis.text.x=element_text(angle=30,hjust=1), legend.position="none")
  tryCatch(ggsave(file.path(OUT, "08c_vascular_by_condition.png"), p,
                 width=12,height=8,dpi=180),
           error=function(e) warning("Falha no grafico: ",conditionMessage(e)))
}
checks <- if (length(check_rows)) rbindlist(check_rows,fill=TRUE) else data.table()
fwrite(checks, file.path(OUT,"08c_status.tsv"), sep="\t")
if (length(failures)) {
  fwrite(rbindlist(failures), file.path(OUT,"08c_errors.tsv"), sep="\t")
  stop("08c terminou com falhas; verifique 08c_status.tsv")
}
cat("\n08c concluido; outputs em ", OUT, "\n", sep="")
cat("Verifique 08c_sensitivity_summary.tsv, 08c_celltype_summary.tsv e 08c_status.tsv\n")
