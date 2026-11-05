#!/usr/bin/env Rscript
# ==============================================================================
# 08d | VALIDACAO SINTETICA DAS 4 CONFIGURACOES DO dtangle
# Estudo AD / PD / HD | Referencia Mathys (2019), pathology-free
# ==============================================================================
# OBJETIVO
# Avaliar a capacidade do dtangle de recuperar misturas de proporcoes
# CONHECIDAS produzidas a partir de um doador inteiramente excluido da
# referencia de treinamento (leave-one-donor-out, LODO).
#
# CONFIGURACOES (as mesmas escolhas metodologicas do 08c):
#   donor_raw        : perfis donor x tipo; log2(CPM+1)
#   donor_quantile   : perfis donor x tipo; log2(CPM+1) + quantilica conjunta
#   pooled_raw       : 1 perfil agregado por tipo; log2(CPM+1)
#   pooled_quantile  : 1 perfil agregado por tipo; + quantilica conjunta
#
# Duas atribuicoes de marcadores:
#   train_only : PRIMARIA: ranking top300 refeito usando APENAS treino
#                (escala de especificidade/deteccao inspirada no passo 04).
#   fixed_08   : SENSIBILIDADE: atribuicoes do 08 original. ATENCAO:
#                marcadores calculados originalmente com todos os 24 doadores,
#                incluindo holdout: vazamento na selecao de features.
#
# VERDADES CONHECIDAS
# A: fracao de NUCLEOS VIRTUAIS em cada mistura (alvo de desenho).
# B: fracao aproximada de RNA / reads por tipo, derivada do perfil do doador.
# As fracoes A e B sao diferentes quando os tipos tem RNA por nucleo desigual.
# A mistura e uma SOMA DETERMINISTICA de expressoes medias por nucleo: nao e
# bulk fisico, nem uma reamostragem independente de 1000 nucleos reais.
# End/Per do doador de teste podem vir de 1-2 nucleos observados.
#
# RESULTADOS NAO DEVEM SER USADOS para escolher metodo com base nas
# associacoes com doenca: aqui nao ha rotulos clinicos nem teste de grupos.
# Sintetico Mathys -> Mathys NAO e validacao externa bulk x snRNA.
#
# EXECUCAO NA RAIZ DO REPOSITORIO:
# Rscript scripts/phase02_deconvolution/08d_validacao_sintetica_dtangle.R
# DTANGLE08D_DONORS=2 Rscript scripts/phase02_deconvolution/08d_validacao_sintetica_dtangle.R
# DTANGLE08D_SKIP_FIXED=1 Rscript scripts/phase02_deconvolution/08d_validacao_sintetica_dtangle.R
#
# NAO instala pacotes, NAO atualiza renv.lock, NAO altera os resultados 08/08c.
# ==============================================================================

options(stringsAsFactors = FALSE)
PKGS <- c("data.table", "Matrix", "dtangle", "limma", "ggplot2")
missing_pkg <- PKGS[!vapply(PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkg)) stop("Pacotes ausentes: ", paste(missing_pkg, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
CONDITIONS <- c("donor_raw", "donor_quantile", "pooled_raw", "pooled_quantile")
SEED <- 62026L
N_HOLDOUT <- as.integer(Sys.getenv("DTANGLE08D_DONORS", "4"))
if (is.na(N_HOLDOUT) || N_HOLDOUT < 1L || N_HOLDOUT > 12L) {
  stop("DTANGLE08D_DONORS deve estar entre 1 e 12")
}
INCLUDE_FIXED <- Sys.getenv("DTANGLE08D_SKIP_FIXED", "0") != "1"
VIRTUAL_NUCLEI <- 1000L
N_MARKERS <- 300L
CPM_DETECT <- 1
SUM_TOL <- 1e-6
REFDIR <- "data/reference/mathys2019/phase02_prepared/reference"
OLD08 <- "results/deconvolution/dtangle"
OUT <- file.path(OLD08, "08d_synthetic_validation")
PLOTS <- file.path(OUT, "plots")
dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)

REF_COUNTS <- file.path(REFDIR, "mathys2019_pathology_free_geneid_counts_sparse.rds")
REF_META <- file.path(REFDIR, "mathys2019_pathology_free_cell_metadata_aligned.tsv")
FIXED_MARKER_FILE <- file.path(OLD08, "08_markers_unique_assignment.tsv")
for (p in c(REF_COUNTS, REF_META)) {
  if (!file.exists(p)) stop("Arquivo de referencia ausente: ", p)
}
if (INCLUDE_FIXED && !file.exists(FIXED_MARKER_FILE)) {
  stop("Marcadores do 08 original nao encontrados: ", FIXED_MARKER_FILE,
       "\nUse DTANGLE08D_SKIP_FIXED=1 para executar apenas train_only.")
}

# Matrizes: 9 misturas por doador, cada uma soma 1000 nucleos virtuais.
COMPOSITIONS <- rbind(
  equal = c(Ex=.125, In=.125, Ast=.125, Oli=.125, Opc=.125, Mic=.125, End=.125, Per=.125),
  cortical_like = c(Ex=.38, In=.17, Ast=.12, Oli=.20, Opc=.05, Mic=.04, End=.02, Per=.02),
  inhibitory_enriched = c(Ex=.22, In=.30, Ast=.12, Oli=.20, Opc=.06, Mic=.05, End=.03, Per=.02),
  opc_enriched = c(Ex=.20, In=.12, Ast=.13, Oli=.22, Opc=.20, Mic=.07, End=.03, Per=.03),
  vascular_enriched = c(Ex=.23, In=.12, Ast=.15, Oli=.18, Opc=.06, Mic=.06, End=.10, Per=.10),
  low_in = c(Ex=.45, In=.01, Ast=.13, Oli=.25, Opc=.06, Mic=.06, End=.02, Per=.02),
  low_opc = c(Ex=.38, In=.16, Ast=.14, Oli=.22, Opc=.01, Mic=.05, End=.02, Per=.02),
  low_vascular = c(Ex=.43, In=.17, Ast=.14, Oli=.19, Opc=.04, Mic=.02, End=.005, Per=.005),
  low_all_rare = c(Ex=.51, In=.01, Ast=.16, Oli=.28, Opc=.01, Mic=.02, End=.005, Per=.005)
)
if (!identical(colnames(COMPOSITIONS), TYPES) ||
    any(!is.finite(COMPOSITIONS)) || any(COMPOSITIONS <= 0) ||
    max(abs(rowSums(COMPOSITIONS)-1)) > 1e-10) {
  stop("Composicoes sinteticas invalidas")
}
VIRTUAL_COUNTS <- COMPOSITIONS * VIRTUAL_NUCLEI
if (any(abs(VIRTUAL_COUNTS - round(VIRTUAL_COUNTS)) > 1e-8)) {
  stop("Contagens virtuais devem ser inteiras")
}

cat("======================================================\n")
cat("08d | Validacao sintetica LODO das condicoes dtangle\n")
cat("======================================================\n")

sc <- readRDS(REF_COUNTS)
if (!inherits(sc, "dgCMatrix")) sc <- as(sc, "dgCMatrix")
if (is.null(rownames(sc)) || is.null(colnames(sc)) ||
    anyDuplicated(rownames(sc)) || anyDuplicated(colnames(sc)) ||
    nrow(sc) != 17394L || ncol(sc) != 35140L) {
  stop("Matriz Mathys fora do esperado 17394 x 35140, ou GeneIDs/TAG invalidos")
}
meta <- fread(REF_META)
if (!all(c("TAG", "donor_id", "cell_type") %in% names(meta))) {
  stop("Metadata sem TAG / donor_id / cell_type")
}
ix <- match(colnames(sc), as.character(meta[["TAG"]]))
if (anyNA(ix) || anyDuplicated(meta[["TAG"]])) stop("TAG sem correspondencia unica")
meta <- meta[ix]
if (!identical(as.character(meta[["TAG"]]), colnames(sc))) stop("Metadata desalinhado")
meta[, `:=`(donor_id=as.character(donor_id), cell_type=as.character(cell_type))]
if (uniqueN(meta$donor_id) != 24L || !setequal(unique(meta$cell_type), TYPES)) {
  stop("Referência nao tem 24 doadores / 8 tipos esperados")
}

# Recriar suporte por doador sem inventar combinacoes nao observadas.
support <- meta[, .N, by=.(donor_id,cell_type)]
support_wide <- dcast(support, donor_id ~ cell_type, value.var="N", fill=0L)
eligible <- support_wide[
  rowSums(as.matrix(support_wide[, ..TYPES]) > 0L) == length(TYPES)
]
if (nrow(eligible) < N_HOLDOUT) stop("Poucos doadores com os oito tipos observados")
set.seed(SEED)
chosen <- eligible[sample(seq_len(nrow(eligible)), N_HOLDOUT, replace=FALSE)]
chosen[, holdout_index := seq_len(.N)]
chosen[, holdout_label := sprintf("holdout_%02d", holdout_index)]
# Privacidade: exportar apenas rotulos locais, nao o identificador original.
fwrite(chosen[, c("holdout_label", TYPES), with=FALSE],
       file.path(OUT,"08d_donor_selection.tsv"), sep="\t")

# Agregacao unica dos 35140 nucleos nos grupos observados donor x tipo.
groups <- unique(meta[, .(donor_id,cell_type)])
setorder(groups, donor_id, cell_type)
full_group_id <- match(paste(meta$donor_id, meta$cell_type, sep="\r"),
                       paste(groups$donor_id, groups$cell_type, sep="\r"))
if (anyNA(full_group_id)) stop("Grupo de pseudobulk ausente")
A <- Matrix::sparseMatrix(i=seq_len(ncol(sc)), j=full_group_id,
                           x=rep(1, ncol(sc)),
                           dims=c(ncol(sc), nrow(groups)))
ref_pb_all <- sc %*% A
rownames(ref_pb_all) <- rownames(sc)
if (ncol(ref_pb_all) != 181L) stop("Esperados 181 grupos de referencia")
lib_pb_all <- as.numeric(Matrix::colSums(ref_pb_all))
if (any(lib_pb_all <= 0)) stop("Pseudobulks vazios")

fixed_markers <- NULL
if (INCLUDE_FIXED) {
  fixed_markers <- fread(FIXED_MARKER_FILE)
  if (!all(c("gene_id", "cell_type") %in% names(fixed_markers))) {
    stop("08_markers_unique_assignment.tsv sem gene_id / cell_type")
  }
  fixed_markers <- unique(fixed_markers[, .(
    gene_id=as.character(gene_id), cell_type=as.character(cell_type))])
  if (!setequal(unique(fixed_markers$cell_type), TYPES) ||
      anyDuplicated(fixed_markers$gene_id)) stop("Marcadores originais nao exclusivos")
}

# ------------------------------------------------------------------------------
# FUNCOES DE AUDITORIA / ESCALAS
# ------------------------------------------------------------------------------
log2cpm <- function(counts, lib) {
  if (length(lib) != ncol(counts) || any(!is.finite(lib)) || any(lib <= 0)) {
    stop("Library size invalida")
  }
  m <- as.matrix(counts)
  storage.mode(m) <- "double"
  log2(sweep(m, 2, lib, FUN="/") * 1e6 + 1)
}

qnorm_both <- function(y_ref, y_mix) {
  if (!identical(colnames(y_ref), colnames(y_mix))) {
    stop("Genes desalinhados antes da normalizacao quantilica")
  }
  nr <- nrow(y_ref)
  # limma espera genes x amostras
  joint <- cbind(t(y_ref), t(y_mix))
  q <- limma::normalizeBetweenArrays(joint, method="quantile")
  list(ref=t(q[,seq_len(nr),drop=FALSE]),
       bulk=t(q[,nr + seq_len(nrow(y_mix)),drop=FALSE]))
}

validate_estimates <- function(raw, sample_ids, label) {
  e <- as.matrix(raw)
  if (nrow(e) != length(sample_ids) || ncol(e) != 8L) {
    stop(label, ": dimensoes inesperadas")
  }
  if (is.null(colnames(e))) {
    # Na API do dtangle a ordem corresponde a pure_samples; algumas versoes
    # nao devolvem nomes. Aceita-se SOMENTE na ordem predefinida TYPES.
    colnames(e) <- TYPES
  } else if (!setequal(colnames(e), TYPES)) {
    stop(label, ": colunas dos tipos nao identificaveis")
  }
  e <- e[, TYPES, drop=FALSE]
  if (!is.null(rownames(e))) {
    if (!setequal(rownames(e), sample_ids)) stop(label, ": amostras divergentes")
    e <- e[sample_ids, , drop=FALSE]
  } else rownames(e) <- sample_ids
  storage.mode(e) <- "double"
  if (any(!is.finite(e)) || min(e) < -1e-8 ||
      max(e) > 1 + 1e-8 || max(abs(rowSums(e)-1)) > SUM_TOL) {
    stop(label, ": estimativas invalidas")
  }
  e
}

# Seleciona top300 DENTRO de cada fold, SEM acessar o doador de teste.
# Metodo exploratorio inspirando-se no score usado na etapa 04:
# especificidade * log1p(media CPM) * fracao de doadores com CPM>1.
# Ranking deterministico por score, especificidade, expressao e GeneID.
make_training_markers <- function(pb_counts, pb_lib, pb_labels, gene_ids) {
  pcpm <- sweep(as.matrix(pb_counts), 2, pb_lib, FUN="/") * 1e6
  means <- matrix(0, nrow=nrow(pcpm), ncol=length(TYPES),
                  dimnames=list(gene_ids,TYPES))
  detect <- means
  for (k in seq_along(TYPES)) {
    id <- which(pb_labels == TYPES[k])
    if (!length(id)) stop("Nenhum pseudobulk de ",TYPES[k])
    means[,k] <- rowMeans(pcpm[,id,drop=FALSE])
    detect[,k] <- rowMeans(pcpm[,id,drop=FALSE] > CPM_DETECT)
  }
  denom <- pmax(rowSums(means), 1e-12)
  out <- vector("list",length(TYPES))
  for (k in seq_along(TYPES)) {
    spec <- means[,k]/denom
    score <- spec*log1p(means[,k])*detect[,k]
    # Todos os genes com escore zero ficam inelegiveis.
    ids <- which(is.finite(score) & score > 0)
    if (length(ids) < N_MARKERS) stop("Menos de 300 genes elegiveis em ",TYPES[k])
    ord <- order(-score[ids], -spec[ids], -means[ids,k], gene_ids[ids])
    take <- ids[ord[seq_len(N_MARKERS)]]
    out[[k]] <- data.table(gene_id=gene_ids[take],cell_type=TYPES[k],
                            rank=seq_len(N_MARKERS),score=score[take],
                            ct_index=k)
  }
  candidates <- rbindlist(out)
  setorder(candidates, rank, -score, ct_index, gene_id)
  unique_assignment <- candidates[!duplicated(gene_id)]
  if (uniqueN(unique_assignment$cell_type) != 8L ||
      anyDuplicated(unique_assignment$gene_id)) stop("Ranking nao gerou 8 tipos exclusivos")
  unique_assignment[, c("rank","score","ct_index") := NULL]
  unique_assignment
}

make_heldout_mixtures <- function(counts_full, meta_full, donor) {
  means <- matrix(0, nrow=nrow(counts_full), ncol=8L,
                  dimnames=list(rownames(counts_full),TYPES))
  for (k in seq_along(TYPES)) {
    jj <- which(meta_full$donor_id == donor & meta_full$cell_type == TYPES[k])
    if (!length(jj)) stop("Doador de teste sem ",TYPES[k])
    means[,k] <- as.numeric(Matrix::rowMeans(counts_full[,jj,drop=FALSE]))
  }
  mix_unrounded <- means %*% t(VIRTUAL_COUNTS)
  mix <- round(mix_unrounded)
  colnames(mix) <- rownames(COMPOSITIONS)
  # Verdade de FRACAO DE RNA capturado: diferenca importante em relacao
  # a fracao de nucleos, dada a variabilidade de contagens por nucleo.
  cell_rna <- colSums(means)
  contrib <- sweep(VIRTUAL_COUNTS, 2, cell_rna, FUN="*")
  true_rna <- sweep(contrib, 1, rowSums(contrib), FUN="/")
  list(counts=mix, true_rna=true_rna, mean_library_by_type=cell_rna)
}

# ------------------------------------------------------------------------------
# AJUSTES, UM DOADOR EXCLUIDO POR FOLD
# ------------------------------------------------------------------------------
predictions <- list()
run_audit <- list()
marker_audit <- list()
truth_tables <- list()
checks <- list()
errors <- list()

for (holdout_i in seq_len(nrow(chosen))) {
  donor <- chosen$donor_id[holdout_i]
  holdout_label <- chosen$holdout_label[holdout_i]
  cat("\n------------------------------------------------------\n")
  cat("08d | ",holdout_label," (23 doadores de treino)\n",sep="")
  tryCatch({
    train_g <- which(groups$donor_id != donor)
    test_g <- which(groups$donor_id == donor)
    if (length(intersect(train_g,test_g)) ||
        uniqueN(groups$donor_id[train_g]) != 23L ||
        !setequal(groups$cell_type[train_g], TYPES)) {
      stop("Falha de separacao dos doadores / tipos")
    }
    g_train <- groups[train_g]
    pb_train <- ref_pb_all[,train_g,drop=FALSE]
    lib_train <- lib_pb_all[train_g]
    mat_pooled <- sapply(TYPES,function(ct) {
      as.numeric(Matrix::rowSums(pb_train[,g_train$cell_type == ct,drop=FALSE]))
    })
    rownames(mat_pooled) <- rownames(ref_pb_all)
    lib_pooled <- colSums(mat_pooled)
    if (any(lib_pooled <= 0)) stop("Pooled vazio")

    synth <- make_heldout_mixtures(sc,meta,donor)
    synthetic <- synth$counts
    sample_ids <- paste(holdout_label, colnames(synthetic), sep="_")
    colnames(synthetic) <- sample_ids
    lib_synth <- colSums(synthetic)
    if (any(lib_synth <= 0)) stop("Bibliotecas sinteticas vazias")
    truth_rna <- synth$true_rna
    rownames(truth_rna) <- sample_ids
    truth_nuc <- COMPOSITIONS
    rownames(truth_nuc) <- sample_ids
    truth_tbl <- melt(as.data.table(truth_nuc,keep.rownames="sample_id"),
                      id.vars="sample_id",variable.name="cell_type",
                      value.name="true_nucleus_fraction")
    rna_tbl <- melt(as.data.table(truth_rna,keep.rownames="sample_id"),
                    id.vars="sample_id",variable.name="cell_type",
                    value.name="true_rna_fraction")
    truth_tbl <- merge(truth_tbl,rna_tbl,by=c("sample_id","cell_type"),sort=FALSE)
    truth_tbl[, `:=`(holdout=holdout_label,
                      scenario=sub(paste0("^", holdout_label, "_"), "", sample_id))]
    if (!setequal(unique(truth_tbl$scenario), rownames(COMPOSITIONS)) ||
        nrow(truth_tbl) != nrow(COMPOSITIONS)*length(TYPES)) {
      stop("Desenho das misturas ficou desalinhado")
    }
    truth_tables[[length(truth_tables)+1L]] <- truth_tbl

    variants <- list(train_only=make_training_markers(
      pb_train,lib_train,g_train$cell_type,rownames(pb_train)))
    if (INCLUDE_FIXED) variants$fixed_08 <- fixed_markers

    for (marker_source in names(variants)) {
      mt <- variants[[marker_source]]
      genes <- rownames(ref_pb_all)[rownames(ref_pb_all) %in% mt$gene_id]
      keep <- rowSums(synthetic[genes,,drop=FALSE]) > 0 &
              as.numeric(Matrix::rowSums(pb_train[genes,,drop=FALSE])) > 0
      genes <- genes[keep]
      if (length(genes) < 500L) stop("Menos de 500 marcadores sobreviventes")
      markers <- setNames(lapply(TYPES,function(ct) {
        pos <- as.integer(match(mt$gene_id[mt$cell_type == ct],genes,nomatch=0L))
        pos[pos > 0L]
      }), TYPES)
      if (any(lengths(markers) < 15L) ||
          anyDuplicated(unlist(markers,use.names=FALSE))) {
        stop("Marcadores fora dos 8 tipos, insuficientes ou sobrepostos")
      }
      marker_audit[[length(marker_audit)+1L]] <- data.table(
        holdout=holdout_label,marker_source=marker_source,
        genes_used=length(genes),
        min_markers_per_type=min(lengths(markers)),
        max_markers_per_type=max(lengths(markers)),
        feature_leakage=(marker_source=="fixed_08"),
        training_donors=23L)
      y_bulk <- t(log2cpm(synthetic[genes,,drop=FALSE],lib_synth))
      y_donor <- t(log2cpm(pb_train[genes,,drop=FALSE],lib_train))
      y_pool <- t(log2cpm(mat_pooled[genes,,drop=FALSE],lib_pooled))
      if (!identical(colnames(y_bulk),colnames(y_donor)) ||
          !identical(colnames(y_bulk),colnames(y_pool))) stop("Genes desalinhados")
      pure_donor <- setNames(lapply(TYPES,function(ct) {
        which(g_train$cell_type == ct)
      }),TYPES)
      pure_pool <- setNames(lapply(seq_along(TYPES),function(i) as.integer(i)),TYPES)
      # Quantilica feita separadamente para a versao donor e pooled,
      # mas sempre incluindo as mesmas 9 misturas de teste do holdout.
      qd <- qnorm_both(y_donor,y_bulk)
      qp <- qnorm_both(y_pool,y_bulk)
      inputs <- list(
        donor_raw=list(bulk=y_bulk,ref=y_donor,pure=pure_donor),
        donor_quantile=list(bulk=qd$bulk,ref=qd$ref,pure=pure_donor),
        pooled_raw=list(bulk=y_bulk,ref=y_pool,pure=pure_pool),
        pooled_quantile=list(bulk=qp$bulk,ref=qp$ref,pure=pure_pool)
      )
      for (condition in CONDITIONS) {
        inp <- inputs[[condition]]
        t0 <- Sys.time()
        result <- dtangle::dtangle(Y=inp$bulk,references=inp$ref,
                                  pure_samples=inp$pure,markers=markers,
                                  data_type="rna-seq")
        seconds <- as.numeric(difftime(Sys.time(),t0,units="secs"))
        est <- validate_estimates(result$estimates,sample_ids,
                   paste(holdout_label,marker_source,condition))
        dat <- melt(as.data.table(est,keep.rownames="sample_id"),
                    id.vars="sample_id",variable.name="cell_type",
                    value.name="estimated_fraction")
        dat <- merge(dat,truth_tbl[,.(sample_id,cell_type,scenario,
                                    true_nucleus_fraction,true_rna_fraction)],
                     by=c("sample_id","cell_type"),sort=FALSE)
        dat[, `:=`(holdout=holdout_label,marker_source=marker_source,
                   condition=condition,
                   error_nucleus_pp=100*(estimated_fraction-true_nucleus_fraction),
                   abs_error_nucleus_pp=100*abs(estimated_fraction-true_nucleus_fraction),
                   error_rna_pp=100*(estimated_fraction-true_rna_fraction),
                   abs_error_rna_pp=100*abs(estimated_fraction-true_rna_fraction))]
        predictions[[length(predictions)+1L]] <- dat
        run_audit[[length(run_audit)+1L]] <- data.table(
          holdout=holdout_label,marker_source=marker_source,condition=condition,
          n_train_donors=23L,n_train_pseudobulks=length(train_g),
          n_train_nuclei=sum(meta$donor_id != donor),
          n_test_nuclei_observed=sum(meta$donor_id == donor),
          genes_used=length(genes),n_mixtures=nrow(est),
          max_sum_error=max(abs(rowSums(est)-1)),
          min_estimate=min(est),max_estimate=max(est),
          elapsed_seconds=seconds,status="PASS")
        cat(sprintf("  %-10s %-17s | MAE nucleos = %.2f pp, End+Per = %.2f%%\n",
                    marker_source,condition,mean(dat$abs_error_nucleus_pp),
                    100*mean(rowSums(est[,c("End","Per"),drop=FALSE]))))
      }
    }
    checks[[length(checks)+1L]] <- data.table(holdout=holdout_label,
        check="no_train_test_nucleus_overlap",status="PASS",value="TRUE")
    # Checkpoints em arquivo separado, sem tocar 08/08c.
    fwrite(rbindlist(predictions,fill=TRUE),file.path(OUT,"08d_predictions_partial.tsv"),sep="\t")
    fwrite(rbindlist(run_audit,fill=TRUE),file.path(OUT,"08d_run_audit_partial.tsv"),sep="\t")
  },error=function(e) {
    msg <- conditionMessage(e)
    message("FALHA: ",holdout_label,": ",msg)
    errors[[length(errors)+1L]] <<- data.table(holdout=holdout_label,error=msg)
    checks[[length(checks)+1L]] <<- data.table(holdout=holdout_label,
      check="fit",status="FAIL",value=msg)
  })
}

# ------------------------------------------------------------------------------
# RESULTADOS, METRICAS, GRAFICOS E AUDITORIA
# ------------------------------------------------------------------------------
status <- if (length(checks)) rbindlist(checks,fill=TRUE) else data.table()
fwrite(status,file.path(OUT,"08d_status.tsv"),sep="\t")
if (length(errors)) fwrite(rbindlist(errors),file.path(OUT,"08d_errors.tsv"),sep="\t")
if (!length(predictions)) stop("08d nao produziu estimativas; confira status")

pred <- rbindlist(predictions,fill=TRUE)
run <- rbindlist(run_audit,fill=TRUE)
marks <- rbindlist(marker_audit,fill=TRUE)
truth <- rbindlist(truth_tables,fill=TRUE)
fwrite(pred,file.path(OUT,"08d_predictions_long.tsv"),sep="\t")
fwrite(run,file.path(OUT,"08d_run_audit.tsv"),sep="\t")
fwrite(marks,file.path(OUT,"08d_marker_audit.tsv"),sep="\t")
fwrite(truth,file.path(OUT,"08d_truths.tsv"),sep="\t")
fwrite(as.data.table(COMPOSITIONS,keep.rownames="scenario"),
       file.path(OUT,"08d_mixture_design.tsv"),sep="\t")

metrics_type <- pred[,.(
  n_predictions=as.integer(.N),
  mae_nucleus_pp=as.numeric(mean(abs_error_nucleus_pp)),
  rmse_nucleus_pp=as.numeric(sqrt(mean(error_nucleus_pp^2))),
  bias_nucleus_pp=as.numeric(mean(error_nucleus_pp)),
  mae_rna_pp=as.numeric(mean(abs_error_rna_pp)),
  bias_rna_pp=as.numeric(mean(error_rna_pp)),
  zero_for_positive_truth=as.numeric(mean(estimated_fraction <= 1e-10 &
                                             true_nucleus_fraction > 0)),
  mean_truth_nucleus=as.numeric(mean(true_nucleus_fraction)),
  mean_estimate=as.numeric(mean(estimated_fraction))
),by=.(marker_source,condition,cell_type)]

metrics_cond <- pred[,.(
  n_predictions=as.integer(.N),
  n_holdouts=as.integer(uniqueN(holdout)),
  n_scenarios=as.integer(uniqueN(scenario)),
  mae_nucleus_pp=as.numeric(mean(abs_error_nucleus_pp)),
  rmse_nucleus_pp=as.numeric(sqrt(mean(error_nucleus_pp^2))),
  mae_rna_pp=as.numeric(mean(abs_error_rna_pp)),
  mean_abs_bias_in_pp=as.numeric(abs(mean(error_nucleus_pp[cell_type == "In"]))),
  mean_abs_bias_opc_pp=as.numeric(abs(mean(error_nucleus_pp[cell_type == "Opc"]))),
  mean_abs_bias_vascular_pp=as.numeric(abs(mean(error_nucleus_pp[cell_type %in% c("End","Per")]))),
  fraction_exact_zero=as.numeric(mean(estimated_fraction <= 1e-10))
),by=.(marker_source,condition)]

# Cada holdout tem os mesmos 9 cenarios; metricas por doador deixam explicita
# a unidade biologica independente (doador), e nao tratam as 9 misturas
# correlacionadas como 9 doadores distintos.
metrics_holdout <- pred[,.(
  mae_nucleus_pp=as.numeric(mean(abs_error_nucleus_pp)),
  mae_rna_pp=as.numeric(mean(abs_error_rna_pp)),
  mae_vascular_pp=as.numeric(mean(abs_error_nucleus_pp[cell_type %in% c("End","Per")])),
  mae_in_pp=as.numeric(mean(abs_error_nucleus_pp[cell_type == "In"])),
  mae_opc_pp=as.numeric(mean(abs_error_nucleus_pp[cell_type == "Opc"]))
),by=.(holdout,marker_source,condition)]

metrics_scenario <- pred[,.(
  mae_nucleus_pp=as.numeric(mean(abs_error_nucleus_pp)),
  mae_rna_pp=as.numeric(mean(abs_error_rna_pp)),
  mae_in_pp=as.numeric(mean(abs_error_nucleus_pp[cell_type == "In"])),
  mae_opc_pp=as.numeric(mean(abs_error_nucleus_pp[cell_type == "Opc"])),
  mae_vascular_pp=as.numeric(mean(abs_error_nucleus_pp[cell_type %in% c("End","Per")]))
),by=.(marker_source,condition,scenario)]

fwrite(metrics_type,file.path(OUT,"08d_metrics_by_type.tsv"),sep="\t")
fwrite(metrics_cond,file.path(OUT,"08d_metrics_by_condition.tsv"),sep="\t")
fwrite(metrics_holdout,file.path(OUT,"08d_metrics_by_holdout.tsv"),sep="\t")
fwrite(metrics_scenario,file.path(OUT,"08d_metrics_by_scenario.tsv"),sep="\t")

# Diferencas pareadas (mesma mistura / mesmo doador) versus donor_raw.
baseline <- pred[condition=="donor_raw",.(holdout,marker_source,sample_id,cell_type,
                                          baseline_abs=abs_error_nucleus_pp)]
paired <- merge(pred,baseline,by=c("holdout","marker_source","sample_id","cell_type"),
                all.x=TRUE,sort=FALSE)
paired[, improvement_pp := baseline_abs-abs_error_nucleus_pp]
paired_summary <- paired[,.(
   mean_improvement_pp=as.numeric(mean(improvement_pp)),
   fraction_improved=as.numeric(mean(improvement_pp > 0)),
   n_compared=as.integer(.N)
),by=.(marker_source,condition)]
fwrite(paired_summary,file.path(OUT,"08d_paired_improvements.tsv"),sep="\t")

tryCatch({
  # Principal: apenas marcadores selecionados nos doadores de treino.
  primary <- metrics_cond[marker_source=="train_only"]
  p1 <- ggplot(primary,aes(x=condition,y=mae_nucleus_pp,fill=condition))+
    geom_col(width=.7)+theme_minimal(base_size=12)+
    labs(title="08d - erro absoluto medio (validacao por doador)",
         subtitle="Marcadores reranqueados sem doador de teste; misturas sinteticas da Mathys",
         x="Configuracao dtangle",y="MAE da fracao de nucleos (pp)")+
    theme(axis.text.x=element_text(angle=30,hjust=1),legend.position="none")
  ggsave(file.path(PLOTS,"08d_mae_by_condition.png"),p1,width=11,height=7,dpi=170)
  p2 <- ggplot(metrics_type[marker_source=="train_only"],
               aes(x=cell_type,y=mae_nucleus_pp,fill=condition))+
    geom_col(position=position_dodge(width=.8),width=.7)+
    theme_minimal(base_size=11)+
    labs(title="08d - erro medio por tipo celular",x="Tipo celular",
         y="MAE nuclear (pp)",fill="Configuracao")
  ggsave(file.path(PLOTS,"08d_mae_by_type.png"),p2,width=12,height=7,dpi=170)
  p3 <- ggplot(pred[marker_source=="train_only"],
               aes(x=100*true_nucleus_fraction,y=100*estimated_fraction,
                   color=condition))+
    geom_point(alpha=.40,size=.75)+geom_abline(slope=1,intercept=0,linetype=2)+
    facet_wrap(~cell_type,scales="free",ncol=4)+
    theme_minimal(base_size=11)+
    labs(title="08d - verdadeiro vs estimado, por tipo",
         subtitle="Pontos sao misturas vinculadas a quatro doadores; linha = recuperacao perfeita",
         x="Fracao virtual verdadeira (%)",y="Fracao estimada (%)")
  ggsave(file.path(PLOTS,"08d_true_vs_estimated.png"),p3,width=15,height=9,dpi=170)
},error=function(e) warning("08d: falha apenas nos graficos: ",conditionMessage(e)))

writeLines(capture.output(sessionInfo()),file.path(OUT,"08d_sessionInfo.txt"))
cat("\nResumo global (marcadores treino, MAE pp):\n")
print(metrics_cond[marker_source=="train_only",
                   .(condition,n_holdouts,n_scenarios,mae_nucleus_pp,rmse_nucleus_pp,mae_rna_pp)])
cat("\nSaidas em: ",OUT,"\n",sep="")
if (length(errors)) stop("08d teve falhas parciais: verifique 08d_errors.tsv")
if (uniqueN(run$holdout) != N_HOLDOUT ||
    any(run[, .N, by=.(holdout,marker_source)]$N != 4L)) {
  stop("Nem todos os 4 ajustes por doador / marcadores foram concluidos")
}
cat("08d CONCLUIDO: sintese exploratoria; nao equivale a validacao externa.\n")
