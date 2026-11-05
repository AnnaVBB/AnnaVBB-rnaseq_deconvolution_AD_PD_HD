#!/usr/bin/env Rscript
# ==============================================================================
# 08f | VALIDACAO SINTETICA LODO DO DTANGLE NA CONFIGURACAO 08e
# Projeto: AD / PD / HD | Mathys 2019 pathology-free | oito tipos celulares
# ==============================================================================
# PERGUNTA
# A configuracao 08e (todos os genes compartilhados + 8 perfis por tipo +
# 100 marcadores automaticos) recupera proporcoes CONHECIDAS de misturas
# sinteticas produzidas a partir de doadores NAO usados no treinamento?
#
# CONFIGURACOES ANALISADAS (sem normalizacao quantilica):
#   pooled_counts     = soma dos counts por tipo, CPM do total, log2(CPM+1)
#   mean_nucleus_cpm  = media dos CPMs individuais por tipo, log2(CPM+1)
#
# PARAMETROS IDENTICOS AO 08e:
#   dtangle(Y=..., references=..., pure_samples=..., markers=NULL,
#           n_markers=100, marker_method='ratio', data_type='rna-seq')
#   Universo de genes = genes da referencia expressos em treino e na mistura.
#   Referencia = 8 perfis (Ex/In/Ast/Oli/Opc/Mic/End/Per).
#
# LEAVE-ONE-DONOR-OUT (LODO):
#   - Selecionar quatro doadores com todos os oito tipos observados.
#   - Para cada fold, excluir completamente o doador de teste da referencia.
#   - Montar as nove misturas a partir das contagens do doador excluido.
#   - Recalcular os dois perfis e os marcadores somente com os outros 23.
#   - Auditar a selecao automatica versus ranking so da referencia de treino.
#
# IMPORTANTES LIMITACOES:
# * As misturas sao somas deterministicas de expressoes medias de nucleos,
#   nao bulk fisico preparado em laboratorio.
# * End e Per por doador podem ter apenas 1-2 nucleos; verdade sintetica
#   nesses tipos possui elevada incerteza biologica.
# * 'truth_nucleus' = proporcao de nucleos virtuais no desenho.
# * 'truth_rna' = contribuicao teorica de CONTAGENS/reads segundo bibliotecas
#   do doador. NAO equivale a massa de RNA validada experimentalmente.
# * A normalizacao quantilica NAO e usada: ela nao consta da configuracao 08e.
# * Resultados nao sao validacao externa de bulk nem inferencia sobre doenca.
# * NAO modifica 08, 08c, 08d, 08e, 07 ou arquivos renv.
#
# EXECUCAO (a partir da raiz do repositorio):
# Rscript scripts/phase02_deconvolution/08f_validacao_sintetica_dtangle_08e.R
#
# Teste rapido (apenas um doador; resultados INCOMPLETOS):
# DTANGLE08F_DONORS=1 Rscript scripts/phase02_deconvolution/08f_validacao_sintetica_dtangle_08e.R
#
# Se preferir testar somente uma construcao de referencia:
# DTANGLE08F_MODE=pooled_counts Rscript scripts/phase02_deconvolution/08f_validacao_sintetica_dtangle_08e.R
#
# Saidas: results/deconvolution/dtangle/08f_synthetic_validation/
# ==============================================================================

options(stringsAsFactors=FALSE)
REQUIRED <- c('data.table','Matrix','dtangle','ggplot2')
MISSING <- REQUIRED[!vapply(REQUIRED, requireNamespace, logical(1), quietly=TRUE)]
if (length(MISSING)) stop('Pacotes ausentes: ',paste(MISSING,collapse=', '),
                         '. Instale no ambiente apropriado sem alterar renv automaticamente.')
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

TYPES <- c('Ex','In','Ast','Oli','Opc','Mic','End','Per')
REF_MODES <- c('pooled_counts','mean_nucleus_cpm')
N_MARKERS <- 100L
MARKER_METHOD <- 'ratio'
DATA_TYPE <- 'rna-seq'
SEED <- 62026L                 # Igual ao 08d, para selecionar os mesmos doadores.
VIRTUAL_NUCLEI <- 1000L
SUM_TOL <- 1e-6
N_HOLDOUT <- suppressWarnings(as.integer(Sys.getenv('DTANGLE08F_DONORS','4')))
if (length(N_HOLDOUT)!=1L || is.na(N_HOLDOUT) || N_HOLDOUT<1L || N_HOLDOUT>12L)
  stop('DTANGLE08F_DONORS deve ser inteiro de 1 a 12')
MODE_ENV <- trimws(Sys.getenv('DTANGLE08F_MODE',''))
MODES <- if (nzchar(MODE_ENV)) trimws(strsplit(MODE_ENV,',',fixed=TRUE)[[1L]]) else REF_MODES
if (!length(MODES) || any(!MODES %in% REF_MODES) || anyDuplicated(MODES))
  stop('DTANGLE08F_MODE invalido: pooled_counts,mean_nucleus_cpm')

REF_ROOT <- 'data/reference/mathys2019/phase02_prepared/reference'
REF_COUNTS <- file.path(REF_ROOT,'mathys2019_pathology_free_geneid_counts_sparse.rds')
REF_META <- file.path(REF_ROOT,'mathys2019_pathology_free_cell_metadata_aligned.tsv')
OUT <- 'results/deconvolution/dtangle/08f_synthetic_validation'
PLOTS <- file.path(OUT,'plots')
dir.create(PLOTS,recursive=TRUE,showWarnings=FALSE)
for (p in c(REF_COUNTS,REF_META)) if (!file.exists(p)) stop('Arquivo ausente: ',p)

# Mesmo desenho das nove composicoes do 08d para comparacao justa.
COMPOSITIONS <- rbind(
  equal               = c(Ex=.125,In=.125,Ast=.125,Oli=.125,Opc=.125,Mic=.125,End=.125,Per=.125),
  cortical_like       = c(Ex=.38,In=.17,Ast=.12,Oli=.20,Opc=.05,Mic=.04,End=.02,Per=.02),
  inhibitory_enriched = c(Ex=.22,In=.30,Ast=.12,Oli=.20,Opc=.06,Mic=.05,End=.03,Per=.02),
  opc_enriched        = c(Ex=.20,In=.12,Ast=.13,Oli=.22,Opc=.20,Mic=.07,End=.03,Per=.03),
  vascular_enriched   = c(Ex=.23,In=.12,Ast=.15,Oli=.18,Opc=.06,Mic=.06,End=.10,Per=.10),
  low_in              = c(Ex=.45,In=.01,Ast=.13,Oli=.25,Opc=.06,Mic=.06,End=.02,Per=.02),
  low_opc             = c(Ex=.38,In=.16,Ast=.14,Oli=.22,Opc=.01,Mic=.05,End=.02,Per=.02),
  low_vascular        = c(Ex=.43,In=.17,Ast=.14,Oli=.19,Opc=.04,Mic=.02,End=.005,Per=.005),
  low_all_rare        = c(Ex=.51,In=.01,Ast=.16,Oli=.28,Opc=.01,Mic=.02,End=.005,Per=.005)
)
if (!identical(colnames(COMPOSITIONS),TYPES) || any(COMPOSITIONS<=0) ||
    any(!is.finite(COMPOSITIONS)) || max(abs(rowSums(COMPOSITIONS)-1))>1e-10)
  stop('Composicoes sinteticas invalidas')
VIRTUAL_COUNTS <- COMPOSITIONS * VIRTUAL_NUCLEI
if (any(abs(VIRTUAL_COUNTS-round(VIRTUAL_COUNTS))>1e-8))
  stop('Misturas devem conter contagens virtuais inteiras')

write_table <- function(tab,filename) {
  data.table::fwrite(tab,file.path(OUT,filename),sep='\t',na='NA',quote=FALSE)
}
verify_ids <- function(ids,label) {
  if (is.null(ids) || !length(ids) || anyNA(ids) || any(!nzchar(as.character(ids))) ||
      anyDuplicated(ids)) stop(label,': identificadores ausentes/duplicados')
}
verify_prop <- function(est,sample_ids,tag) {
  x <- as.matrix(est)
  if (nrow(x)!=length(sample_ids) || ncol(x)!=length(TYPES))
    stop(tag,': dimensoes estimadas inesperadas')
  if (!is.null(rownames(x)) && !identical(rownames(x),sample_ids)) {
    if (!setequal(rownames(x),sample_ids)) stop(tag,': amostras nao correspondem')
    x <- x[sample_ids,,drop=FALSE]
  }
  if (!is.null(colnames(x)) && !identical(colnames(x),TYPES)) {
    if (!setequal(colnames(x),TYPES)) stop(tag,': tipos celulares nao correspondem')
    x <- x[,TYPES,drop=FALSE]
  }
  rownames(x) <- sample_ids
  colnames(x) <- TYPES
  storage.mode(x) <- 'double'
  if (any(!is.finite(x)) || any(x < -1e-8) || any(x>1+1e-8) ||
      max(abs(rowSums(x)-1))>SUM_TOL) stop(tag,': proporcoes invalidas')
  if (any(x<0)) stop(tag,': proporcoes negativas; sem truncagem silenciosa')
  x
}
make_truth <- function(pseudomean,sample_ids) {
  # Fração de contagens teóricas em todos os 17.394 genes anotados.
  per_type_reads <- colSums(pseudomean)
  if (any(per_type_reads<=0) || any(!is.finite(per_type_reads)))
    stop('Media de counts por nucleo invalida no doador excluido')
  weighted <- sweep(VIRTUAL_COUNTS,2,per_type_reads,'*')
  truth_rna <- sweep(weighted,1,rowSums(weighted),'/')
  rownames(truth_rna) <- sample_ids
  truth_nuc <- COMPOSITIONS
  rownames(truth_nuc) <- sample_ids
  list(nuc=truth_nuc,rna=truth_rna)
}

cat('======================================================\n')
cat('08f | Validacao sintética LODO da configuracao dtangle 08e\n')
cat('======================================================\n')
cat('Pacote dtangle: ',as.character(utils::packageVersion('dtangle')),'\n',sep='')
cat('Modos: ',paste(MODES,collapse=', '),' | doadores: ',N_HOLDOUT,'\n',sep='')

# ------------------------------------------------------------------------------
# 1. Leitura e integridade da referencia
# ------------------------------------------------------------------------------
sc <- readRDS(REF_COUNTS)
if (!inherits(sc,'dgCMatrix')) sc <- methods::as(sc,'dgCMatrix')
if (is.null(rownames(sc)) || is.null(colnames(sc))) stop('Mathys sem GeneIDs/TAG')
verify_ids(rownames(sc),'GeneIDs Mathys')
verify_ids(colnames(sc),'TAGs Mathys')
if (nrow(sc)!=17394L || ncol(sc)!=35140L)
  stop('Referência inesperada: esperado 17394 genes x 35140 nucleos')
if (any(!is.finite(sc@x)) || any(sc@x<0) || any(abs(sc@x-round(sc@x))>1e-7))
  stop('Counts da referencia invalidos')

meta <- data.table::fread(REF_META)
if (!all(c('TAG','donor_id','cell_type') %in% names(meta)))
  stop('Metadata sem TAG, donor_id ou cell_type')
meta[,`:=`(TAG=as.character(TAG),donor_id=as.character(donor_id),
           cell_type=as.character(cell_type))]
verify_ids(meta$TAG,'TAG metadata')
idx <- match(colnames(sc),meta$TAG)
if (anyNA(idx)) stop('TAG sem correspondencia')
meta <- meta[idx]
if (!identical(meta$TAG,colnames(sc)) || data.table::uniqueN(meta$donor_id)!=24L ||
    !setequal(unique(meta$cell_type),TYPES)) stop('Metadata desalinhado ou referencia distinta da main24')
lib_nuc <- as.numeric(Matrix::colSums(sc))
if (any(!is.finite(lib_nuc)) || any(lib_nuc<=0)) stop('Nucleo com biblioteca nula')

# Escolher os MESMOS quatro doadores utilizados no 08d, desde que
# a ordem/alinhamento dos arquivos de referencia seja a mesma.
support <- meta[, .N, by=.(donor_id,cell_type)]
support_wide <- data.table::dcast(support,donor_id ~ cell_type,value.var='N',fill=0L)
eligible <- support_wide[rowSums(as.matrix(support_wide[,..TYPES])>0L)==length(TYPES)]
if (nrow(eligible)<N_HOLDOUT) stop('Menos doadores elegiveis que o solicitado')
set.seed(SEED)
chosen <- eligible[sample(seq_len(nrow(eligible)),N_HOLDOUT,replace=FALSE)]
chosen[,fold_index:=seq_len(.N)]
chosen[,fold:=sprintf('holdout_%02d',fold_index)]
# Exportamos apenas nomes locais e suporte celular, nao IDs dos doadores.
write_table(chosen[,c('fold',TYPES),with=FALSE],'08f_holdout_support.tsv')
write_table(data.table::as.data.table(COMPOSITIONS,keep.rownames='scenario'),
            '08f_mixture_design.tsv')

# Teste opcional de comparabilidade com os doadores escolhidos pelo 08d.
old_sel <- 'results/deconvolution/dtangle/08d_synthetic_validation/08d_donor_selection.tsv'
comp_08d <- 'not_available'
if (file.exists(old_sel)) {
  prev <- data.table::fread(old_sel)
  now <- chosen[,c('fold',TYPES),with=FALSE]
  if (all(c('holdout_label',TYPES) %in% names(prev)) && nrow(prev)>=N_HOLDOUT) {
    now_df <- as.data.frame(now[,..TYPES]); prev_df <- as.data.frame(prev[seq_len(N_HOLDOUT),..TYPES])
    comp_08d <- if (isTRUE(all.equal(now_df,prev_df,check.attributes=FALSE))) 'same_support' else 'DIFFERENT_SUPPORT'
    if (comp_08d=='DIFFERENT_SUPPORT')
      warning('Selecao dos doadores difere de 08d; nao fazer comparacao pareada sem auditar.')
  }
}

pred_list <- list(); run_list <- list(); marker_list <- list(); fold_list <- list()
status_list <- list(); truth_list <- list(); errors <- list()

# ------------------------------------------------------------------------------
# 2. Uma pessoa excluida por fold; reconstruir a referencia de zero
# ------------------------------------------------------------------------------
for (ii in seq_len(nrow(chosen))) {
  donor <- chosen$donor_id[ii]
  fold <- chosen$fold[ii]
  cat('\n------------------------------------------------------\n')
  cat('08f | ',fold,' (treino: 23 doadores)\n',sep='')
  tryCatch({
    train <- which(meta$donor_id!=donor)
    test <- which(meta$donor_id==donor)
    if (length(intersect(train,test)) || data.table::uniqueN(meta$donor_id[train])!=23L ||
        !setequal(unique(meta$cell_type[train]),TYPES) ||
        !setequal(unique(meta$cell_type[test]),TYPES))
      stop('Separacao dos doadores / oito tipos invalida')

    train_type <- match(meta$cell_type[train],TYPES)
    test_type <- match(meta$cell_type[test],TYPES)
    n_train_type <- tabulate(train_type,nbins=length(TYPES))
    n_test_type <- tabulate(test_type,nbins=length(TYPES))
    if (any(n_train_type<=0) || any(n_test_type<=0)) stop('Tipo ausente em treino/teste')

    # REF 1: somar counts por tipo -> CPM pelo total de cada tipo.
    # Todas as operacoes de agregacao partem EXCLUSIVAMENTE do treino.
    agg_train <- Matrix::sparseMatrix(
      i=seq_along(train),j=train_type,x=rep(1,length(train)),
      dims=c(length(train),length(TYPES)))
    pool_counts <- sc[,train,drop=FALSE] %*% agg_train
    pool_lib <- as.numeric(Matrix::colSums(pool_counts))
    if (any(pool_lib<=0)) stop('Biblioteca agregada de treino nula')
    ref_pool_cpm <- sweep(as.matrix(pool_counts),2,pool_lib,'/')*1e6
    dimnames(ref_pool_cpm) <- list(rownames(sc),TYPES)

    # REF 2: cada nucleo recebe CPM da biblioteca INTEIRA; depois media.
    # IMPORTANTE: pesos excluem completamente os nucleos do fold de teste.
    w <- 1/(lib_nuc[train]*n_train_type[train_type])
    if (any(!is.finite(w))) stop('Pesos de CPM por nucleo invalidos')
    avg_train <- Matrix::sparseMatrix(
      i=seq_along(train),j=train_type,x=w,
      dims=c(length(train),length(TYPES)))
    ref_mean_cpm <- as.matrix(sc[,train,drop=FALSE] %*% avg_train)*1e6
    dimnames(ref_mean_cpm) <- list(rownames(sc),TYPES)

    # Misturas de um doador EXCLUIDO.
    # Media das contagens de cada gene/nucleo por tipo, sem incorporar
    # outros doadores e sem usar o doador de teste no ranking de marcadores.
    test_means <- sapply(seq_along(TYPES),function(k) {
      jj <- test[test_type==k]
      as.numeric(Matrix::rowMeans(sc[,jj,drop=FALSE]))
    })
    dimnames(test_means) <- list(rownames(sc),TYPES)
    test_counts <- round(test_means %*% t(VIRTUAL_COUNTS))
    sample_ids <- paste(fold,rownames(COMPOSITIONS),sep='_')
    colnames(test_counts) <- sample_ids
    rownames(test_counts) <- rownames(sc)
    test_lib <- as.numeric(colSums(test_counts))
    if (any(test_lib<=0) || any(!is.finite(test_lib))) stop('Mistura sintetica vazia')

    truth <- make_truth(test_means,sample_ids)
    if (!identical(rownames(truth$nuc),sample_ids) ||
        !identical(rownames(truth$rna),sample_ids)) stop('Verdade sintética desalinhada')

    # Universo de genes: somente genes de referencia com expressao de treino
    # e de pelo menos uma das misturas; equivalente ao filtro 'allgenes' 08e.
    # A filtragem das misturas nao utiliza labels clinicos (nao existem).
    keep <- rowSums(ref_pool_cpm)>0 & rowSums(test_counts)>0
    genes <- rownames(sc)[keep]
    if (length(genes)<1000L || anyDuplicated(genes)) stop('Universo de genes insuficiente')
    mix_log <- t(log2(sweep(test_counts[genes,,drop=FALSE],2,test_lib,'/')*1e6+1))
    if (!identical(rownames(mix_log),sample_ids) || !identical(colnames(mix_log),genes))
      stop('Genes ou amostras do bulk sintetico desalinhados')
    pure <- setNames(lapply(seq_along(TYPES),function(k) as.integer(k)),TYPES)

    t_truth <- data.table::melt(as.data.table(truth$nuc,keep.rownames='sample_id'),
      id.vars='sample_id',variable.name='cell_type',value.name='truth_nucleus')
    t_rna <- data.table::melt(as.data.table(truth$rna,keep.rownames='sample_id'),
      id.vars='sample_id',variable.name='cell_type',value.name='truth_rna')
    t_truth <- merge(t_truth,t_rna,by=c('sample_id','cell_type'),sort=FALSE)
    t_truth[,`:=`(fold=fold,scenario=sub(paste0('^',fold,'_'),'',sample_id))]
    if (nrow(t_truth)!=nrow(COMPOSITIONS)*length(TYPES))
      stop('Verdade synthetica nao tem 72 linhas por fold')
    truth_list[[length(truth_list)+1L]] <- t_truth

    ref_by_mode <- list(pooled_counts=ref_pool_cpm,mean_nucleus_cpm=ref_mean_cpm)
    for (mode in MODES) {
      cat('  ',mode,' | genes=',length(genes),'\n',sep='')
      ref_log <- t(log2(ref_by_mode[[mode]][genes,,drop=FALSE]+1))
      if (!identical(rownames(ref_log),TYPES) ||
          !identical(colnames(ref_log),colnames(mix_log)) ||
          any(!is.finite(ref_log))) stop('Matriz referencia desalinhada')

      # AUDITORIA SEM VAZAMENTO: ranquear marcadores somente com as 8 linhas
      # da referencia criada com 23 doadores. Depois conferir que os 100
      # marcadores escolhidos automaticamente pelo dtangle sao identicos.
      ranked <- dtangle::find_markers(
        Y=ref_log,pure_samples=pure,
        data_type=DATA_TYPE,marker_method=MARKER_METHOD)
      if (is.null(ranked$L) || length(ranked$L)!=length(TYPES))
        stop('Ranking find_markers nao retornou oito listas')
      train_markers <- lapply(ranked$L,function(m) head(as.integer(m),N_MARKERS))
      if (any(lengths(train_markers)!=N_MARKERS))
        stop('Referencia treino nao selecionou 100 marcadores por tipo')
      t0 <- proc.time()[['elapsed']]
      fit <- dtangle::dtangle(
        Y=mix_log,references=ref_log,pure_samples=pure,
        markers=NULL,n_markers=N_MARKERS,
        marker_method=MARKER_METHOD,data_type=DATA_TYPE)
      seconds <- as.numeric(proc.time()[['elapsed']]-t0)
      est <- verify_prop(fit$estimates,sample_ids,paste(fold,mode))
      if (length(fit$markers)!=length(TYPES)) stop('fit sem 8 listas de marcadores')
      for (k in seq_along(TYPES)) {
        ind <- as.integer(fit$markers[[k]])
        if (length(ind)!=N_MARKERS || anyNA(ind) || any(ind<1L | ind>length(genes)))
          stop('Indice de marcador invalido em ',TYPES[k])
        if (!identical(ind,train_markers[[k]]))
          stop('ATENCAO: marcadores automaticos divergem da selecao apenas treino em ',
               fold,'/',mode,'/',TYPES[k],'. Nao interpretar ate resolver o leakage.')
        marker_list[[length(marker_list)+1L]] <- data.table(
          fold=fold,reference_mode=mode,cell_type=TYPES[k],
          marker_rank=seq_len(N_MARKERS),gene_id=genes[ind],
          gene_column_index=ind)
      }

      est_long <- data.table::melt(as.data.table(est,keep.rownames='sample_id'),
        id.vars='sample_id',variable.name='cell_type',value.name='estimated')
      est_long <- merge(est_long,t_truth,
        by=c('sample_id','cell_type'),all.x=TRUE,sort=FALSE)
      if (nrow(est_long)!=nrow(est)*length(TYPES) ||
          anyNA(est_long$truth_nucleus) || anyNA(est_long$truth_rna))
        stop('Estimativas e verdades nao casaram exatamente')
      est_long[,reference_mode:=mode]
      # Evitar dependencia de colunas novas no MESMO data.table :=.
      est_long[,`:=`(
        error_nucleus_pp=100*(estimated-truth_nucleus),
        error_rna_pp=100*(estimated-truth_rna))]
      est_long[,`:=`(
        abs_error_nucleus_pp=abs(error_nucleus_pp),
        abs_error_rna_pp=abs(error_rna_pp))]
      pred_list[[length(pred_list)+1L]] <- est_long
      run_list[[length(run_list)+1L]] <- data.table(
        fold=fold,reference_mode=mode,n_train_donors=23L,
        train_nuclei=as.integer(length(train)),test_nuclei=as.integer(length(test)),
        genes_used=as.integer(length(genes)),n_markers_per_type=N_MARKERS,
        n_mixtures=as.integer(nrow(est)),gamma=as.numeric(fit$gamma),
        max_sum_deviation=as.numeric(max(abs(rowSums(est)-1))),
        mean_vascular_percent=as.numeric(100*mean(est[,'End']+est[,'Per'])),
        mae_nucleus_pp=as.numeric(mean(est_long$abs_error_nucleus_pp)),
        mae_rna_pp=as.numeric(mean(est_long$abs_error_rna_pp)),
        elapsed_seconds=seconds,marker_leakage_check='PASS',status='PASS')
      cat(sprintf('    MAE nucleo %.2f pp | MAE contagens %.2f pp | End+Per %.2f%%\n',
                  mean(est_long$abs_error_nucleus_pp),
                  mean(est_long$abs_error_rna_pp),
                  100*mean(est[,'End']+est[,'Per'])))
    }
    fold_list[[length(fold_list)+1L]] <- data.table(
      fold=fold,train_donors=23L,train_nuclei=as.integer(length(train)),
      test_nuclei=as.integer(length(test)),genes_used=as.integer(length(genes)),
      min_train_nuclei_by_type=as.integer(min(n_train_type)),
      min_heldout_nuclei_by_type=as.integer(min(n_test_type)))
    status_list[[length(status_list)+1L]] <- data.table(
      fold=fold,check='full_fold',status='PASS',info=paste(length(MODES),'modos'))

    # Checkpoint parcial dos resultados calculados (sem sobrescrever 08e).
    write_table(data.table::rbindlist(run_list,fill=TRUE),'08f_run_audit_partial.tsv')
    write_table(data.table::rbindlist(pred_list,fill=TRUE),'08f_predictions_partial.tsv')
  },error=function(e) {
    msg <- conditionMessage(e)
    message('FALHA em ',fold,': ',msg)
    status_list[[length(status_list)+1L]] <<- data.table(
      fold=fold,check='full_fold',status='FAIL',info=msg)
    errors[[length(errors)+1L]] <<- data.table(fold=fold,error=msg)
  })
}

# ------------------------------------------------------------------------------
# 3. Auditoria final e calculo de desempenho
# ------------------------------------------------------------------------------
status <- data.table::rbindlist(status_list,fill=TRUE)
write_table(status,'08f_status.tsv')
if (length(errors)) write_table(data.table::rbindlist(errors),'08f_ERRORS.tsv')
if (length(errors) || length(run_list)!=N_HOLDOUT*length(MODES))
  stop('08f INCOMPLETO. Confira 08f_status.tsv e 08f_ERRORS.tsv. NAO interprete resumos parciais.')

pred <- data.table::rbindlist(pred_list,fill=TRUE)
runs <- data.table::rbindlist(run_list,fill=TRUE)
marks <- data.table::rbindlist(marker_list,fill=TRUE)
truths <- data.table::rbindlist(truth_list,fill=TRUE)
fold_info <- data.table::rbindlist(fold_list,fill=TRUE)
expected_rows <- as.integer(N_HOLDOUT*length(MODES)*nrow(COMPOSITIONS)*length(TYPES))
if (nrow(pred)!=expected_rows || nrow(marks)!=N_HOLDOUT*length(MODES)*8L*N_MARKERS ||
    anyDuplicated(pred[,.(fold,reference_mode,sample_id,cell_type)]) ||
    any(pred$estimated<0 | pred$estimated>1) ||
    any(!is.finite(pred$estimated))) stop('Falha na auditoria final de linhas/estimativas')

metrics_global <- pred[,.(
  n_holdouts=as.integer(data.table::uniqueN(fold)),
  n_scenarios=as.integer(data.table::uniqueN(scenario)),
  n_predictions=as.integer(.N),
  mae_nucleus_pp=as.numeric(mean(abs_error_nucleus_pp)),
  rmse_nucleus_pp=as.numeric(sqrt(mean(error_nucleus_pp^2))),
  bias_nucleus_pp=as.numeric(mean(error_nucleus_pp)),
  mae_rna_pp=as.numeric(mean(abs_error_rna_pp)),
  rmse_rna_pp=as.numeric(sqrt(mean(error_rna_pp^2))),
  fraction_exact_zero=as.numeric(mean(estimated==0))
),by=reference_mode]

metrics_type <- pred[,.(
  n_predictions=as.integer(.N),
  mean_true_nucleus=as.numeric(mean(truth_nucleus)),
  mean_true_rna=as.numeric(mean(truth_rna)),
  mean_estimated=as.numeric(mean(estimated)),
  mae_nucleus_pp=as.numeric(mean(abs_error_nucleus_pp)),
  rmse_nucleus_pp=as.numeric(sqrt(mean(error_nucleus_pp^2))),
  bias_nucleus_pp=as.numeric(mean(error_nucleus_pp)),
  mae_rna_pp=as.numeric(mean(abs_error_rna_pp)),
  bias_rna_pp=as.numeric(mean(error_rna_pp)),
  fraction_zero_positive_truth=as.numeric(mean(estimated==0 & truth_nucleus>0)),
  spearman_nucleus=as.numeric(suppressWarnings(stats::cor(estimated,truth_nucleus,
                       method='spearman',use='complete.obs')))
),by=.(reference_mode,cell_type)]

metrics_holdout <- pred[,.(
  mae_nucleus_pp=as.numeric(mean(abs_error_nucleus_pp)),
  mae_rna_pp=as.numeric(mean(abs_error_rna_pp)),
  mae_In_pp=as.numeric(mean(abs_error_nucleus_pp[cell_type=='In'])),
  mae_Opc_pp=as.numeric(mean(abs_error_nucleus_pp[cell_type=='Opc'])),
  mae_End_Per_pp=as.numeric(mean(abs_error_nucleus_pp[cell_type %in% c('End','Per')])),
  mean_End_Per_percent=as.numeric(100*(mean(estimated[cell_type=='End'])+mean(estimated[cell_type=='Per'])))
),by=.(fold,reference_mode)]

metrics_scenario <- pred[,.(
  mae_nucleus_pp=as.numeric(mean(abs_error_nucleus_pp)),
  mae_rna_pp=as.numeric(mean(abs_error_rna_pp)),
  mae_In_pp=as.numeric(mean(abs_error_nucleus_pp[cell_type=='In'])),
  mae_Opc_pp=as.numeric(mean(abs_error_nucleus_pp[cell_type=='Opc'])),
  mae_End_Per_pp=as.numeric(mean(abs_error_nucleus_pp[cell_type %in% c('End','Per')]))
),by=.(reference_mode,scenario)]

# Comparacao PAREADA por tipo/mistura entre os dois modos. Nenhum p-valor:
# unidade biologica independente = doador (apenas quatro doadores por padrao).
paired <- NULL
if (setequal(MODES,REF_MODES)) {
  a <- pred[reference_mode=='pooled_counts',
            .(fold,scenario,sample_id,cell_type,
              error_pooled=abs_error_nucleus_pp)]
  b <- pred[reference_mode=='mean_nucleus_cpm',
            .(fold,scenario,sample_id,cell_type,
              error_mean=abs_error_nucleus_pp)]
  paired <- merge(a,b,by=c('fold','scenario','sample_id','cell_type'),all=TRUE)
  if (nrow(paired)!=N_HOLDOUT*nrow(COMPOSITIONS)*length(TYPES) ||
      anyNA(paired$error_pooled) || anyNA(paired$error_mean))
    stop('Comparacao pareada incompleta')
  paired[, delta_mae_pp:=error_mean-error_pooled]
  write_table(paired,'08f_paired_errors.tsv')
  write_table(paired[,.(
    n_pairs=as.integer(.N),
    mean_delta_pp=as.numeric(mean(delta_mae_pp)),
    fraction_pooled_better=as.numeric(mean(delta_mae_pp>0))
  ),by=cell_type],'08f_paired_summary_by_type.tsv')
}

write_table(pred,'08f_predictions_long.tsv')
write_table(runs,'08f_run_audit.tsv')
write_table(marks,'08f_marker_assignments.tsv')
write_table(truths,'08f_truths.tsv')
write_table(fold_info,'08f_fold_audit.tsv')
write_table(metrics_global,'08f_metrics_global.tsv')
write_table(metrics_type,'08f_metrics_by_type.tsv')
write_table(metrics_holdout,'08f_metrics_by_holdout.tsv')
write_table(metrics_scenario,'08f_metrics_by_scenario.tsv')
write_table(data.table(
  dataset='Mathys2019',reference='pathology_free_main24',
  dtangle_version=as.character(utils::packageVersion('dtangle')),
  fold_count=N_HOLDOUT,scenario_count=nrow(COMPOSITIONS),
  n_markers=N_MARKERS,marker_method=MARKER_METHOD,data_type=DATA_TYPE,
  modes=paste(MODES,collapse=','),
  reference_normalization='log2(CPM+1)',
  quantile_normalization=FALSE,
  holdout_support_match_08d=comp_08d,
  note='Validacao sintética INTRA-plataforma; nao teste com bulk real conhecido'
),'08f_parameters.tsv')

# ------------------------------------------------------------------------------
# 4. Figuras: a falha grafica nao apaga TSVs nem invalida os ajustes
# ------------------------------------------------------------------------------
plot_status <- 'PASS'
tryCatch({
  g1 <- ggplot(metrics_global,aes(x=reference_mode,y=mae_nucleus_pp,fill=reference_mode))+
    geom_col(width=.62)+geom_text(aes(label=sprintf('%.2f pp',mae_nucleus_pp)),
                 vjust=-.55,size=4)+
    scale_y_continuous(expand=expansion(mult=c(0,.15)))+
    theme_minimal(base_size=13)+theme(legend.position='none',
      plot.margin=margin(14,20,18,15))+
    labs(title='08f | dtangle: erro absoluto medio da validacao LODO',
         subtitle='Oito perfis, todos os genes elegiveis e 100 marcadores automaticos/tipo',
         x='Representacao da referencia',y='MAE por fracao de nucleos (pp)')
  ggsave(file.path(PLOTS,'08f_mae_global.png'),g1,width=11,height=6.5,dpi=180)

  g2 <- ggplot(metrics_type,aes(x=cell_type,y=mae_nucleus_pp,fill=reference_mode))+
    geom_col(position=position_dodge(width=.8),width=.7)+
    theme_minimal(base_size=12)+theme(plot.margin=margin(12,20,18,15))+
    labs(title='08f | Erro por tipo celular',subtitle='MAE em pontos percentuais',
         x='Tipo celular',y='MAE da fracao de nucleos (pp)',fill='Referencia')
  ggsave(file.path(PLOTS,'08f_mae_by_type.png'),g2,width=13,height=7,dpi=180)

  g3 <- ggplot(pred[cell_type %in% c('In','Opc','End','Per')],
               aes(x=100*truth_nucleus,y=100*estimated,color=reference_mode))+
    geom_point(alpha=.58,size=1.5)+geom_abline(slope=1,intercept=0,linetype=2)+
    facet_wrap(~cell_type,scales='free',ncol=2)+
    theme_minimal(base_size=12)+theme(plot.margin=margin(12,20,18,15))+
    labs(title='08f | Recuperacao dos tipos celulares criticos',
         subtitle='Linha pontilhada: estimativa igual a verdade; cores: referencia',
         x='Fracao nuclear conhecida (%)',y='Proporcao estimada (%)',color='Referencia')
  ggsave(file.path(PLOTS,'08f_truth_vs_estimated_critical.png'),g3,
         width=12,height=9,dpi=180)

  g4 <- ggplot(metrics_type,aes(x=cell_type,y=bias_nucleus_pp,fill=reference_mode))+
    geom_hline(yintercept=0,linetype=2,color='grey45')+
    geom_col(position=position_dodge(width=.8),width=.7)+
    theme_minimal(base_size=12)+theme(plot.margin=margin(12,20,18,15))+
    labs(title='08f | Vies sistematico por tipo celular',
         subtitle='Acima de zero = superestimacao; abaixo = subestimacao',
         x='Tipo celular',y='Vies medio (pp)',fill='Referencia')
  ggsave(file.path(PLOTS,'08f_bias_by_type.png'),g4,width=13,height=7,dpi=180)
},error=function(e) {
  plot_status <<- paste('FAIL:',conditionMessage(e))
  warning('Falha apenas nas figuras: ',conditionMessage(e))
})
write_table(data.table(check='plots',status=plot_status),'08f_plot_status.tsv')
capture.output(sessionInfo(),file=file.path(OUT,'08f_sessionInfo.txt'))
cat('\n=== 08f CONCLUIDO ===\n')
print(metrics_global[,.(reference_mode,n_holdouts,n_scenarios,n_predictions,
                        mae_nucleus_pp,rmse_nucleus_pp,mae_rna_pp)])
cat('\nResumo In/Opc/Oli/End/Per:\n')
print(metrics_type[cell_type %in% c('In','Opc','Oli','End','Per'),
  .(reference_mode,cell_type,mae_nucleus_pp,bias_nucleus_pp,mae_rna_pp)])
cat('\nSaidas: ',OUT,'\n',sep='')
cat('NOTA: simulacao Mathys->Mathys nao valida exatidao nos bulks reais.\n')
