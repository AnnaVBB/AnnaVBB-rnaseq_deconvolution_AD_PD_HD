#!/usr/bin/env Rscript
# ==============================================================================
# 07c_auditar_BisqueRNA.R — AUDITORIA DE CALIBRAÇÃO/COMPARABILIDADE
# ============================================================================== 
# NÃO altera a referência, a Fase 06, os outputs originais do Bisque (07), nem
# renv.lock. A comparação de cenários é exploratória, não seleção por p-valor.
#
# ANALISE A (sem recálculo):
#   1. Para 4 coortes, comparar proporção média estimada pelo Bisque à
#      proporção média por doador de 'sc.props' usada internamente pelo pacote.
#   2. Comparar também com a frequência de núcleos do snRNA-seq (pooling).
#   3. Quantificar se as médias são praticamente iguais entre coortes.
#   4. Checar amostras de controle compartilhadas HD/PD, quando presentes.
#
# ANALISE B (opcional, 1 coorte; RECALCULA APENAS UM AJUSTE):
#   Comparar o ajuste v4/v5 principal (old.cpm=FALSE, denominador completo)
#   com normalização apenas pelos genes de top300 (old.cpm=TRUE), que é a
#   opção default de várias versões do BisqueRNA. NÃO usar esta análise
#   para escolher a condição que parece mais plausível biologicamente.
#
# EXECUÇÃO:
#   Rscript scripts/phase02_deconvolution/07c_auditar_BisqueRNA.R
#   BISQUE_RUN_LEGACY=1 Rscript scripts/phase02_deconvolution/07c_auditar_BisqueRNA.R
#   BISQUE_RUN_LEGACY=1 BISQUE_LEGACY_DATASET=ROSMAP Rscript ...
#
# OUTPUTS: results/deconvolution/BisqueRNA_Audit/07c_*.tsv, .txt
# ==============================================================================

ROOT <- file.path('results','deconvolution','BisqueRNA')
OUT <- file.path('results','deconvolution','BisqueRNA_Audit')
DATASETS <- c('ROSMAP','GSE125583','GSE64810','GSE68719')
CT <- c('Ex','In','Ast','Oli','Opc','Mic','End','Per')
RUN_LEGACY <- identical(trimws(Sys.getenv('BISQUE_RUN_LEGACY','0')),'1')
LEGACY_DS <- trimws(Sys.getenv('BISQUE_LEGACY_DATASET','GSE64810'))
if (!LEGACY_DS %in% DATASETS) stop('BISQUE_LEGACY_DATASET inválido.')
req <- c('data.table','Biobase','Matrix','BisqueRNA','limSolve')
missing <- req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop('Pacotes ausentes: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(data.table);library(Biobase);library(Matrix)})
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(x,name) data.table::fwrite(x,file.path(OUT,name),sep='\t',na='NA')
report <- function(...) message(sprintf(...))

# ------------------------------------------------------------------------------
# ANALISE A — Verificação de compatibilidade e ancoragem nas médias da referência
# ------------------------------------------------------------------------------
prior_tables <- list(); sample_tables <- list(); run_tables <- list()
for (ds in DATASETS) {
  ddir <- file.path(ROOT,ds,'main24','top300')
  pth <- file.path(ddir,'bisque_result.rds')
  wide_file <- file.path(ddir,'proportions_wide.tsv')
  if (!all(file.exists(c(pth,wide_file)))) {
    stop('Arquivos do 07 ausentes para ',ds,': ',paste(c(pth,wide_file),collapse=', '))
  }
  fit <- readRDS(pth)
  if (is.null(fit$sc.props) || is.null(fit$bulk.props)) {
    stop(ds,': resultado sem bulk.props ou sc.props.')
  }
  sc <- as.matrix(fit$sc.props)
  bp <- as.matrix(fit$bulk.props)
  if (!setequal(rownames(sc),CT) || !setequal(rownames(bp),CT))
    stop(ds,': tipos celulares incompatíveis com o plano.')
  sc <- sc[CT,,drop=FALSE]
  bp <- bp[CT,,drop=FALSE]
  if (any(!is.finite(sc)) || any(!is.finite(bp))) stop(ds,': valores não finitos.')
  if (max(abs(colSums(bp)-1))>1e-6) stop(ds,': proporções não somam 1.')
  wide <- fread(wide_file)
  if (!all(c('sample_id','group',CT) %in% names(wide)))
    stop(ds,': proportions_wide.tsv incompleto.')
  if (!setequal(wide$sample_id,colnames(bp))) stop(ds,': IDs divergentes entre fit e TSV.')
  for (ct in CT) {
    # A ordem dos samples pode variar; corrigir por chave, não posição.
    if (max(abs(wide[[ct]][match(colnames(bp),wide$sample_id)]-bp[ct,]))>1e-9)
      stop(ds,': TSV e fit diferem para ',ct)
  }
  means <- rowMeans(bp)
  prior <- rowMeans(sc)
  prior_tables[[ds]] <- data.table(dataset=ds,cell_type=CT,
                                  mean_bulk=as.numeric(means),
                                  mean_reference_donor=as.numeric(prior),
                                  difference_pp=100*as.numeric(means-prior),
                                  absolute_difference_pp=100*abs(as.numeric(means-prior)),
                                  n_bulk=ncol(bp),n_ref_donors=ncol(sc))
  sample_dt <- as.data.table(t(bp))
  setnames(sample_dt,CT)
  sample_dt[,`:=`(dataset=ds,sample_id=colnames(bp),
                  group=wide$group[match(colnames(bp),wide$sample_id)])]
  setcolorder(sample_dt,c('dataset','sample_id','group',CT))
  sample_tables[[ds]] <- sample_dt
  run_tables[[ds]] <- data.table(dataset=ds,n_samples=ncol(bp),
                                 n_reference_donors=ncol(sc),
                                 bisque_version=as.character(packageVersion('BisqueRNA')),
                                 min_value=min(bp),max_value=max(bp),
                                 max_sum_error=max(abs(colSums(bp)-1)),
                                 max_abs_mean_minus_ref_pp=max(abs(100*(means-prior))))
  rm(fit,sc,bp); invisible(gc(FALSE))
}
anchor <- rbindlist(prior_tables)
ref_meta_file <- file.path('data','reference','mathys2019','phase02_prepared','reference',
                          'mathys2019_pathology_free_cell_metadata_aligned.tsv')
if (file.exists(ref_meta_file)) {
  m <- fread(ref_meta_file)
  if ('cell_type' %in% names(m)) {
    counts <- m[,.(n_nuclei=.N),by=cell_type]
    counts[,pooled_fraction:=n_nuclei/sum(n_nuclei)]
    anchor <- merge(anchor,counts,by='cell_type',all.x=TRUE,sort=FALSE)
    anchor[,difference_to_pooled_pp:=100*(mean_bulk-pooled_fraction)]
  }
}
setorder(anchor,dataset,cell_type)
write_tsv(anchor,'07c_reference_anchor_audit.tsv')
write_tsv(rbindlist(run_tables),'07c_run_integrity.tsv')

stability <- anchor[,.(min_cohort_mean=min(mean_bulk),
                       max_cohort_mean=max(mean_bulk),
                       range_cohort_mean_pp=100*(max(mean_bulk)-min(mean_bulk))),
                    by=cell_type]
setorder(stability,-range_cohort_mean_pp)
write_tsv(stability,'07c_cross_cohort_mean_stability.tsv')

all_samples <- rbindlist(sample_tables,fill=TRUE)
# Observação: as mesmas amostras controle aparecem em HD/PD; estimativas
# podem mudar porque use.overlap=FALSE padroniza sobre TODA A COORTE.
a <- all_samples[dataset=='GSE64810' & group=='Control']
b <- all_samples[dataset=='GSE68719' & group=='Control']
shared <- merge(a,b,by='sample_id',suffixes=c('_HD','_PD'))
if (nrow(shared)) {
  per_type <- rbindlist(lapply(CT,function(ct) {
    dif <- abs(shared[[paste0(ct,'_HD')]]-shared[[paste0(ct,'_PD')]])
    data.table(cell_type=ct,n_shared_controls=nrow(shared),
               mean_absolute_diff_pp=100*mean(dif),max_absolute_diff_pp=100*max(dif))
  }))
  write_tsv(per_type,'07c_shared_HD_PD_controls.tsv')
}

report('=== AUDITORIA 07c A (sem novo ajuste) ===')
print(anchor[,.(dataset,cell_type,mean_bulk,mean_reference_donor,difference_pp)])
report('Amplitude entre médias das 4 coortes, por tipo (pp):')
print(stability)
report('Atenção: proximidade à média da referência é compatível com a transformação')
report('semissupervisionada do Bisque; não valida exatidão da composição.')

# ------------------------------------------------------------------------------
# ANALISE B — Diferença da implementação old.cpm, apenas se solicitada.
# Para economizar tempo/memória, executar uma coorte de cada vez.
# ------------------------------------------------------------------------------
if (RUN_LEGACY) {
  report('=== AUDITORIA 07c B: old.cpm=TRUE na coorte %s ===',LEGACY_DS)
  ddir <- file.path(ROOT,LEGACY_DS,'main24','top300')
  ref_path <- file.path('data','reference','mathys2019','phase02_prepared','reference',
                        'mathys2019_pathology_free_geneid_counts_sparse.rds')
  meta_path <- file.path('data','reference','mathys2019','phase02_prepared','reference',
                         'mathys2019_pathology_free_cell_metadata_aligned.tsv')
  bulk_path <- file.path('data','processed','deconvolution',LEGACY_DS,'bulk_counts_top300.rds')
  bmeta_path <- file.path('data','processed','deconvolution',LEGACY_DS,'metadata_deconvolution.tsv')
  needed <- c(ref_path,meta_path,bulk_path,bmeta_path)
  if (any(!file.exists(needed))) stop('Arquivos ausentes: ',paste(needed[!file.exists(needed)],collapse='; '))
  ref <- readRDS(ref_path)
  if (!inherits(ref,'dgCMatrix')) ref <- as(ref,'dgCMatrix')
  meta <- fread(meta_path)
  ix <- match(colnames(ref),meta$TAG)
  if (anyNA(ix)) stop('Metadados Mathys não alinham núcleos.')
  meta <- meta[ix]
  if (anyNA(meta$cell_type) || anyNA(meta$donor_id)) stop('Metadados Mathys incompletos.')
  bulk <- as.matrix(readRDS(bulk_path))
  bmeta <- fread(bmeta_path)
  idx <- match(colnames(bulk),bmeta$sample_id)
  if (anyNA(idx)) stop('Metadados bulk não alinham colunas.')
  bmeta <- bmeta[idx]
  genes <- rownames(bulk)[rownames(bulk) %in% rownames(ref)]
  pos <- rowSums(bulk[genes,,drop=FALSE])>0 & Matrix::rowSums(ref[genes,,drop=FALSE])>0
  genes <- genes[pos]
  if (length(genes)<500) stop('Menos de 500 genes em comum.')
  report('Genes: %d | Núcleos: %d | Bulk: %d',length(genes),ncol(ref),ncol(bulk))
  sc_mat <- as.matrix(ref[genes,,drop=FALSE])
  sc_pd <- data.frame(cellType=factor(meta$cell_type,levels=CT),
                      SubjectName=factor(meta$donor_id),row.names=colnames(ref),check.names=FALSE)
  bl_pd <- data.frame(group=bmeta$group,row.names=colnames(bulk),check.names=FALSE)
  sc_eset <- Biobase::ExpressionSet(assayData=sc_mat,
                       phenoData=Biobase::AnnotatedDataFrame(sc_pd))
  bl_eset <- Biobase::ExpressionSet(assayData=bulk[genes,,drop=FALSE],
                       phenoData=Biobase::AnnotatedDataFrame(bl_pd))
  rm(sc_mat,ref,bulk); invisible(gc(FALSE))
  t0 <- proc.time()[['elapsed']]
  legacy <- BisqueRNA::ReferenceBasedDecomposition(
    bulk.eset=bl_eset,sc.eset=sc_eset,markers=genes,
    cell.types='cellType',subject.names='SubjectName',
    use.overlap=FALSE,verbose=TRUE,old.cpm=TRUE)
  duration <- proc.time()[['elapsed']]-t0
  old <- legacy$bulk.props[CT,,drop=FALSE]
  original <- readRDS(file.path(ddir,'bisque_result.rds'))$bulk.props[CT,,drop=FALSE]
  if (!setequal(colnames(old),colnames(original))) stop('Conjunto de IDs distinto entre ajustes.')
  old <- old[,colnames(original),drop=FALSE]
  if (any(!is.finite(old)) || max(abs(colSums(old)-1))>1e-6)
    stop('Ajuste legacy numericamente inválido.')
  diff <- old-original
  comparison <- data.table(dataset=LEGACY_DS,cell_type=CT,
    mean_full_denominator=as.numeric(rowMeans(original)),
    mean_marker_denominator=as.numeric(rowMeans(old)),
    mean_difference_pp=as.numeric(100*rowMeans(diff)),
    mean_absolute_difference_pp=as.numeric(100*rowMeans(abs(diff))),
    max_absolute_difference_pp=as.numeric(100*apply(abs(diff),1,max)),
    fraction_zero_full=as.numeric(rowMeans(original<=1e-10)),
    fraction_zero_markers=as.numeric(rowMeans(old<=1e-10)),
    spearman=vapply(CT,function(ct) {
      if (sd(old[ct,])<1e-12 || sd(original[ct,])<1e-12) return(NA_real_)
      suppressWarnings(cor(old[ct,],original[ct,],method='spearman'))
    },numeric(1)))
  write_tsv(comparison,paste0('07c_oldCPM_sensitivity_',LEGACY_DS,'.tsv'))
  wide <- as.data.table(t(old))
  setnames(wide,CT)
  wide[,sample_id:=colnames(old)]
  setcolorder(wide,c('sample_id',CT))
  write_tsv(wide,paste0('07c_legacy_proportions_',LEGACY_DS,'.tsv'))
  write_tsv(data.table(dataset=LEGACY_DS,genes_used=length(legacy$genes.used),
                       n_bulk=ncol(old),elapsed_seconds=duration,
                       method='use.overlap=FALSE; old.cpm=TRUE; no library remainder'),
            paste0('07c_legacy_run_audit_',LEGACY_DS,'.tsv'))
  report('Diferenças por tipo (pontos percentuais):')
  print(comparison[,.(cell_type,mean_full_denominator,mean_marker_denominator,
                       mean_absolute_difference_pp,spearman)])
  report('Não substituir o 07 principal pela condição legacy sem decisão metodológica.')
}

report('Auditoria 07c concluída. Saídas: %s',OUT)
