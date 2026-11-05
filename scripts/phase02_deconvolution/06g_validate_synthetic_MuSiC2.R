#!/usr/bin/env Rscript
# ==============================================================================
# 06g_validate_synthetic_MuSiC2.R | REVISAO 03: counts()/assay() namespace correto
# Validação interna MuSiC2: composição conhecida + DE celular simulada
# Complementa 06c–06f. NÃO altera nenhuma saída dos scripts 06–06f, 09 ou 12.
#
# EXECUÇÃO (a partir da raiz do projeto):
#   Rscript scripts/phase02_deconvolution/06g_validate_synthetic_MuSiC2.R --smoke
#   Rscript scripts/phase02_deconvolution/06g_validate_synthetic_MuSiC2.R --full
#   Rscript scripts/phase02_deconvolution/06g_validate_synthetic_MuSiC2.R --full --top300
#
# --smoke : teste estrutural mínimo, não interpretável biologicamente.
# --full  : quatro doadores, cinco perfis, quatro cenários (ver abaixo).
# --top300: sensibilidade adicional. Marcadores reconstruídos NO TREINO.
# O tempo de execução pode ser elevado; cada ajuste é salvo isoladamente.
#
# CENÁRIOS: baseline, composition_only, expression_only, both.
# Controle/caso são rótulos SINTÉTICOS, não diagnósticos reais.
# As frações nucleares são conhecidas, mas a contribuição de contagens RNA
# não é igual a elas: registramos AMBAS as verdades de referência.
#
# O procedimento do MuSiC2 é condicionado à API xuranw/MuSiC::music2_prop.
# Não há qualquer promessa de acurácia para bulk real. Instabilidade de
# convergência é registrada, NÃO ocultada nem corrigida artificialmente.
# ==============================================================================

args <- commandArgs(trailingOnly=TRUE)
if (!any(args %in% c('--smoke','--full'))) stop('Use --smoke ou --full.')
if (all(c('--smoke','--full') %in% args)) stop('Escolha só um modo.')
MODE <- if ('--smoke' %in% args) 'smoke' else 'full'
WITH_TOP300 <- '--top300' %in% args
N_FOLDS <- if (MODE == 'smoke') 1L else 4L
N_PER_GROUP <- if (MODE == 'smoke') 4L else 20L
PROFILE_NAMES <- if (MODE == 'smoke') 'cortical_like' else c('equal','cortical_like','inhibitory_enriched','opc_enriched','vascular_enriched')
SCENARIOS <- if (MODE == 'smoke') 'baseline' else c('baseline','composition_only','expression_only','both')
GENE_SETS <- if (WITH_TOP300) c('allgenes','top300_train_only') else 'allgenes'
SEED <- 62027L
VIRTUAL_N <- 1000
DE_GENE_FRACTION <- 0.10 # 10% dos genes elegíveis no tipo escolhido (apenas caso sintético)
FOLD_CHANGE <- 2.0       # cenário DE sintético: expressão duplicada
COMP_SHIFT <- 0.05       # deslocar 5 pontos percentuais de Ex para Ast no grupo caso
DIR <- file.path('results','deconvolution','MuSiC2_Synthetic_Validation_06g', MODE)
dir.create(DIR,recursive=TRUE,showWarnings=FALSE)

# MuSiC usa o accessor counts(), cujo generico pertence a BiocGenerics.
# Reproduzir as dependencias do script 06 e realizar verificacoes com
# SummarizedExperiment::assay(x, 'counts') para evitar namespaces incorretos.
for (pkg in c('data.table','Matrix','BiocGenerics','SingleCellExperiment',
              'SummarizedExperiment','S4Vectors','MuSiC')) {
 if (!requireNamespace(pkg, quietly=TRUE)) stop('Pacote R ausente: ', pkg)
}
suppressPackageStartupMessages({
 library(data.table)
 library(Matrix)
 library(BiocGenerics)
 library(SingleCellExperiment)
 library(SummarizedExperiment)
 library(S4Vectors)
 library(MuSiC)
})
# 'counts()' e um accessor generico exportado por BiocGenerics, e NAO por
# SummarizedExperiment. O accessor com namespace explicitamente correto e:
# SummarizedExperiment::assay(sce, 'counts').
if (!'counts' %in% getNamespaceExports('BiocGenerics')) {
 stop('BiocGenerics::counts indisponivel; verifique a instalacao do Bioconductor.')
}
if (!'assay' %in% getNamespaceExports('SummarizedExperiment')) {
 stop('SummarizedExperiment::assay indisponivel; verifique a instalacao.')
}
message('Pacotes carregados: BiocGenerics::counts() e ',
        'SummarizedExperiment::assay() disponiveis.')
if (!exists('music2_prop', envir=asNamespace('MuSiC'),inherits=FALSE) ||
    !exists('music_prop', envir=asNamespace('MuSiC'),inherits=FALSE)) stop('Instale a versão de MuSiC usada nos scripts 06/06f.')
req2 <- c('bulk.control.mtx','bulk.case.mtx','sc.sce','clusters','samples','select.ct','expr_low','prop_r','eps_c','eps_r','n_resample','sample_prop','cutoff_expr','cutoff_c','cutoff_r','maxiter','markers','ct.cov','centered','normalize')
if (!all(req2 %in% names(formals(MuSiC::music2_prop)))) stop('API music2_prop incompatível com script 06. Verifique formals(MuSiC::music2_prop).')
CELL <- c('Ex','In','Ast','Oli','Opc','Mic','End','Per')
ROOT <- file.path('data','reference','mathys2019','phase02_prepared','reference')
REF_FILE <- file.path(ROOT,'mathys2019_pathology_free_geneid_counts_sparse.rds')
META_FILE <- file.path(ROOT,'mathys2019_pathology_free_cell_metadata_aligned.tsv')
for(f in c(REF_FILE,META_FILE)) if (!file.exists(f)) stop('Arquivo ausente: ',f)

compositions <- rbind(
 equal=c(Ex=.125,In=.125,Ast=.125,Oli=.125,Opc=.125,Mic=.125,End=.125,Per=.125),
 cortical_like=c(Ex=.38,In=.17,Ast=.12,Oli=.20,Opc=.05,Mic=.04,End=.02,Per=.02),
 inhibitory_enriched=c(Ex=.22,In=.30,Ast=.12,Oli=.20,Opc=.06,Mic=.05,End=.03,Per=.02),
 opc_enriched=c(Ex=.20,In=.12,Ast=.13,Oli=.22,Opc=.20,Mic=.07,End=.03,Per=.03),
 vascular_enriched=c(Ex=.23,In=.12,Ast=.15,Oli=.18,Opc=.06,Mic=.06,End=.10,Per=.10))
stopifnot(identical(colnames(compositions),CELL),all(abs(rowSums(compositions)-1)<1e-10))

ref <- readRDS(REF_FILE)
if (!inherits(ref,'dgCMatrix')) ref <- as(ref,'dgCMatrix')
meta <- fread(META_FILE)
if (!all(c('TAG','donor_id','cell_type') %in% names(meta))) stop('Metadados incompletos.')
ix <- match(colnames(ref),meta$TAG)
if (anyNA(ix) || anyDuplicated(meta$TAG)) stop('Núcleos sem correspondência 1:1.')
meta <- meta[ix]
stopifnot(identical(as.character(meta$TAG),colnames(ref)),nrow(ref)==17394L,ncol(ref)==35140L,uniqueN(meta$donor_id)==24L)
if(anyDuplicated(rownames(ref))) stop('GeneIDs duplicados')

support <- dcast(meta[, .N,by=.(donor_id,cell_type)],donor_id~cell_type,value.var='N',fill=0)
eligible <- support[rowSums(as.matrix(support[,..CELL])>0)==8]
if(nrow(eligible)<N_FOLDS)stop('Sem doadores suficientes com os oito tipos.')
set.seed(SEED)
holdouts <- as.character(sample(eligible$donor_id,N_FOLDS))
# Do not export identifying donor labels. Each fold is an independent unit.
fwrite(data.table(fold=seq_along(holdouts),n_test_nuclei=vapply(holdouts,function(d)sum(meta$donor_id==d),integer(1))),file.path(DIR,'06g_folds.tsv'),sep='\t')

# Helpers: explicit sample matching and numerical integrity.
standardize <- function(obj,ids,what) {
 if(is.null(obj))stop(what,': missing proportions')
 x <- as.matrix(obj)
 if(!is.null(rownames(x)) && all(ids %in% rownames(x))) x<-x[ids,,drop=FALSE]
 else if(!is.null(colnames(x)) && all(ids %in% colnames(x))) x<-t(x[,ids,drop=FALSE])
 else stop(what,': sample names not found')
 if(!setequal(colnames(x),CELL))stop(what,': cell types mismatch')
 x<-x[,CELL,drop=FALSE];storage.mode(x)<-'double'
 if(any(!is.finite(x))||any(x < -1e-8)||max(abs(rowSums(x)-1))>1e-5)stop(what,': invalid proportions')
 x
}

# Marker score is reconstructed on training donors only, reproducing core
# formula from script 04: specificity * log1p(mean_CPM) * detection_fraction.
training_top300 <- function(x,m) {
 donors <- unique(as.character(m$donor_id))
 means <- matrix(0,nrow(x),length(CELL),dimnames=list(rownames(x),CELL))
 det <- means
 for(k in seq_along(CELL)) {
  ct <- CELL[k]
  who <- unique(as.character(m$donor_id[m$cell_type==ct]))
  prof <- vapply(who, function(d) {
   j <- which(m$donor_id==d & m$cell_type==ct)
   cnt <- Matrix::rowSums(x[,j,drop=FALSE]); as.numeric(cnt)*1e6/max(1,sum(cnt))
  },numeric(nrow(x)))
  if(is.null(dim(prof)))prof<-matrix(prof,ncol=1)
  means[,k]<-rowMeans(prof);det[,k]<-rowMeans(prof>0)
 }
 sp <- means / pmax(rowSums(means),.Machine$double.eps)
 scored <- sp * log1p(means) * det
 genes <- character()
 for(k in seq_along(CELL)) {
  ord<-order(-scored[,k],-sp[,k],-means[,k],rownames(x))
  genes<-c(genes,rownames(x)[ord[seq_len(min(300,length(ord)))]] )
 }
 unique(genes)
}

# Build average count-per-nucleus profiles for an entirely excluded test donor.
means_for_donor <- function(d) {
 ans<-matrix(0,nrow(ref),8,dimnames=list(rownames(ref),CELL))
 for(k in seq_along(CELL)) {
  j<-which(meta$donor_id==d & meta$cell_type==CELL[k]);if(!length(j)) stop('Missing test cell type')
  ans[,k]<-as.numeric(Matrix::rowMeans(ref[,j,drop=FALSE]))
 }
 ans
}

# Differential expression genes drawn without inspecting downstream estimates.
# Vary affected cell types and direction; both up/down to avoid trivial bias.
make_DE <- function(means,fold_i,profile) {
 affected <- c('Ast','In','Opc','Mic')
 target <- if (profile=='inhibitory_enriched') 'In' else if(profile=='opc_enriched') 'Opc' else if(profile=='vascular_enriched') 'Mic' else 'Ast'
 k<-match(target,CELL)
 elig<-which(means[,k]>0)
 if(length(elig)<20)stop('Too few genes for DE simulation.')
 n<-max(20L,as.integer(length(elig)*DE_GENE_FRACTION))
 selected<-sample(elig,min(n,length(elig)))
 factor<-rep(1,nrow(means))
 factor[selected]<-rep(c(FOLD_CHANGE,1/FOLD_CHANGE),length.out=length(selected))
 list(type=target,genes=rownames(means)[selected],factor=factor)
}

# For each condition and sample: multinomial-style jitter of fractions around a
# fixed center, then independent Poisson counts at expected mean values.
# Here the RNA contribution is computed as the expected contribution to counts.
make_bulk <- function(means,composition,ids,modified=FALSE,DE=NULL) {
 n<-length(ids);out<-matrix(0,nrow(means),n,dimnames=list(rownames(means),ids))
 truth<-vector('list',n)
 for(j in seq_len(n)) {
  frac<-rgamma(length(CELL),shape=pmax(composition*200,0.1),rate=1)
  frac<-frac/sum(frac);names(frac)<-CELL
  m<-means
  if(modified && !is.null(DE)) m[,DE$type]<-m[,DE$type]*DE$factor
  contributions<-sweep(m,2,frac*VIRTUAL_N,'*')
  expected<-rowSums(contributions)
  out[,j]<-stats::rpois(length(expected),lambda=pmax(expected,0))
  rna<-colSums(contributions)/sum(contributions)
  truth[[j]]<-data.table(sample_id=ids[j],cell_type=CELL,true_nuclear=as.numeric(frac),true_rna=as.numeric(rna))
 }
 list(counts=out,truth=rbindlist(truth))
}

PROPS <- list();RUNS <- list();DESIGN <- list();DEGENES<-list();run_index<-0L
for(fi in seq_along(holdouts)) {
 donor<-holdouts[fi];label<-sprintf('F%02d',fi)
 train_idx<-which(meta$donor_id!=donor)
 trainmeta<-meta[train_idx]
 if(uniqueN(trainmeta$donor_id)!=23L||!setequal(unique(trainmeta$cell_type),CELL))stop('Train/test failure')
 test_means<-means_for_donor(donor)
 train_counts<-ref[,train_idx,drop=FALSE]
 for(gene_set in GENE_SETS) {
  genes<-if(gene_set=='allgenes')rownames(ref) else training_top300(train_counts,trainmeta)
  for(profile in PROFILE_NAMES) for(scenario in SCENARIOS) {
   run_index<-run_index+1L
   code<-paste(label,gene_set,profile,scenario,sep='__')
   output_run<-file.path(DIR,'runs',code);dir.create(output_run,recursive=TRUE,showWarnings=FALSE)
   if(file.exists(file.path(output_run,'DONE.tsv'))) {message('Skipping completed: ',code);next}
   # Fixed seeds per fold/scenario, reproducible for resumed runs.
   set.seed(SEED+fi*100000L+match(profile,rownames(compositions))*1000L+match(scenario,c('baseline','composition_only','expression_only','both'))*10L)
   base<-compositions[profile,]
   cas<-base
   has_shift<-scenario %in% c('composition_only','both')
   has_DE<-scenario %in% c('expression_only','both')
   if(has_shift) {cas['Ex']<-cas['Ex']-COMP_SHIFT;cas['Ast']<-cas['Ast']+COMP_SHIFT}
   if(any(cas<=0))stop('Negative or zero composition requested')
   de<-if(has_DE)make_DE(test_means,fi,profile) else NULL
   controls<-make_bulk(test_means,base,paste0(code,'__CTRL_',seq_len(N_PER_GROUP)))
   cases<-make_bulk(test_means,cas,paste0(code,'__CASE_',seq_len(N_PER_GROUP)),modified=has_DE,DE=de)
   design<-rbindlist(list(controls$truth[,group:='Control'],cases$truth[,group:='Case']))
   design[,`:=`(fold=fi,gene_set=gene_set,profile=profile,scenario=scenario)]
   fwrite(design,file.path(output_run,'06g_truth_long.tsv'),sep='\t')
   if(has_DE)fwrite(data.table(gene_id=de$genes,affected_cell_type=de$type,fold=fi,profile=profile,scenario=scenario),file.path(output_run,'06g_simulated_DE_genes.tsv'),sep='\t')
   # Subset train AND synthetic bulk to common, nonzero genes.
   cmat<-controls$counts[genes,,drop=FALSE]
   amat<-cases$counts[genes,,drop=FALSE]
   sc<-train_counts[genes,,drop=FALSE]
   keep <- Matrix::rowSums(sc)>0 & (rowSums(cmat)+rowSums(amat)>0)
   sc<-sc[keep,,drop=FALSE];cmat<-cmat[keep,,drop=FALSE];amat<-amat[keep,,drop=FALSE]
   if(nrow(sc)<100)stop('Less than 100 eligible genes in ',code)
   sce<-SingleCellExperiment::SingleCellExperiment(
     assays=list(counts=sc),colData=S4Vectors::DataFrame(
       cellType=factor(trainmeta$cell_type,levels=CELL),
       sampleID=factor(trainmeta$donor_id),row.names=trainmeta$TAG))
   colnames(sce)<-trainmeta$TAG;rownames(sce)<-rownames(sc)
   # Preflight: verificar que a SCE contem matriz de contagens e metadados.
   stopifnot('counts' %in% SummarizedExperiment::assayNames(sce))
   stopifnot(identical(dim(SummarizedExperiment::assay(sce, 'counts')),dim(sc)))
   stopifnot(identical(as.character(SummarizedExperiment::colData(sce)$cellType),
                       as.character(trainmeta$cell_type)))
   # Os dois algoritmos sao executados independentemente. Se um falhar,
   # registramos os resultados do outro, sem chamar isso de comparacao aprovada.
   t0<-Sys.time();fitM<-fit2<-NULL
   method_errors<-character()
   run_method <- function(label, fn) {
     tryCatch(fn(), error=function(e) {
       method_errors <<- c(method_errors, paste0(label, ': ', conditionMessage(e)))
       NULL
     })
   }
   fitM <- run_method('MuSiC', function() MuSiC::music_prop(
     bulk.mtx=cbind(cmat,amat),sc.sce=sce,clusters='cellType',samples='sampleID',
     select.ct=CELL,markers=NULL,cell_size=NULL,ct.cov=FALSE,
     iter.max=1000,nu=0.0001,eps=0.01,centered=FALSE,normalize=FALSE,
     verbose=FALSE))
   fit2 <- run_method('MuSiC2', function() MuSiC::music2_prop(
     bulk.control.mtx=cmat,bulk.case.mtx=amat,sc.sce=sce,clusters='cellType',
     samples='sampleID',select.ct=CELL,expr_low=20,prop_r=.10,eps_c=.05,
     eps_r=.01,n_resample=20,sample_prop=.50,cutoff_expr=.05,cutoff_c=.05,
     cutoff_r=.01,maxiter=200,markers=NULL,cell_size=NULL,ct.cov=FALSE,
     centered=FALSE,normalize=FALSE))
   # Extrair resultados inclusive quando um dos metodos falhar.
   ids<-c(colnames(cmat),colnames(amat))
   extract_fit <- function(obj, slot, name) {
     if (is.null(obj)) return(NULL)
     run_method(paste0('Extracao ',name), function() standardize(obj[[slot]],ids,name))
   }
   estimates <- list(
     MuSiC=extract_fit(fitM,'Est.prop.weighted','MuSiC'),
     MuSiC2=extract_fit(fit2,'Est.prop','MuSiC2'))
   completed_methods <- names(estimates)[!vapply(estimates,is.null,logical(1))]
   # NULL sem erro explicito nao pode ser registrado como ajuste bem sucedido.
   absent_methods <- setdiff(c('MuSiC','MuSiC2'),completed_methods)
   if (length(absent_methods) && !length(method_errors)) {
     method_errors <- paste0(absent_methods, ': retorno NULL sem estimativas.')
   }
   for(method in completed_methods) {
     run_method(paste0('Exportacao ',method),function() {
       zz<-as.data.table(estimates[[method]],keep.rownames='sample_id')
       zz<-melt(zz,id.vars='sample_id',variable.name='cell_type',value.name='estimate')
       zz[,cell_type:=as.character(cell_type)]
       z<-merge(zz,design[,.(sample_id,cell_type,group,true_nuclear,true_rna)],
                by=c('sample_id','cell_type'))
       if (nrow(z)!=length(ids)*length(CELL)) stop('Grades de amostras divergentes')
       z[,`:=`(fold=fi,gene_set=gene_set,profile=profile,scenario=scenario,
                method=method,nuclear_error_pp=100*(estimate-true_nuclear),
                rna_error_pp=100*(estimate-true_rna))]
       fwrite(z,file.path(output_run,paste0('06g_',method,'_estimates.tsv')),sep='\t')
     })
   }
   if ('MuSiC2' %in% completed_methods) {
     # Ausencia de campos de convergencia e tratada como NA, nao como sucesso.
     cv<-fit2$convergence
     converged<-if(is.null(cv) || all(is.na(cv))) NA else
       all(as.logical(cv),na.rm=TRUE)
     n_nonconv<-if(is.null(fit2$id.not.converge)) NA_integer_ else
       length(unique(fit2$id.not.converge))
     n_iter<-if(is.null(fit2$n.iter)) NA_integer_ else as.integer(fit2$n.iter[1])
     removed<-if(is.null(fit2$DE.genes)) character() else
       unique(as.character(unlist(fit2$DE.genes,use.names=FALSE)))
     removed<-removed[!is.na(removed)&nzchar(removed)]
     fwrite(data.table(gene_id=removed),file.path(output_run,'06g_removed_by_MuSiC2.tsv'),sep='\t')
     if (has_DE) {
       actual<-intersect(de$genes,rownames(sc))
       found<-intersect(removed,rownames(sc))
       tp<-length(intersect(actual,found))
       fwrite(data.table(n_true_DE=length(actual),n_removed=length(found),
         true_positive=tp,recall=if(length(actual))tp/length(actual) else NA_real_,
         precision=if(length(found))tp/length(found) else NA_real_),
         file.path(output_run,'06g_DE_detection.tsv'),sep='\t')
     }
   } else {
     converged<-NA;n_nonconv<-NA_integer_;n_iter<-NA_integer_
     removed<-character()
   }
   if(length(method_errors)) {
     fwrite(data.table(run=code,status='FAIL',error=method_errors),
            file.path(output_run,'ERROR.tsv'),sep='\t')
   } else if(file.exists(file.path(output_run,'ERROR.tsv'))) {
     unlink(file.path(output_run,'ERROR.tsv')) # sucesso apos tentativa anterior
   }
   status<-if(length(method_errors)) 'FAIL' else 'PASS'
   fwrite(data.table(run=code,status=status,n_genes=nrow(sc),
        n_controls=ncol(cmat),n_cases=ncol(amat),
        n_methods_exported=length(completed_methods),
        methods_exported=paste(completed_methods,collapse=';'),
        converged=converged,n_not_converged=n_nonconv,n_iterations=n_iter,
        n_removed_genes=length(removed),
        elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units='secs'))),
        file.path(output_run,'STATUS.tsv'),sep='\t')
   # Apenas criamos DONE quando AMBOS os resultados estao corretos.
   if(status=='PASS') {
     file.copy(file.path(output_run,'STATUS.tsv'),file.path(output_run,'DONE.tsv'),
               overwrite=TRUE)
   } else if(file.exists(file.path(output_run,'DONE.tsv'))) {
     unlink(file.path(output_run,'DONE.tsv'))
   }
   message(code,' : ',status,if(length(method_errors))
           paste0(' : ',paste(method_errors,collapse=' | ')) else '')
   rm(cmat,amat,sc,sce,controls,cases,fitM,fit2);gc(verbose=FALSE)
  }
 }
}

# Registrar TODOS os status antes da consolidacao, inclusive falhas.
status_paths<-list.files(file.path(DIR,'runs'),pattern='^STATUS.tsv$',
                         recursive=TRUE,full.names=TRUE)
if (length(status_paths)) {
  all_status<-rbindlist(lapply(status_paths,fread),fill=TRUE)
  fwrite(all_status,file.path(DIR,'06g_all_run_status.tsv'),sep='\t')
  message('Execucoes registradas: ',nrow(all_status),
          ' | PASS: ',sum(all_status$status=='PASS'),
          ' | FAIL: ',sum(all_status$status=='FAIL'))
}
# Consolidar SOMENTE ajustes completos, sem esconder os que falharam.
paths<-list.files(file.path(DIR,'runs'),pattern='^DONE.tsv$',
                  recursive=TRUE,full.names=TRUE)
if(!length(paths)) stop('Nenhum ajuste completo. Veja 06g_all_run_status.tsv e ERROR.tsv em runs/.')
summary<-rbindlist(lapply(paths,fread),fill=TRUE)
fwrite(summary,file.path(DIR,'06g_run_summary.tsv'),sep='\t')
longfiles<-list.files(file.path(DIR,'runs'),pattern='06g_(MuSiC|MuSiC2)_estimates.tsv$',recursive=TRUE,full.names=TRUE)
if (!length(longfiles)) stop('Ausencia de estimativas em execucoes marcadas como DONE.')
# Evitar agregar arquivos de estimativas parciais de execucoes com status FAIL.
# Somente os arquivos nas pastas com DONE.tsv entram no resumo.
complete_dirs<-normalizePath(dirname(paths),winslash='/',mustWork=TRUE)
longfiles<-longfiles[normalizePath(dirname(longfiles),winslash='/',mustWork=TRUE) %in% complete_dirs]
if (!length(longfiles)) stop('Nenhuma estimativa disponivel em execucoes completas.')
estimates_all<-rbindlist(lapply(longfiles,fread),fill=TRUE)
fwrite(estimates_all,file.path(DIR,'06g_estimates_long.tsv'),sep='\t')
metrics<-estimates_all[,.(n=.N,MAE_nuclear_pp=mean(abs(nuclear_error_pp)),RMSE_nuclear_pp=sqrt(mean(nuclear_error_pp^2)),
 MAE_RNA_pp=mean(abs(rna_error_pp)),bias_nuclear_pp=mean(nuclear_error_pp),
 frac_zero_positive=mean(estimate<=1e-10 & true_nuclear>0)),
 by=.(fold,gene_set,profile,scenario,method,group,cell_type)]
fwrite(metrics,file.path(DIR,'06g_metrics_by_cell.tsv'),sep='\t')
aggregate<-metrics[,.(MAE_nuclear_pp=mean(MAE_nuclear_pp),MAE_RNA_pp=mean(MAE_RNA_pp),
 RMSE_nuclear_pp=sqrt(mean(RMSE_nuclear_pp^2)),frac_zero_positive=mean(frac_zero_positive)),
 by=.(gene_set,scenario,method,group,cell_type)]
fwrite(aggregate,file.path(DIR,'06g_summary_by_scenario.tsv'),sep='\t')
# Differences must be paired at the FOLD (held-out donor) level; individual
# synthetic libraries from the same donor are NOT independent replicates.
pair<-dcast(metrics[,.(fold,gene_set,profile,scenario,group,cell_type,method,MAE_nuclear_pp)],
 fold+gene_set+profile+scenario+group+cell_type~method,value.var='MAE_nuclear_pp')
if(base::all(c('MuSiC','MuSiC2') %in% names(pair)))pair[,delta_MAE_MuSiC2_minus_MuSiC:=MuSiC2-MuSiC]
fwrite(pair,file.path(DIR,'06g_paired_methods.tsv'),sep='\t')
fwrite(data.table(R=R.version.string,MuSiC=as.character(utils::packageVersion('MuSiC')),mode=MODE,top300=WITH_TOP300,seed=SEED,
 n_folds=N_FOLDS,n_per_group=N_PER_GROUP,foldchange=FOLD_CHANGE,DE_fraction=DE_GENE_FRACTION,version='06g_v3_accessor_fix',
 note='Simulação baseada em Mathys; não prova acurácia em bulk real.'),file.path(DIR,'06g_parameters.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(DIR,'06g_sessionInfo.txt'))
cat('\nResultados em: ',DIR,'\nAjustes concluídos: ',nrow(summary),'\n',sep='')
print(aggregate[cell_type %in% c('In','Opc')])
if(any(summary$converged %in% FALSE,na.rm=TRUE))message('ATENÇÃO: MuSiC2 apresentou execuções não convergidas. Não interpretar como validação final.')
