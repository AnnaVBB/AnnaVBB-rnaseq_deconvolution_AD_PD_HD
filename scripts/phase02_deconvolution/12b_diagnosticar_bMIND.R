#!/usr/bin/env Rscript
# ==============================================================================
# 12 | bMIND — PREPARACAO E PILOTO TECNICO COM PROPORCOES AUDITADAS
# CORRECOES: parser R na selecao de fracoes e verificacao de IDs duplicados.
# Metodologia, parametros e modo de auditoria padrao preservados.
# TCC AD / PD / HD. MIND::bMIND NAO E ALGORITMO DE PROPORCOES.
# Ele recebe bulk log2(CPM+1) (genes x amostras) e frac (amostras x tipos)
# e estima expressao especifica de tipo (gene x tipo x amostra).
#
# MODO SEGURO PADRAO: apenas auditoria de preparo (sem MCMC).
#   Rscript scripts/phase02_deconvolution/12_piloto_bMIND.R
# PARA RODAR O MCMC APOS A REVISAO DAS FRACOES E DOS TIPOS <5%:
#   BMIND12_RUN=YES BMIND12_ALLOW_LOW=YES Rscript .../12_piloto_bMIND.R
# O override nao corrige o problema de baixa abundancia; registra a ressalva.
#
# Parametros:
#   BMIND12_DATASET=ROSMAP (default: apenas piloto; AD ainda nao definido)
#   BMIND12_FRACTIONS=dtangle_pooled_counts (default)
#       opcoes: dtangle_pooled_counts, dtangle_mean_nucleus_cpm,
#               BisqueRNA, MuSiC2
#   BMIND12_NGENES=30 (genes tecnicos escolhidos sem grupos/valores-p)
#   BMIND12_NCORE=1
#   BMIND12_RUN=NO/YES
#   BMIND12_ALLOW_LOW=NO/YES
# ==============================================================================
options(stringsAsFactors=FALSE)
REQ <- c('data.table','Matrix')
miss <- REQ[!vapply(REQ,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss))stop('Pacotes ausentes: ',paste(miss,collapse=', '),' (renv).')
suppressPackageStartupMessages({library(data.table);library(Matrix)})
ROOT <- normalizePath(Sys.getenv('DECONV_ROOT','.'),mustWork=TRUE)
P <- function(...) file.path(ROOT,...)
DS <- Sys.getenv('BMIND12_DATASET','ROSMAP')
if(!DS%in%c('ROSMAP','GSE125583','GSE64810','GSE68719'))stop('Coorte nao reconhecida')
SOURCE <- Sys.getenv('BMIND12_FRACTIONS','dtangle_pooled_counts')
if(!SOURCE%in%c('dtangle_pooled_counts','dtangle_mean_nucleus_cpm','BisqueRNA','MuSiC2'))stop('Fonte de proporcoes desconhecida')
RUN <- toupper(Sys.getenv('BMIND12_RUN','NO'))=='YES'
LOW_OK <- toupper(Sys.getenv('BMIND12_ALLOW_LOW','NO'))=='YES'
NGENES <- suppressWarnings(as.integer(Sys.getenv('BMIND12_NGENES','30')))
NCORE <- suppressWarnings(as.integer(Sys.getenv('BMIND12_NCORE','1')))
if(is.na(NGENES)||NGENES<5||NGENES>100)stop('BMIND12_NGENES entre 5 e 100 para piloto')
if(is.na(NCORE)||NCORE<1||NCORE>4)stop('BMIND12_NCORE entre 1 e 4')
CELLS <- c('Ex','In','Ast','Oli','Opc','Mic','End','Per')
OUT <- P('results','deconvolution','bMIND_pilot','12b',DS,SOURCE)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
W <- function(z,n){fwrite(z,file.path(OUT,n),sep='\t',quote=FALSE,na='NA');invisible(z)}
ERR <- function(...) stop(paste0(...,collapse=''))
# Conteudos originais da etapa 05 para a matriz de expressao.
cf <- P('data','processed','deconvolution',DS,'bulk_counts_reference_intersection.rds')
lf <- P('data','processed','deconvolution',DS,'library_sizes.tsv')
mf <- P('data','processed','deconvolution',DS,'metadata_deconvolution.tsv')
# IMPORTANT: em R, `else` deve fazer parte da mesma expressao do `if`.
# Blocos explicitos evitam o erro de parse: unexpected 'else'.
ff <- if (SOURCE %in% c('dtangle_pooled_counts', 'dtangle_mean_nucleus_cpm')) {
  P('results', 'deconvolution', 'dtangle',
    '08e_replicacao_orientador', '08e_proportions_long_all.tsv')
} else if (SOURCE == 'BisqueRNA') {
  P('results', 'deconvolution', 'BisqueRNA', '07_proportions_long_all.tsv')
} else {
  P('results', 'deconvolution', 'MuSiC2', '06_proportions_long_all.tsv')
}
for(path in c(cf,lf,mf,ff))if(!file.exists(path))ERR('Arquivo ausente: ',path)
meta <- fread(mf); libs <- fread(lf); fr <- fread(ff)
for(nm in c('sample_id','group'))if(!nm%in%names(meta))ERR('Metadata missing ',nm)
for(nm in c('sample_id','library_size_full_input'))if(!nm%in%names(libs))ERR('Library sizes missing ',nm)
for(nm in c('dataset','sample_id','group','cell_type','proportion'))if(!nm%in%names(fr))ERR('Fractions missing ',nm)
# Verifica os vetores de IDs diretamente. O segundo argumento posicional de
# anyDuplicated() e 'incomparables', nao o nome da coluna (erro anterior).
if(anyDuplicated(meta$sample_id)||anyDuplicated(libs$sample_id))ERR('Duplicate sample IDs')
fr <- fr[dataset==DS]
if(SOURCE=='MuSiC2') {
 if(!all(c('marker_set','reference_version')%in%names(fr)))ERR('MuSiC2 source labels missing')
 fr <- fr[marker_set=='top300'&reference_version=='main24']
} else if(SOURCE=='BisqueRNA') {
 if('marker_set'%in%names(fr))fr <- fr[marker_set=='top300']
 if('reference_version'%in%names(fr))fr <- fr[reference_version=='main24']
} else {
 if(!'reference_mode'%in%names(fr))ERR('dtangle reference_mode missing')
 mode <- if(SOURCE=='dtangle_pooled_counts')'pooled_counts' else 'mean_nucleus_cpm'
 fr <- fr[reference_mode==mode]
}
if(!nrow(fr)||anyDuplicated(fr,by=c('sample_id','cell_type')))ERR('Fractions missing or duplicated')
if(!setequal(fr$sample_id,meta$sample_id)||!setequal(fr$cell_type,CELLS))ERR('Fractions do not match 05 metadata/8 types')
checks <- merge(fr,meta[,.(sample_id,correct_group=group)],by='sample_id',all.x=TRUE)
if(anyNA(checks$correct_group)||any(checks$group!=checks$correct_group))ERR('Group mismatch')
wide <- dcast(fr,sample_id~cell_type,value.var='proportion')
wide <- wide[match(meta$sample_id,wide$sample_id)]
if(!identical(as.character(wide$sample_id),as.character(meta$sample_id)))ERR('Sample order mismatch')
F <- as.matrix(wide[,..CELLS]);storage.mode(F)<-'double'
rownames(F) <- wide$sample_id
if(any(!is.finite(F))||any(F<0)||any(F>1)||max(abs(rowSums(F)-1))>1e-6)ERR('Invalid fractions')
counts <- readRDS(cf)
if(is.null(rownames(counts))||is.null(colnames(counts)))ERR('Counts without dimnames')
if(anyDuplicated(rownames(counts))||anyDuplicated(colnames(counts)))ERR('Duplicate gene or sample IDs')
if(!setequal(colnames(counts),meta$sample_id))ERR('Counts not aligned with metadata')
counts <- counts[,meta$sample_id,drop=FALSE]
lib <- as.numeric(libs$library_size_full_input[match(colnames(counts),libs$sample_id)])
if(any(!is.finite(lib))||any(lib<=0))ERR('Invalid complete library sizes')
if(any(as.numeric(Matrix::colSums(counts))>lib+1e-5))ERR('Full libraries below intersected counts')
if(inherits(counts,'sparseMatrix')) {
 if(any(counts@x<0))ERR('Negative counts')
} else if(any(counts<0)) ERR('Negative counts')
if(!all(colnames(counts)==rownames(F)))ERR('Final sample order mismatch')
# Escolha TECNICA de genes sem usar grupo, DEGs, FDR ou algoritmo como verdade.
# Top-expressos dentre genes nao constantes; podem ser substituidos em futuro
# estudo preespecificado, mas este painel nao fundamenta interpretacao biologica.
mean_cpm <- as.numeric(Matrix::rowMeans(t(t(counts)/lib)*1e6))
ix <- which(is.finite(mean_cpm)&mean_cpm>1)
if(length(ix)<NGENES)ERR('Insufficient detected genes for pilot')
ix <- ix[order(-mean_cpm[ix],rownames(counts)[ix])][seq_len(NGENES)]
genes <- rownames(counts)[ix]
b <- as.matrix(counts[genes,,drop=FALSE])
L <- log2(sweep(b,2,lib,'/')*1e6+1)
if(any(!is.finite(L))||max(L)>50)ERR('Bulk log2CPM transformation suspicious')
source_audit <- data.table(dataset=DS,source=SOURCE,n_samples=ncol(L),n_genes=nrow(L),
                           n_types=ncol(F),mean_fraction=NA_real_,run_requested=RUN)
W(source_audit,'12_input_summary.tsv')
W(data.table(gene_id=genes,mean_expression_rank=seq_along(genes),
    purpose='TECHNICAL PILOT ONLY; not a disease/DEG-selected biological panel'),
   '12_gene_panel.tsv')
quality <- data.table(cell_type=CELLS,mean_fraction=colMeans(F),median_fraction=apply(F,2,median),
                      sd_pp=100*apply(F,2,sd),fraction_exact_zero=colMeans(F<=1e-12))
quality[,low_mean_below_5pct:=mean_fraction<0.05]
quality[,high_zero_fraction:=fraction_exact_zero>0.9]
W(quality,'12_fraction_readiness.tsv')
W(data.table(input_name=c('bulk_counts','libraries','metadata','fractions'),
             path=c(cf,lf,mf,ff)),'12_input_manifest.tsv')
# MIND source-code constraints: names without punctuation / periods; maintain ID map.
newid <- paste0('s',seq_len(ncol(L)))
W(data.table(original_sample_id=colnames(L),bmind_safe_sample_id=newid),
  '12_sample_id_mapping.tsv')
colnames(L) <- rownames(F) <- newid

# 12b: diagnostico por gene, usando a funcao interna instalada.
# Nao modifica namespace/pacote. Guarda Sol e VCV antes do resumo do MIND.
for(pkg in c("MIND","coda")) if(!requireNamespace(pkg,quietly=TRUE)) stop(pkg," ausente")
fun <- getFromNamespace("lme_mc2","MIND")
writeLines(deparse(fun),file.path(OUT,"12b_original_lme_mc2.R"))
bb <- as.list(body(fun))
seed_idx <- which(vapply(bb,function(z) is.call(z)&&identical(z[[1]],as.name("set.seed")),logical(1)))
fit_idx <- which(vapply(bb,function(z) is.call(z)&&as.character(z[[1]])[1]%in%c("=","<-")&&identical(z[[2]],as.name("lme2")),logical(1)))
if(length(seed_idx)!=1||length(fit_idx)!=1) stop("Funcao interna diferente: revisar 12b_original_lme_mc2.R")
# Remove seed fixa e retorna o objeto MCMCglmm completo imediatamente apos ajuste.
bb[[seed_idx]] <- quote(NULL)
bb <- bb[seq_len(fit_idx)]
bb[[length(bb)+1L]] <- quote(return(lme2))
body(fun) <- as.call(bb)
writeLines(deparse(fun),file.path(OUT,"12b_instrumented_lme_mc2.R"))
writeLines(capture.output(sessionInfo()),file.path(OUT,"12b_sessionInfo.txt"))
# Mesmos priors padrao do bmind: media bulk, V_fe=.5I, nu=50,
# V_re=I*var(bulk)/sum(mean(frac)^2). Painel tecnico, sem selecao por grupo.
seeds <- c(17011L,28019L,39023L)
schedules <- list(short=c(1300L,300L,1L),long=c(13000L,3000L,10L))
records <- list(); stability <- list(); failures <- list()
W(data.table(status="RUNNING",note="Current run in progress"),"12b_status.tsv")
for(g in rownames(L)) {
 for(stage in names(schedules)) {
  pars <- schedules[[stage]]; fits <- list()
  for(ch in seq_along(seeds)) {
   set.seed(seeds[ch])
   fit <- tryCatch(fun(x=L[g,],W=F,sample_id=rownames(F),mu=rep(mean(L[g,]),ncol(F)),
     V_fe=diag(.5,ncol(F)),V_re=diag(ncol(F))*var(L[g,])/sum(colMeans(F)^2),
     nu=50,nitt=pars[1],burnin=pars[2],thin=pars[3]), error=function(e) e)
   if(inherits(fit,"error")) {
    failures[[length(failures)+1L]] <- data.table(gene=g,stage=stage,chain=ch,seed=seeds[ch],message=conditionMessage(fit))
    W(rbindlist(failures),"12b_failures.tsv")
    cat("12b falha registrada:",g,stage,ch,conditionMessage(fit),"\n")
    next
   }
   saveRDS(list(gene=g,stage=stage,seed=seeds[ch],Sol=fit$Sol,VCV=fit$VCV,
     parameters=pars,samples=rownames(F),types=colnames(F)),
     file.path(OUT,paste0("12b_chain_",g,"_",stage,"_",ch,".rds")))
   fits[[as.character(ch)]] <- fit
  }
  if(length(fits)!=length(seeds)) {
   cat("12b diagnostico conjunto omitido: cadeias incompletas para",g,stage,"\n")
   next
  }
  for(component in c("Sol","VCV")) {
   chains <- lapply(fits,function(f) f[[component]])
   if(!all(vapply(chains,function(z) identical(colnames(z),colnames(chains[[1]])),logical(1)))) stop("Chain names mismatch")
   ml <- coda::mcmc.list(chains)
   # Rhat classico, nao rank-normalized. Sem declaracao automatica de convergencia.
   rh <- tryCatch(coda::gelman.diag(ml,autoburnin=FALSE,multivariate=FALSE)$psrf[,1],error=function(e) rep(NA_real_,ncol(chains[[1]])))
   ess <- coda::effectiveSize(ml)
   records[[length(records)+1L]] <- data.table(gene=g,stage=stage,component=component,
     parameter=colnames(chains[[1]]),rhat_classic=rh,ess_combined=as.numeric(ess))
   for(ch in seq_along(chains)) stability[[length(stability)+1L]] <- data.table(
     gene=g,stage=stage,component=component,chain=ch,parameter=colnames(chains[[ch]]),
     posterior_mean=colMeans(chains[[ch]]),posterior_sd=apply(chains[[ch]],2,sd))
   # Graficos dos efeitos fixos e covariancias; todas as cadeias ficam no RDS.
   selected <- seq_len(min(if(component=="Sol") ncol(F) else 9L,ncol(chains[[1]])))
   pdf(file.path(OUT,paste0("12b_trace_",g,"_",stage,"_",component,".pdf")))
   plot(coda::mcmc.list(lapply(chains,function(z) z[,selected,drop=FALSE])))
   dev.off()
  }
  W(rbindlist(records),"12b_diagnostics.tsv")
  W(rbindlist(stability),"12b_posterior_summary.tsv")
  cat("12b concluido:",g,stage,"\n")
 }
}
W(data.table(status=if(length(failures)) "COMPLETED_WITH_FAILED_CHAINS_REVIEW_REQUIRED" else "DIAGNOSTICS_COMPLETED_REVIEW_REQUIRED",note="Classic Rhat/ESS; inspect chains; no biological DE claims"),"12b_status.tsv")

if(!length(failures)) W(data.table(gene=character(),stage=character(),chain=integer(),seed=integer(),message=character()),"12b_failures.tsv")
