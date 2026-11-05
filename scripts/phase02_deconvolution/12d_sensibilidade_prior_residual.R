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
OUT <- P('results','deconvolution','bMIND_pilot','12d',DS,SOURCE)
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

# SENSIBILIDADE EXPLORATORIA: nao substitui automaticamente MIND original.
# R: inverse-Wishart escalar; V usa var bulk como escala empirica.
# Comparar duas escolhas de nu; nao selecionar pela significancia biologica.
for(pkg in c("MCMCglmm","coda")) if(!requireNamespace(pkg,quietly=TRUE)) stop(pkg," ausente")
OUT <- file.path(OUT,format(Sys.time(),"run_%Y%m%d_%H%M%S"))
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
configs <- data.table(model=c("original","R_nu0002","R_nu002"),nu_R=c(NA,0.002,0.02))
W(configs,"12d_configs.tsv")
seeds <- c(17011L,28019L,39023L)
status <- list(); diag_rows <- list(); estimates <- list(); residual <- list(); compact <- list()
K <- ncol(F); ids <- rownames(F)
flush <- function() {
 if(length(status)) W(rbindlist(status),"12d_runs.tsv")
 if(length(diag_rows)) W(rbindlist(diag_rows),"12d_diagnostics.tsv")
 if(length(estimates)) W(rbindlist(estimates),"12d_expression.tsv")
 if(length(residual)) W(rbindlist(residual),"12d_residual.tsv")
}
W(data.table(status="RUNNING"),"12d_status.tsv")
for(g in rownames(L)[1:3]) for(ci in seq_len(nrow(configs))) {
 model <- configs$model[ci]; chains <- list()
 for(ch in seq_along(seeds)) {
  x <- as.numeric(L[g,]); sample_id <- ids
  dat <- data.frame(F,sample_id=sample_id,x=x,check.names=FALSE)
  prior <- list(B=list(mu=rep(mean(x),K),V=diag(.5,K)),
    G=list(G1=list(nu=50,V=diag(K)*var(x)/sum(colMeans(F)^2))))
  if(!is.na(configs$nu_R[ci])) prior$R <- list(V=var(x),nu=configs$nu_R[ci])
  set.seed(seeds[ch])
  fit <- tryCatch(MCMCglmm::MCMCglmm(
   fixed=as.formula(paste("x ~ -1 +",paste(CELLS,collapse="+"))),
   random=as.formula(paste("~us(",paste(CELLS,collapse="+"),"):sample_id")),
   data=dat,prior=prior,verbose=FALSE,pr=TRUE,nitt=13000,burnin=3000,thin=10),error=function(e)e)
  ok <- !inherits(fit,"error")
  status[[length(status)+1L]] <- data.table(gene=g,model=model,chain=ch,seed=seeds[ch],
   status=if(ok) "RETURNED" else "FAILED",message=if(ok) "" else conditionMessage(fit),
   residual_V=if(model=="original") NA_real_ else var(x),nu_R=configs$nu_R[ci])
  if(!ok) {flush();next}
  # Expressao por tipo/amostra = efeito fixo + efeito aleatorio correspondente.
  # Resolver nomes explicitamente; nao assumir ordem alfabetica dos IDs.
  sol <- as.matrix(fit$Sol); out <- matrix(NA_real_,nrow(sol),K*length(ids))
  labels <- character(ncol(out));j <- 0L
  for(cell in CELLS) for(id in ids) {
   j <- j+1L; rn <- paste0(cell,".sample_id.",id)
   if(!all(c(cell,rn)%in%colnames(sol))) stop("Nome de coeficiente inesperado: ",rn)
   out[,j] <- sol[,cell]+sol[,rn];labels[j] <- paste(cell,id,sep="|")
  }
  colnames(out) <- labels
  estimates[[length(estimates)+1L]] <- data.table(gene=g,model=model,chain=ch,
   cell_type=rep(CELLS,each=length(ids)),sample_id=rep(ids,K),
   mean=colMeans(out),sd=apply(out,2,sd))
  u <- fit$VCV[,"units"]
  residual[[length(residual)+1L]] <- data.table(gene=g,model=model,chain=ch,
   minimum=min(u),median=median(u),maximum=max(u),first_quarter=mean(head(u,250)),last_quarter=mean(tail(u,250)))
  # Guardar apenas VCV e efeitos fixos: compactos e suficientes para revisar residuos.
  compact[[paste(g,model,ch,sep="_")]] <- list(VCV=fit$VCV,fixed=fit$Sol[,CELLS,drop=FALSE],seed=seeds[ch],prior=prior)
  chains[[as.character(ch)]] <- list(VCV=fit$VCV,CTS=coda::mcmc(out,start=3001,thin=10))
  flush()
 }
 if(length(chains)==3L) for(component in c("VCV","CTS")) {
  ml <- coda::mcmc.list(lapply(chains,function(z) z[[component]]))
  rh <- tryCatch(coda::gelman.diag(ml,multivariate=FALSE,autoburnin=FALSE)$psrf[,1],error=function(e)rep(NA_real_,ncol(ml[[1]])))
  es <- coda::effectiveSize(ml)
  diag_rows[[length(diag_rows)+1L]] <- data.table(gene=g,model=model,component=component,
    parameter=colnames(ml[[1]]),rhat_classic=rh,ess_combined=as.numeric(es))
 }
 flush();saveRDS(compact,file.path(OUT,"12d_compact_chains.rds"))
 cat("12d:",g,model,"concluido;",length(chains),"cadeias retornadas\n")
}
writeLines(capture.output(sessionInfo()),file.path(OUT,"12d_sessionInfo.txt"))
W(data.table(status="SENSITIVITY_COMPLETED_REVIEW_REQUIRED",note="Empirical residual scales; no automatic model selection or convergence claim"),"12d_status.tsv")
cat("Saidas:",OUT,"\n")
