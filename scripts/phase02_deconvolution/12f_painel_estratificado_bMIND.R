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
OUT <- P('results','deconvolution','bMIND_pilot','12f',DS,SOURCE)
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
# Painel TECNICO: proposta predefinida, sem usar grupos ou DEGs.
# Elegibilidade: CPM >=1 em pelo menos 20% das amostras; variancia logCPM >0.
# Calculo em blocos para evitar densificar toda a matriz de uma vez.
metrics <- list()
for(start in seq.int(1L,nrow(counts),by=1000L)) {
 idx <- start:min(start+999L,nrow(counts))
 cp <- sweep(as.matrix(counts[idx,,drop=FALSE]),2,lib,"/")*1e6
 if(any(!is.finite(cp))||any(cp<0)) stop("Counts invalidos")
 lc <- log2(cp+1)
 metrics[[length(metrics)+1L]] <- data.table(gene_id=rownames(counts)[idx],
  mean_cpm=rowMeans(cp),detection=rowMeans(cp>=1),variance_logCPM=apply(lc,1,var))
}
universe <- rbindlist(metrics)
universe[,eligible:=detection>=0.20 & is.finite(variance_logCPM) & variance_logCPM>0]
# Excluir os tres genes usados na confirmacao anterior; manter registro.
universe[,previous_confirmation:=gene_id %in% c("2670","6506","3798")]
pool <- copy(universe[eligible & !previous_confirmation])
if(nrow(pool)<27L) stop("Universo elegivel pequeno demais para estratificacao")
# Tercis por ranks; empates desemp atados pelo ID, explicitamente reproduzivel.
setorder(pool,mean_cpm,gene_id)
pool[,expression_bin:=pmin(3L,ceiling(seq_len(.N)*3/.N))]
# Variabilidade estratificada DENTRO de cada faixa de expressao.
setorder(pool,expression_bin,variance_logCPM,gene_id)
pool[,variability_bin:=pmin(3L,ceiling(seq_len(.N)*3/.N)),by=expression_bin]
if(nrow(unique(pool[,.(expression_bin,variability_bin)]))!=9L) stop("Estratos incompletos")
set.seed(20261024)
panel <- pool[, .SD[sample.int(.N,1L)],by=.(expression_bin,variability_bin)]
setorder(panel,expression_bin,variability_bin)
genes <- panel$gene_id
NGENES <- length(genes)
W(universe,"12f_selection_universe.tsv")
W(pool,"12f_selection_strata.tsv")
W(panel,"12f_stratified_panel.tsv")
b <- as.matrix(counts[genes,,drop=FALSE])
L <- log2(sweep(b,2,lib,"/")*1e6+1)
source_audit <- data.table(dataset=DS,source=SOURCE,n_samples=ncol(L),n_genes=nrow(L),
                           n_types=ncol(F),mean_fraction=NA_real_,run_requested=RUN)
W(source_audit,'12_input_summary.tsv')
W(data.table(gene_id=genes,panel_order=seq_along(genes),
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
for(pkg in c("MCMCglmm","coda","posterior")) if(!requireNamespace(pkg,quietly=TRUE)) stop(pkg," ausente")
OUT <- file.path(OUT,format(Sys.time(),"run_%Y%m%d_%H%M%S"))
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
configs <- data.table(model="R_nu002",nu_R=0.02,nitt=50000L,burnin=10000L,thin=10L)
W(configs,"12f_configs.tsv")
W(data.table(chain=seq_along(seeds_placeholder <- c(17011L,28019L,39023L,47017L)),seed=seeds_placeholder),"12f_seeds.tsv")
seeds <- c(17011L,28019L,39023L,47017L)
if(toupper(Sys.getenv("BMIND12F_RUN","NO"))!="YES") {
 W(data.table(status="PANEL_ONLY",note="Panel selected; MCMC not run. BMIND12F_RUN=YES to run."),"12f_status.tsv")
 cat("Painel selecionado. Saidas:",OUT,"\n")
 quit(save="no",status=0)
}
status <- list(); diag_rows <- list(); estimates <- list(); residual <- list(); compact <- list()
K <- ncol(F); ids <- rownames(F)
flush <- function() {
 if(length(status)) W(rbindlist(status),"12f_runs.tsv")
 if(length(diag_rows)) W(rbindlist(diag_rows),"12f_diagnostics.tsv")
 if(length(estimates)) W(rbindlist(estimates),"12f_expression.tsv")
 if(length(residual)) W(rbindlist(residual),"12f_residual.tsv")
}
W(data.table(status="RUNNING"),"12f_status.tsv")
for(g in rownames(L)) for(ci in seq_len(nrow(configs))) {
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
   data=dat,prior=prior,verbose=FALSE,pr=TRUE,nitt=50000,burnin=10000,thin=10),error=function(e)e)
  ok <- !inherits(fit,"error")
  status[[length(status)+1L]] <- data.table(gene=g,model=model,chain=ch,seed=seeds[ch],
   status=if(ok) "RETURNED" else "FAILED",message=if(ok) "" else gsub("[\r\n\t]+"," ",conditionMessage(fit)),
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
   minimum=min(u),median=median(u),maximum=max(u),first_quarter=mean(head(u,length(u)%/%4)),last_quarter=mean(tail(u,length(u)%/%4)))
  # Guardar apenas VCV e efeitos fixos: compactos e suficientes para revisar residuos.
  compact[[paste(g,model,ch,sep="_")]] <- list(VCV=fit$VCV,fixed=fit$Sol[,CELLS,drop=FALSE],seed=seeds[ch],prior=prior)
  chains[[as.character(ch)]] <- list(VCV=fit$VCV,CTS=coda::mcmc(out,start=10001,thin=10))
  flush()
 }
 if(length(chains)==4L) for(component in c("VCV","CTS")) {
  mats <- lapply(chains,function(z) as.matrix(z[[component]]))
  stopifnot(all(vapply(mats,function(z) identical(colnames(z),colnames(mats[[1]])),logical(1))))
  dd <- lapply(seq_len(ncol(mats[[1]])),function(j) {
   draws <- do.call(cbind,lapply(mats,function(z) z[,j]))
   sd_post <- sd(as.vector(draws)); mc <- posterior::mcse_mean(draws)
   data.table(gene=g,model=model,component=component,parameter=colnames(mats[[1]])[j],
    posterior_mean=mean(draws),posterior_sd=sd_post,
    rhat=posterior::rhat(draws),ess_bulk=posterior::ess_bulk(draws),
    ess_tail=posterior::ess_tail(draws),mcse_mean=mc,mcse_over_sd=mc/sd_post)
  })
  diag_rows[[length(diag_rows)+1L]] <- rbindlist(dd)

 }
 flush();saveRDS(compact,file.path(OUT,"12f_compact_chains.rds"))
 cat("12f:",g,model,"concluido;",length(chains),"cadeias retornadas\n")
}
writeLines(capture.output(sessionInfo()),file.path(OUT,"12f_sessionInfo.txt"))
W(data.table(status="SENSITIVITY_COMPLETED_REVIEW_REQUIRED",note="Empirical residual scales; no automatic model selection or convergence claim"),"12f_status.tsv")
cat("Saidas:",OUT,"\n")

# Criteres de triagem predefinidos: nao equivalem a validacao biologica.
if(length(diag_rows)) {
 d <- rbindlist(diag_rows)
 d[,flag:=!is.finite(rhat)|!is.finite(ess_bulk)|!is.finite(ess_tail)|
    !is.finite(mcse_over_sd)|rhat>=1.01|ess_bulk<400|ess_tail<400|mcse_over_sd>0.05]
 W(d[,.(n_parameters=.N,n_flagged=sum(flag)),by=.(gene,component)],"12f_summary.tsv")
 W(d[flag==TRUE],"12f_flagged.tsv")
}
