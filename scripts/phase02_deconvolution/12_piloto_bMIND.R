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
OUT <- P('results','deconvolution','bMIND_pilot','12',DS,SOURCE)
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
set.seed(20261023)
if(!RUN) {
 W(data.table(status='PREFLIGHT_ONLY',message=paste('bMIND not executed. ',
   sum(quality$low_mean_below_5pct),' cell types have mean <5%; ',
   'review fractions before opt-in.')), '12_status.tsv')
 cat('12 PREFLIGHT PASS | ',DS,' / ',SOURCE,'; bMIND NAO executado. Saidas: ',OUT,'\n',sep='')
 quit(save='no',status=0)
}
if(any(quality$low_mean_below_5pct)||any(quality$high_zero_fraction)) {
 if(!LOW_OK)ERR('Low-abundance/zero-fraction types found. ',
    'bMIND not executed. Use BMIND12_ALLOW_LOW=YES ONLY after methodological review.')
}
if(!requireNamespace('MIND',quietly=TRUE))ERR('MIND package unavailable. Install in renv separately, no automatic updates.')
# Pilot WITHOUT supplied informative scRNA priors: documented default MIND
# priors computed by bMIND, not proof that priors are correct.
start <- proc.time()[['elapsed']]
fit <- tryCatch(MIND::bMIND(bulk=L,frac=F,ncore=NCORE,nitt=1300,burnin=300,thin=1),
                error=function(e)e)
if(inherits(fit,'error')) {
 W(data.table(status='FAIL_MCMC',message=conditionMessage(fit)),'12_status.tsv')
 stop('bMIND failed: ',conditionMessage(fit))
}
if(!is.list(fit)||is.null(fit$A)||is.null(fit$SE))ERR('MIND::bMIND output missing A/SE')
a <- fit$A; se <- fit$SE
if(length(dim(a))!=3L||length(dim(se))!=3L||!identical(dim(a),dim(se)))ERR('A / SE unexpected dimensions')
if(dim(a)[2]!=length(CELLS)||dim(a)[3]!=nrow(F))ERR('A mismatch types or samples')
if(any(!is.finite(a))||any(!is.finite(se))||any(se<0))ERR('A/SE nonfinite or negative SE')
# Critically, MCMC return A/SE is not itself a convergence diagnostic.
saveRDS(list(A=a,SE=se,Sigma_c=fit$Sigma_c,mu=fit$mu,
             fraction_source=SOURCE,dataset=DS,parameters=list(ncore=NCORE,nitt=1300,
             burnin=300,thin=1),original_sample_id=colnames(counts)),
        file.path(OUT,'12_bMIND_pilot_result.rds'),compress=TRUE)
W(data.table(status='MCMC_COMPLETED_NOT_CONVERGENCE_VALIDATED',
             dataset=DS,source=SOURCE,elapsed_seconds=proc.time()[['elapsed']]-start,
             n_input_genes=nrow(L),n_output_genes=dim(a)[1],
             n_cell_types=dim(a)[2],n_samples=dim(a)[3],
             note='No claim of MCMC convergence or biological CTS-DE'),
   '12_status.tsv')
writeLines(capture.output(sessionInfo()),file.path(OUT,'12_sessionInfo.txt'))
cat('12 | MCMC pilot returned A/SE. Review low types, output genes, MCMC diagnostics; NO DE claims.\n')
