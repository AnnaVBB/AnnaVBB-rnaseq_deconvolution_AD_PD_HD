#!/usr/bin/env Rscript
# ==============================================================================
# 11 | ANALISE DE PROPORCOES CASO-CONTROLE COM COVARIAVEIS APROVADAS
# TCC AD / PD / HD. Entradas prontas: 09 + auditoria 10.
#
# DESENHO
#   - Modelo PRINCIPAL exploratorio: OLS em 100*fracao com SE robusto HC3.
#     Reporta delta ajustado caso-controle em pontos percentuais (pp).
#     Nao modela explicitamente o simplex; efeitos relativos, NAO n celular.
#   - SENSIBILIDADE composicional: CLR(log(fracao+eps)) com dois epsilons.
#     Coeficiente e LOG-RATIO, nao pontos percentuais. Zeros impedem inferencia
#     robusta em tipos com colapso; nao tratamos esses zeros como ausencia.
#   - Correcoes BH separadas por coorte/metodo/escala.
#   - NUNCA combina HD e PD em uma unica analise (44 controles em comum).
#   - Covariaveis nao sao selecionadas por valor-p; exige PLANO APROVADO.
#   - Amostras completas: as MESMAS para os oito tipos e tres metodos da coorte.
#
# USO: Rscript scripts/phase02_deconvolution/11_modelos_ajustados_deconvolucao.R
# Configuracao obrigatoria para inferencia:
#   data/processed/deconvolution/10_approved_model_plan.tsv
# Colunas:
#   dataset,approved,numeric_covariates,categorical_covariates,
#   cluster_column,justification,provenance,notes
#   Separe covariaveis por virgula; nunca use variaveis pos-desfecho sem desenho.
# SAIDAS: results/deconvolution/Integrated_Comparison/11_adjusted/
# ==============================================================================
options(stringsAsFactors=FALSE)
REQ <- c('data.table','ggplot2')
missing <- REQ[!vapply(REQ,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Pacotes ausentes: ',paste(missing,collapse=','))
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
ROOT <- normalizePath(Sys.getenv('DECONV_ROOT','.'),mustWork=TRUE)
P <- function(...) file.path(ROOT,...)
OUT <- P('results','deconvolution','Integrated_Comparison','11_adjusted')
PLOTS <- file.path(OUT,'plots'); dir.create(PLOTS,recursive=TRUE,showWarnings=FALSE)
INPUT <- P('results','deconvolution','Integrated_Comparison','09','09_proportions_long_harmonized.tsv')
PLAN <- P('data','processed','deconvolution','10_approved_model_plan.tsv')
META_DIR <- P('results','deconvolution','Integrated_Comparison','10_metadata')
DATASETS <- c('ROSMAP','GSE125583','GSE64810','GSE68719')
DISEASE <- c(ROSMAP='AD',GSE125583='AD',GSE64810='HD',GSE68719='PD')
METHODS <- c('MuSiC2','BisqueRNA','dtangle')
CELLS <- c('Ex','In','Ast','Oli','Opc','Mic','End','Per')
ZERO_WARN <- 0.90; SD_MIN_PP <- 0.20; CLR_EPS <- c(1e-4,1e-5)
WR <- function(x,name) {fwrite(x,file.path(OUT,name),sep='\t',quote=FALSE,na='NA');invisible(x)}
blank <- function(z) is.na(z)|trimws(as.character(z))==''
fields <- function(s) {if(length(s)==0||is.na(s)||!nzchar(trimws(s)))return(character()); unique(trimws(strsplit(s,',',fixed=TRUE)[[1]]))}
rows <- list(); qc <- list(); cov_bal <- list(); status <- list(); samples <- list()
logstat <- function(ds,st,msg){status[[length(status)+1L]]<<-data.table(dataset=ds,status=st,message=msg);cat(ds,': ',st,' - ',msg,'\n',sep='')}
if (!file.exists(INPUT)) stop('Rode primeiro o 09: ',INPUT)
D <- fread(INPUT,na.strings=c('NA','','NaN'))
required <- c('dataset','method','cell_type','sample_id','group','proportion')
if(!all(required%in%names(D))) stop('09_long incompleto: ',paste(setdiff(required,names(D)),collapse=','))
if(anyDuplicated(D,by=c('dataset','method','cell_type','sample_id')))stop('09_long possui chaves duplicadas')
if(anyNA(D$proportion)||any(!is.finite(D$proportion))||any(D$proportion<0|D$proportion>1))stop('Proporcoes fora de [0,1]')
if(!setequal(unique(D$dataset),DATASETS)||!setequal(unique(D$method),METHODS)||!setequal(unique(D$cell_type),CELLS)) stop('Conjuntos inesperados no 09')
chk <- D[,.(n_ct=.N,total=sum(proportion)),by=.(dataset,method,sample_id)]
if(any(chk$n_ct!=8L)|any(abs(chk$total-1)>1e-6))stop('Deconvolucao 09 incompleta ou soma diferente de 1')
if(!file.exists(PLAN)) {
 logstat('ALL','PENDING_APPROVAL',paste('Crie um plano aprovado em ',PLAN, '. Modelo nao foi executado.',sep=''))
 WR(rbindlist(status),'11_status.tsv')
 WR(data.table(note='No inference: plan unavailable. Use 10_MODEL_PLAN_TEMPLATE.tsv'),'11_README_STATUS.tsv')
 quit(save='no',status=0)
}
plan <- fread(PLAN,na.strings=c('NA','NaN'))
needplan <- c('dataset','approved','numeric_covariates','categorical_covariates','cluster_column','justification','provenance')
if(!all(needplan%in%names(plan)))stop('Plano sem: ',paste(setdiff(needplan,names(plan)),collapse=','))
if(nrow(plan)!=4L||!setequal(plan$dataset,DATASETS)||anyDuplicated(plan,'dataset'))stop('Plano exige exatamente quatro coortes distintas')
# Ajuste robusto HC3, nao requer pacote sandwich; erro robusto assintotico.
fit_hc3 <- function(y,X,index) {
 n <- nrow(X); p <- ncol(X)
 if(n <= p+5L || qr(X)$rank != p) stop('Matriz de desenho rank-deficiente ou n insuficiente')
 Q <- chol2inv(chol(crossprod(X)))
 beta <- as.numeric(Q%*%crossprod(X,y)); names(beta) <- colnames(X)
 resid <- as.numeric(y-X%*%beta)
 lev <- rowSums((X %*% Q)*X)
 if(any(lev>=(1-1e-8))) stop('Alavancagem 1: HC3 indefinido')
 meat <- crossprod(X, X*as.numeric((resid/(1-lev))^2))
 vc <- Q%*%meat%*%Q
 se <- sqrt(diag(vc)); df <- n-p
 b <- beta[index]; s <- se[index]
 if(!is.finite(b)||!is.finite(s)||s<=0) stop('Coeficiente ou SE nao estimavel')
 crit <- qt(0.975,df)
 c(effect=b,ci_low=b-crit*s,ci_high=b+crit*s,p=2*pt(-abs(b/s),df),
   se=s,resid_sd=sd(resid),max_leverage=max(lev),df_resid=df)
}
for(ds in DATASETS) {
 pp <- plan[dataset==ds]
 if(toupper(trimws(as.character(pp$approved)))!='YES') {
    logstat(ds,'PENDING_APPROVAL','Plan approved != YES; no adjusted inference')
    next
 }
 num <- fields(pp$numeric_covariates); fac <- fields(pp$categorical_covariates)
 covars <- c(num,fac)
 if(!length(covars)) {logstat(ds,'FAIL_NO_COVARS','An adjusted model requires approved covariates');next}
 if(anyDuplicated(covars)||any(covars%in%c('sample_id','group','dataset','disease','method','cell_type','proportion','case'))) {
    logstat(ds,'FAIL_INVALID_COVARS','Duplicated or reserved covariate names');next
 }
 if(!is.na(pp$cluster_column)&&nzchar(trimws(pp$cluster_column))) {
    logstat(ds,'PENDING_CLUSTER_MODEL','Cluster/repeated measures requested; HC3 is not cluster-robust. Do not silently ignore cluster.');next
 }
 if(blank(pp$justification)||blank(pp$provenance)) {logstat(ds,'FAIL_NO_JUSTIFICATION','Missing covariate justification or provenance');next}
 mf <- file.path(META_DIR,paste0('10_metadata_',ds,'_audited.tsv'))
 if(!file.exists(mf)){logstat(ds,'FAIL_METADATA_MISSING','Run script 10 first');next}
 md <- fread(mf,na.strings=c('NA','','NaN','NULL'))
 if(!all(c('sample_id','group',covars)%in%names(md))) {
    logstat(ds,'FAIL_COVARIATE_NOT_FOUND',paste(setdiff(covars,names(md)),collapse=','));next
 }
 if(anyDuplicated(md,'sample_id'))stop(ds,': duplicated metadata ID')
 md <- md[,c('sample_id','group',covars),with=FALSE]
 for(nm in num) if(!is.numeric(md[[nm]])||any(!is.finite(md[[nm]][!is.na(md[[nm]])]))) {
    stop(ds,': ',nm,' nao numerica/finita. Harmonize explicitamente no metadata_extension.')
 }
 for(nm in fac) md[[nm]] <- factor(trimws(as.character(md[[nm]])))
 # listwise deletion is mandatory and identical for all methods; audit excluded IDs.
 complete <- complete.cases(md)
 complete <- complete & apply(as.data.frame(md[,..covars]),1,function(z)all(!blank(z)))
 lost <- md[!complete,.(sample_id,group)]
 kept <- copy(md[complete])
 samples[[ds]] <- data.table(dataset=ds,sample_id=md$sample_id,group=md$group,
                          used_in_model=complete,reason=ifelse(complete,'COMPLETE','MISSING_COVARIATE'))
 if(nrow(kept)<30||sum(kept$group=='Control')<10||sum(kept$group==DISEASE[[ds]])<10) {
    logstat(ds,'FAIL_SAMPLE_SIZE',paste('after complete case n=',nrow(kept)));next
 }
 # all categorical fields must vary and have >= 2 observed levels (remove unused).
 for(nm in fac) {
    kept[[nm]] <- droplevels(kept[[nm]])
    if(nlevels(kept[[nm]])<2) stop(ds,': categorical covariate ',nm,' has only 1 level')
 }
 for(nm in num) {
    if(sd(kept[[nm]])==0)stop(ds,': numeric covariate ',nm,' constant')
    kept[[nm]] <- as.numeric(kept[[nm]])
 }
 kept[,case := as.integer(group==DISEASE[[ds]])]
 # Form a shared design; do NOT fit each method with a different N.
 f <- reformulate(c('case',covars),response=NULL)
 X <- model.matrix(f,data=as.data.frame(kept))
 if(!'case'%in%colnames(X)||qr(X)$rank!=ncol(X)) {
   logstat(ds,'FAIL_RANK','Group or batch fully confounded, or singular design matrix');next
 }
 if(nrow(X) < max(30L,10L*ncol(X))) {
   logstat(ds,'FAIL_OVERFIT',paste('n=',nrow(X),' < 10 parameters per predictor; p=',ncol(X)));next
 }
 conf <- kept[,.(n=.N),by=.(group)]
 qc[[ds]] <- data.table(dataset=ds,n_total=nrow(md),n_complete=nrow(kept),n_dropped=nrow(lost),
   n_control=conf[group=='Control',n],n_cases=conf[group==DISEASE[[ds]],n],
   model_rank=qr(X)$rank,n_parameters=ncol(X),covariates=paste(covars,collapse=';'),
   condition_number=kappa(X),status='PASS_DESIGN')
 for(nm in covars) {
   a <- kept[case==0][[nm]]; b <- kept[case==1][[nm]]
   if(nm%in%num) {
     for(g in c('Control',DISEASE[[ds]])) {
       z <- kept[group==g][[nm]]
       cov_bal[[length(cov_bal)+1L]] <- data.table(
         dataset=ds,covariate=nm,kind='numeric',group=g,
         n=length(z),mean=mean(z),sd=sd(z),levels=NA_character_)
     }
   } else {
     for(g in c('Control',DISEASE[[ds]])) {
       vals <- kept[group==g][[nm]]
       cov_bal[[length(cov_bal)+1L]] <- data.table(
         dataset=ds,covariate=nm,kind='categorical',group=g,
         n=length(vals),mean=NA_real_,sd=NA_real_,
         levels=paste(paste(names(table(vals)),as.integer(table(vals)),sep='='),collapse=';'))
     }
   }
 }
 # Wide fraction matrix; preserves exact sample order and all eight cell types.
 for(mt in METHODS) {
   long <- D[dataset==ds & method==mt]
   wide <- dcast(long,sample_id~cell_type,value.var='proportion')
   if(!setequal(wide$sample_id,md$sample_id))stop(ds,' ',mt,': sample mapping mismatch')
   wide <- wide[match(kept$sample_id,wide$sample_id)]
   if(!identical(wide$sample_id,kept$sample_id))stop('Sample ordering broken')
   M <- as.matrix(wide[,..CELLS]); storage.mode(M) <- 'double'
   if(anyNA(M)||max(abs(rowSums(M)-1))>1e-6)stop('Invalid fractions')
   zero <- colMeans(M<=1e-12)
   sds <- apply(M,2,sd)*100
   # No silent transformation of fractions for the main outcome.
   for(transform in c('percentage','clr_1e-4','clr_1e-5')) {
     eps <- if(transform=='percentage') NA_real_ else if(transform=='clr_1e-4') 1e-4 else 1e-5
     if(transform=='percentage') Y <- 100*M else {
       logm <- log(M+eps)
       Y <- logm-rowMeans(logm)
     }
     for(ct in CELLS) {
       stat <- 'PASS'
       if(transform=='percentage'&&(zero[[ct]]>=ZERO_WARN||sds[[ct]]<SD_MIN_PP)) stat <- 'SKIP_COLLAPSED_TYPE'
       if(transform!='percentage'&&any(zero>=ZERO_WARN)) stat <- 'SKIP_CLR_STRUCTURAL_ZERO_IN_METHOD'
       est <- c(effect=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p=NA_real_,
                se=NA_real_,resid_sd=NA_real_,max_leverage=NA_real_,df_resid=NA_real_)
       if(stat=='PASS') {
          fit <- tryCatch(fit_hc3(Y[,ct],X,which(colnames(X)=='case')),
                  error=function(e) e)
          if(inherits(fit,'error'))stat <- paste0('FAIL_FIT:',conditionMessage(fit)) else est[names(fit)] <- fit
       }
       rows[[length(rows)+1L]] <- data.table(dataset=ds,disease=DISEASE[[ds]],method=mt,
         cell_type=ct,scale=transform,n_total=nrow(md),n_used=nrow(kept),
         n_control=sum(kept$case==0),n_case=sum(kept$case==1),
         effect=est['effect'],ci_low=est['ci_low'],ci_high=est['ci_high'],
         robust_se=est['se'],p_hc3=est['p'],df_resid=est['df_resid'],
         mean_zero_fraction=zero[[ct]],sd_pp=sds[[ct]],
         epsilon=eps,status=stat,covariates=paste(covars,collapse=';'),
         interpretation=if(transform=='percentage')'Adjusted pp, exploratory relative signal; OLS HC3' else
                      'Centered log ratio; NOT percentage points; pseudocount sensitive')
     }
   }
 }
 logstat(ds,'PASS_ADJUSTED_EXPLORATORY',paste('n=',nrow(kept),' p=',ncol(X),
      '; OLS(HC3) + CLR sensitivities. Not biological validation.',sep=''))
}
if(length(samples))WR(rbindlist(samples,fill=TRUE),'11_included_excluded_samples.tsv')
if(length(qc))WR(rbindlist(qc,fill=TRUE),'11_design_audit.tsv')
if(length(cov_bal))WR(rbindlist(cov_bal,fill=TRUE),'11_covariate_balance.tsv')
if(length(rows)) {
 res <- rbindlist(rows,fill=TRUE)
 res[,fdr_bh_8types:=p.adjust(p_hc3,method='BH'),by=.(dataset,method,scale)]
 WR(res,'11_adjusted_contrasts_all.tsv')
 WR(res[scale=='percentage'],'11_adjusted_pp_HC3.tsv')
 WR(res[scale!='percentage'],'11_clr_sensitivity.tsv')
 comp <- dcast(res[scale!='percentage',.(dataset,method,cell_type,scale,effect,status)],
      dataset+method+cell_type~scale,value.var='effect')
 if(all(c('clr_1e-4','clr_1e-5')%in%names(comp)))
   comp[,direction_same := sign(`clr_1e-4`)==sign(`clr_1e-5`)]
 WR(comp,'11_clr_epsilon_stability.tsv')
 plotdt <- res[scale=='percentage' & status=='PASS' & is.finite(effect)]
 if(nrow(plotdt)) {
   g <- ggplot(plotdt,aes(x=effect,y=cell_type,color=method))+
      geom_vline(xintercept=0,linetype='dashed',color='grey60')+
      geom_segment(aes(x=ci_low,xend=ci_high,yend=cell_type),position=position_dodge(width=.6),alpha=.75)+
      geom_point(position=position_dodge(width=.6),size=2)+
      facet_wrap(~dataset,scales='free_x',ncol=2)+
      labs(title='11 | Diferencas caso-controle ajustadas (OLS HC3)',
           subtitle='Estimativas relativas; 95% CI robusto; apenas termos aprovados; NA nao plotados',
           x='Efeito caso vs controle (pontos percentuais)',y='Tipo celular')+
      theme_minimal(base_size=11)+theme(legend.position='bottom',strip.text=element_text(face='bold'))
   ggsave(file.path(PLOTS,'11_adjusted_effects_HC3.png'),g,width=13,height=9,dpi=170,bg='white')
 }
}
WR(rbindlist(status,fill=TRUE),'11_status.tsv')
writeLines(capture.output(sessionInfo()),file.path(OUT,'11_sessionInfo.txt'))
cat('11 terminado. Confira STATUS; PASS do modelo nao significa conclusao biologica. Saidas: ',OUT,'\n')
