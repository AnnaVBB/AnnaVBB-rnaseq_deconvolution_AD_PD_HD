#!/usr/bin/env Rscript
# =============================================================================
# 11b | Sensibilidade por agrupamentos celulares - AD, HD, PD
# Usa exatamente os metadados/covariaveis aprovados na etapa 11 (v3)
# e as MESMAS fracoes da matriz harmonizada da etapa 09 (dtangle 08e).
# Nao soma p-valores nem erros-padrao. Reajusta os modelos com fracoes somadas.
# OLS (%) com SE HC3; CLR dos CINCO compartimentos como sensibilidade.
# =============================================================================
options(stringsAsFactors = FALSE)
if (!requireNamespace('data.table', quietly=TRUE)) stop('Pacote data.table necessario')
suppressPackageStartupMessages(library(data.table))
ROOT <- normalizePath(Sys.getenv('DECONV_ROOT','.'), mustWork=TRUE)
P <- function(...) file.path(ROOT,...)
INPUT <- P('results','deconvolution','Integrated_Comparison','09','09_proportions_long_harmonized.tsv')
PLAN <- P('data','processed','deconvolution','10_approved_model_plan.tsv')
META <- P('results','deconvolution','Integrated_Comparison','10_metadata')
OUT <- P('results','deconvolution','Integrated_Comparison','11b_grouped')
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
WR <- function(x,filename) fwrite(x,file.path(OUT,filename),sep='\t',quote=FALSE,na='NA')
DS <- c('GSE125583','ROSMAP','GSE64810','GSE68719')
DISEASE <- c(GSE125583='AD',ROSMAP='AD',GSE64810='HD',GSE68719='PD')
METHODS <- c('dtangle','BisqueRNA','MuSiC2')
CELLS <- c('Ex','In','Ast','Oli','Opc','Mic','End','Per')
GROUPS <- list(Neuronios='Ex+In',Oligodendroglia='Oli+Opc',Vascular='End+Per')
GROUP_MEMBERS <- list(Neuronios=c('Ex','In'),Oligodendroglia=c('Oli','Opc'),Vascular=c('End','Per'))
PARTS <- c('Neuronios','Oligodendroglia','Vascular','Ast','Mic')
fields <- function(s) { if(length(s)!=1L||is.na(s)||!nzchar(trimws(s))) return(character()); trimws(strsplit(s,',',fixed=TRUE)[[1]]) }
nonblank <- function(z) !is.na(z)&nzchar(trimws(as.character(z)))
fit_hc3 <- function(y,X) {
 n<-nrow(X); p<-ncol(X)
 if(n<=p+5L||qr(X)$rank!=p) stop('Matriz singular ou amostra insuficiente')
 Q<-chol2inv(chol(crossprod(X)))
 beta<-as.numeric(Q%*%crossprod(X,y)); names(beta)<-colnames(X)
 rr<-as.numeric(y-X%*%beta)
 h<-rowSums((X%*%Q)*X)
 if(any(h>=1-1e-8)) stop('HC3: alavancagem unidade')
 v<-Q%*%crossprod(X,X*as.numeric((rr/(1-h))^2))%*%Q
 ss<-sqrt(diag(v)); names(ss)<-colnames(X)
 b<-unname(beta['case']); se<-unname(ss['case'])
 if(!is.finite(b)||!is.finite(se)||se<=0) stop('Coeficiente/SE nao estimavel')
 df<-n-p; tcrit<-qt(.975,df)
 c(effect=b,ci_low=b-tcrit*se,ci_high=b+tcrit*se,robust_se=se,
   p_hc3=2*pt(-abs(b/se),df),df_resid=df,max_leverage=max(h))
}
if(!file.exists(INPUT)||!file.exists(PLAN))stop('Execute etapas 09/10/11; entrada 09 ou plano ausente')
D<-fread(INPUT,na.strings=c('NA','','NaN'))
need<-c('dataset','method','cell_type','sample_id','group','proportion')
if(!all(need%in%names(D)))stop('Faltam colunas 09')
if(anyDuplicated(D,by=c('dataset','method','cell_type','sample_id')))stop('Chaves 09 duplicadas')
if(anyNA(D$proportion)||any(!is.finite(D$proportion))||any(D$proportion<0|D$proportion>1))stop('Fracoes invalidas')
chk<-D[,.(n=.N,total=sum(proportion)),by=.(dataset,method,sample_id)]
if(any(chk$n!=8L)||any(abs(chk$total-1)>1e-6))stop('Matriz 09 incompleta')
plan<-fread(PLAN,na.strings=c('NA','','NaN'))
required_plan<-c('dataset','approved','numeric_covariates','categorical_covariates','cluster_column','justification','provenance')
if(!all(required_plan%in%names(plan)))stop('Plano incompleto')
if(nrow(plan)!=4L||!setequal(plan$dataset,DS)||anyDuplicated(plan,by='dataset'))stop('Plano das 4 coortes invalido')
results<-list(); checks<-list()
for(ds in DS) {
 pp<-plan[dataset==ds]
 if(toupper(trimws(pp$approved))!='YES')stop('Coorte nao aprovada: ',ds)
 if(nonblank(pp$cluster_column))stop('Plano exige modelo clusterizado, HC3 comum nao basta: ',ds)
 if(!nonblank(pp$justification)||!nonblank(pp$provenance))stop('Justificativa/procedencia ausente: ',ds)
 nums<-fields(pp$numeric_covariates); cats<-fields(pp$categorical_covariates); cv<-c(nums,cats)
 if(!length(cv)||anyDuplicated(cv))stop('Covariaveis invalidas: ',ds)
 path<-file.path(META,paste0('10_metadata_',ds,'_audited.tsv'))
 if(!file.exists(path))stop('Execute o 10: ',path)
 md<-fread(path,na.strings=c('NA','','NULL','NaN'))
 if(!all(c('sample_id','group',cv)%in%names(md))||anyDuplicated(md,by='sample_id'))stop('Metadados invalidos: ',ds)
 md<-md[,c('sample_id','group',cv),with=FALSE]
 for(cc in nums) if(!is.numeric(md[[cc]])) stop(ds,': numerica nao numerica ',cc)
 for(cc in cats) md[[cc]]<-factor(trimws(as.character(md[[cc]])))
 keep<-complete.cases(md)
 for(cc in cv)keep<-keep&nonblank(md[[cc]])
 kept<-copy(md[keep]); kept[,case:=as.integer(group==DISEASE[[ds]])]
 if(nrow(kept)<30L||sum(kept$case==0L)<10L||sum(kept$case==1L)<10L)stop('Sem amostra suficiente: ',ds)
 for(cc in cats)kept[[cc]]<-droplevels(kept[[cc]])
 if(any(vapply(cats,function(cc)nlevels(kept[[cc]])<2L,logical(1))))stop('Categoria constante: ',ds)
 X<-model.matrix(reformulate(c('case',cv)),data=as.data.frame(kept))
 if(!'case'%in%colnames(X)||qr(X)$rank!=ncol(X)||nrow(X)<max(30L,10L*ncol(X)))stop('Desenho inadequado: ',ds)
 for(mt in METHODS) {
  long<-D[dataset==ds&method==mt]
  wide<-dcast(long,sample_id~cell_type,value.var='proportion')
  if(!setequal(wide$sample_id,md$sample_id))stop('Mapeamento divergente: ',ds,'/',mt)
  wide<-wide[match(kept$sample_id,wide$sample_id)]
  if(!identical(wide$sample_id,kept$sample_id))stop('Ordem de amostras corrompida')
  M<-as.matrix(wide[,..CELLS]);storage.mode(M)<-'double'
  if(anyNA(M)||any(abs(rowSums(M)-1)>1e-6))stop('Frações inválidas em ',ds,'/',mt)
  G<-sapply(GROUP_MEMBERS,function(z) rowSums(M[,z,drop=FALSE]))
  F<-cbind(G,Ast=M[,'Ast'],Mic=M[,'Mic'])
  if(any(abs(rowSums(F)-1)>1e-6))stop('Cinco compartimentos nao somam 1')
  for(scale in c('percentage','clr_1e-4','clr_1e-5')) {
   eps<-if(scale=='percentage')NA_real_ else if(scale=='clr_1e-4')1e-4 else 1e-5
   Y<-if(scale=='percentage')100*F else {Z<-log(F+eps);Z-rowMeans(Z)}
   any_collapsed<-any(colMeans(F<=1e-12)>=0.90)
   for(g in names(GROUP_MEMBERS)) {
    zero<-mean(F[,g]<=1e-12); sdpp<-sd(F[,g])*100
    state<-'PASS'
    if(scale=='percentage'&&(zero>=.90||sdpp<.20))state<-'SKIP_COLLAPSED_GROUP'
    if(scale!='percentage'&&any_collapsed)state<-'SKIP_CLR_STRUCTURAL_ZERO'
    est<-rep(NA_real_,7L);names(est)<-c('effect','ci_low','ci_high','robust_se','p_hc3','df_resid','max_leverage')
    if(state=='PASS'){
      z<-tryCatch(fit_hc3(Y[,g],X),error=function(e)e)
      if(inherits(z,'error'))state<-paste0('FAIL_FIT:',conditionMessage(z)) else est<-z
    }
    results[[length(results)+1L]]<-data.table(dataset=ds,disease=DISEASE[[ds]],method=mt,
     group_type=g,members=paste(GROUP_MEMBERS[[g]],collapse='+'),scale=scale,
     n_total=nrow(md),n_used=nrow(kept),n_control=sum(kept$case==0),n_case=sum(kept$case==1),
     effect=unname(est['effect']),ci_low=unname(est['ci_low']),ci_high=unname(est['ci_high']),
     robust_se=unname(est['robust_se']),p_hc3=unname(est['p_hc3']),
     df_resid=unname(est['df_resid']),max_leverage=unname(est['max_leverage']),
     zero_fraction=zero,sd_pp=sdpp,epsilon=eps,status=state,
     covariates=paste(cv,collapse=';'),
     caution=if(g=='Vascular')'End+Per: low abundance and marker sensitivity' else if(mt=='MuSiC2')
       'MuSiC2 In/Opc collapsed individually; aggregate does not validate components' else '')
   }
  }
  checks[[length(checks)+1L]]<-data.table(dataset=ds,method=mt,n=nrow(F),n_parts=ncol(F),
      max_sum_error=max(abs(rowSums(F)-1)),n_parameters=ncol(X),status='PASS_INPUT')
 }
}
r<-rbindlist(results,fill=TRUE)
r[,fdr_bh_3groups:=p.adjust(p_hc3,method='BH'),by=.(dataset,method,scale)]
if(any(r[status=='PASS',!is.finite(effect)|!is.finite(p_hc3)]))stop('PASS com coeficiente/p NA')
WR(r,'11b_grouped_all.tsv')
WR(r[scale=='percentage'],'11b_grouped_adjusted_pp_HC3.tsv')
WR(r[scale!='percentage'],'11b_grouped_clr_sensitivity.tsv')
WR(rbindlist(checks),'11b_input_checks.tsv')
writeLines(c('Agrupamentos: Ex+In, Oli+Opc, End+Per; modelos reestimados, nunca soma de p/SE.',
 'CLR refeito na composicao de CINCO partes: Neuronios, Oligodendroglia, Vascular, Ast, Mic.',
 'FDR BH por coorte/metodo/escala entre os tres contrastes agrupados.',
 'Fracoes harmonizadas 09: verificar documentacao que dtangle e 08e pooled_counts.',
 'HD e PD possuem 44 controles compartilhados: sem teste conjunto ingenuo.',
 'ROSMAP validacao AD; GSE125583 principal. Sem inferencia conjunta automatica.'),file.path(OUT,'11b_LEIA_ANTES.txt'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'11b_sessionInfo.txt'))
cat('11b CONCLUIDO: ',OUT,'\n',sep='')
