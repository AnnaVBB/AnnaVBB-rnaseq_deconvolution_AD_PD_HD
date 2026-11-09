#!/usr/bin/env Rscript
# ==============================================================================
# 10e | PREPARAR E AUDITAR EXTENSOES HD, PD E ROSMAP PARA MODELOS 11
# TCC deconvolucao AD / HD / PD | nao sobrescreve entradas nem aprova modelos
# Requer: 10c e 10d concluidos; executado na raiz do repositorio.
# ==============================================================================
options(stringsAsFactors = FALSE)
if (!requireNamespace('data.table', quietly=TRUE)) stop('Pacote data.table ausente')
suppressPackageStartupMessages(library(data.table))
ROOT <- normalizePath(Sys.getenv('DECONV_ROOT','.'), mustWork=TRUE)
P <- function(...) file.path(ROOT,...)
BASE <- P('data','processed','deconvolution')
D10 <- P('results','deconvolution','Integrated_Comparison','10_metadata')
OUT <- file.path(D10,'10e_preparacao_modelos')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
WR <- function(x,nm) fwrite(x,file.path(OUT,nm),sep='\t',na='NA',quote=FALSE)
READ <- function(path,cols,label) {
  if (!file.exists(path)) stop(label, ': ausente: ',path)
  d <- fread(path,na.strings=c('','NA','NaN','NULL'))
  if (!all(cols %in% names(d))) stop(label, ': colunas ausentes: ',paste(setdiff(cols,names(d)),collapse=', '))
  if (!('sample_id' %in% names(d))) stop(label,': faltou sample_id')
  d[,sample_id:=trimws(as.character(sample_id))]
  if (anyNA(d$sample_id)||any(d$sample_id=='')||anyDuplicated(d,by='sample_id')) stop(label, ': IDs invalidos')
  d
}
meta <- function(ds, n, caso) {
  d <- READ(file.path(BASE,ds,'metadata_deconvolution.tsv'),c('sample_id','group'),paste0('metadata ',ds))
  if(nrow(d)!=n||!setequal(unique(d$group),c('Control',caso)))stop(ds, ': n ou grupos inesperados')
  d[,.(sample_id,group)]
}
num <- function(x,label) {
  s <- trimws(as.character(x)); na <- is.na(x)|s==''
  ok <- na|grepl('^[0-9]+(\\.[0-9]+)?$',s)
  if(any(!ok))stop(label,': valores nao numericos: ',paste(head(unique(s[!ok]),5),collapse=';'))
  y <- suppressWarnings(as.numeric(s)); y[na] <- NA_real_
  if(any(!is.finite(y[!is.na(y)]))||any(y[!is.na(y)]<0))stop(label, ': valores invalidos')
  y
}
prepare <- function(ds,n,caso,source,numeric_cols,categorical_cols=character()) {
  m <- meta(ds,n,caso)
  src <- READ(source,c('sample_id',numeric_cols,categorical_cols),paste0('candidato ',ds))
  if(!setequal(m$sample_id,src$sample_id))stop(ds, ': IDs do candidato diferem do metadata base')
  src <- src[match(m$sample_id,src$sample_id)]
  stopifnot(identical(src$sample_id,m$sample_id))
  z <- data.table(sample_id=src$sample_id)
  for(v in numeric_cols) z[,(v):=num(src[[v]],paste(ds,v))]
  for(v in categorical_cols) {
    s <- trimws(as.character(src[[v]])); s[is.na(src[[v]])|s=='']<-NA_character_
    z[,(v):=s]
  }
  if(ds=='ROSMAP') {
    if(any(!is.na(z$msex)&!z$msex %in% c('0','1'))) stop('ROSMAP: msex deve ser codificacao literal 0/1')
    if(anyNA(z$sequencingBatch)) stop('ROSMAP: batch ausente')
  }
  allvars <- c(numeric_cols,categorical_cols)
  merged <- merge(m,z,by='sample_id',all.x=TRUE,sort=FALSE)
  merged <- merged[match(m$sample_id,merged$sample_id)]
  stopifnot(identical(merged$sample_id,m$sample_id))
  complete <- complete.cases(merged[, ..allvars])
  detail <- merged[, .(n=.N,n_complete=sum(complete.cases(.SD))),by=group,.SDcols=allvars]
  detail[,`:=`(dataset=ds,model=paste('case+',paste(allvars,collapse='+')))]
  if(any(detail$n_complete<10L))stop(ds,': um grupo tem menos de 10 amostras completas')
  df <- as.data.frame(merged[complete])
  df$case <- as.integer(df$group==caso)
  for(v in categorical_cols)df[[v]]<-factor(df[[v]])
  fm <- reformulate(c('case',allvars))
  X <- model.matrix(fm,df)
  rank <- qr(X)$rank
  # Mantem a mesma verificacao de n/p do script 11.
  if(rank!=ncol(X))stop(ds,': desenho singular; revisar covariaveis e lotes')
  if(nrow(X)<max(30,10*ncol(X)))stop(ds,': amostra insuficiente para numero de parametros')
  # Leverage e condicionamento sao diagnosticos, nao prova de boa inferencia.
  leverage <- rowSums(qr.Q(qr(X))^2)
  if(any(leverage>=1-1e-8))stop(ds,': observacao com alavancagem 1; HC3 indefinido')
  des <- data.table(dataset=ds,n_total=nrow(m),n_complete=sum(complete),
                    n_control=sum(df$case==0),n_case=sum(df$case==1),
                    n_parameters=ncol(X),design_rank=rank,
                    condition_number=kappa(X),max_leverage=max(leverage),
                    formula=paste(deparse(fm),collapse=''),status='PASS_DESIGN_CANDIDATE')
  WR(z,paste0('10e_',ds,'_covariates_READY_FOR_REVIEW.tsv'))
  WR(merged[,c('sample_id','group',allvars),with=FALSE],paste0('10e_',ds,'_preview.tsv'))
  list(coverage=detail,design=des,ext=z)
}
D10c <- file.path(D10,'10c_HD_PD')
D10d <- file.path(D10,'10d_preparacao')
hd_path <- file.path(D10c,'10c_GSE64810_covariates_REVIEW_FIRST.tsv')
pd_path <- file.path(D10d,'10d_GSE68719_covariates_REVIEW_FIRST.tsv')
ro_path <- file.path(D10d,'10d_ROSMAP_covariates_REVIEW_FIRST.tsv')
# Valida integralmente o crosswalk previamente confirmado.
cross <- READ(file.path(D10c,'10c_PD_HD_shared_controls_crosswalk.tsv'),
              c('sample_id','status_match','age_conflict','sex_conflict'),'crosswalk HD/PD')
pd_meta <- meta('GSE68719',73L,'PD'); hd_meta <- meta('GSE64810',69L,'HD')
if(nrow(cross)!=44L || !setequal(cross$sample_id,pd_meta[group=='Control',sample_id]) ||
   !all(cross$sample_id %in% hd_meta[group=='Control',sample_id]) ||
   any(cross$status_match!='EXACT_CONTROL_IN_BOTH') ||
   any(cross$age_conflict %in% TRUE) || any(cross$sex_conflict %in% TRUE))
  stop('Cruzamento dos 44 controles NAO passou na auditoria')
# Confere que o 10d preservou exatamente a idade dos controles do 10c/HD.
hd_all <- READ(hd_path,c('sample_id','age','rin','pmi'),'10c HD')
pd_all <- READ(pd_path,c('sample_id','age'),'10d PD')
ids <- pd_meta[group=='Control',sample_id]
a_hd <- num(hd_all$age[match(ids,hd_all$sample_id)],'controle HD age')
a_pd <- num(pd_all$age[match(ids,pd_all$sample_id)],'controle PD age')
if(anyNA(a_hd)||anyNA(a_pd)||any(a_hd!=a_pd)) stop('PD: idade dos 44 controles nao coincide com HD')
cat('Crosswalk: 44/44 controles com idade transferida por ID exato.\n')
A <- list(
  prepare('GSE64810',69L,'HD',hd_path,c('age','rin','pmi')),
  prepare('GSE68719',73L,'PD',pd_path,c('age')),
  prepare('ROSMAP',144L,'AD',ro_path,c('rin','pmi'),c('msex','sequencingBatch'))
)
WR(rbindlist(lapply(A,`[[`,'coverage')),'10e_coverage.tsv')
WR(rbindlist(lapply(A,`[[`,'design')),'10e_design_check.tsv')
# Proposta: mantem a decisao sobre GSE125583 e NAO altera o plano definitivo.
plan_path <- file.path(BASE,'10_approved_model_plan.tsv')
if(file.exists(plan_path)){
 p <- fread(plan_path,na.strings=c('NA','NaN'))
 if(nrow(p)!=4L||anyDuplicated(p,by='dataset'))stop('Plano atual tem coortes invalidas')
 p[dataset=='GSE64810',`:=`(approved='NO',numeric_covariates='age,rin,pmi',
       categorical_covariates='',
       justification='Proposta: idade ao obito, integridade de RNA e intervalo post mortem; n=68 completo; revisar fonte/unidade',
       provenance='GEO GSE64810 SOFT via 10b e 10c, controle compartilhado com PD')]
 p[dataset=='GSE68719',`:=`(approved='NO',numeric_covariates='age',
       categorical_covariates='',
       justification='Proposta: idade ao obito presente em ambos grupos; nenhum ajuste possivel por sexo/RIN/PMI',
       provenance='PD idade ao obito GSE68719 GEO; controles compartilhados do GSE64810; 10c/10d')]
 p[dataset=='ROSMAP',`:=`(approved='NO',numeric_covariates='rin,pmi',
       categorical_covariates='msex,sequencingBatch',
       justification='Proposta: qualidade RNA, PMI, sexo e lote; NYGC2 pequeno; idade indisponivel nesta fonte',
       provenance='ROSMAP metadata da fase 01c, IDs correspondidos pelo 10d')]
 WR(p,'10e_model_plan_PROPOSTA_NAO_APROVADA.tsv')
}
writeLines(c(
 '10e nao altera data/processed nem executa regressao.',
 'Revise as unidades e definicoes de age_of_death, age_at_death, RIN e PMI nas fontes originais.',
 'Verifique a tabela de lotes ROSMAP: NYGC2 tem apenas quatro amostras (3 AD, 1 controle).',
 'Aprovacao exige revisar 10e_design_check e 10e_coverage e conferir a origem do dtangle 08e na etapa 09.',
 'Depois copie 10e_<DATASET>_covariates_READY_FOR_REVIEW.tsv para data/processed/deconvolution/metadata_extension/<DATASET>_covariates.tsv.',
 'Reexecute script 10, edite manualmente o plano aprovado (YES somente para modelos aceitos), e depois execute o script 11 v3.',
 'PD e HD compartilham 44 controles. Nao combine seus efeitos como independentes na integracao.',
 'Resultados estatisticos a partir de dados bulk sao efeitos relativos, nao contagens absolutas de celulas.'
),file.path(OUT,'10e_LEIA_ANTES.txt'))
cat('\n10e terminou sem alterar os dados ou plano oficiais. Resultados: ',OUT,'\n',sep='')
