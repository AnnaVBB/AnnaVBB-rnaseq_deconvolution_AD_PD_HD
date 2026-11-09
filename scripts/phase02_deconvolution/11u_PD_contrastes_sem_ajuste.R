#!/usr/bin/env Rscript
# ==============================================================================
# 11u | PARKINSON - CONTRASTES EXPLORATORIOS SEM AJUSTE (NAO INFERENCIA AJUSTADA)
# Apenas quando idade/sexo/RIN/PMI nao estao observados em ambos os grupos.
# dtangle 08e pooled_counts e o metodo principal, conforme fluxo 09.
# Preserva as mesmas regras de aviso para tipos colapsados utilizadas no script 11.
# Executar: Rscript scripts/phase02_deconvolution/11u_PD_contrastes_sem_ajuste.R
# ==============================================================================
options(stringsAsFactors=FALSE)
if (!requireNamespace('data.table',quietly=TRUE)) stop('Pacote data.table ausente')
suppressPackageStartupMessages(library(data.table))
ROOT <- normalizePath(Sys.getenv('DECONV_ROOT','.'),mustWork=TRUE)
P <- function(...) file.path(ROOT,...)
IN <- P('results','deconvolution','Integrated_Comparison','09','09_proportions_long_harmonized.tsv')
META <- P('data','processed','deconvolution','GSE68719','metadata_deconvolution.tsv')
NOTES <- P('results','deconvolution','Integrated_Comparison','09','09_LEIA_ANTES_DE_INTERPRETAR.txt')
OUT <- P('results','deconvolution','Integrated_Comparison','11_unadjusted_PD')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
WR <- function(x,name) fwrite(x,file.path(OUT,name),sep='\t',na='NA',quote=FALSE)
if(!file.exists(IN) || !file.exists(META))stop('Entradas faltantes: etapa 09 ou metadata etapa 05')
if(file.exists(NOTES)){
  note_text <- paste(readLines(NOTES,warn=FALSE),collapse=' ')
  if (!grepl('08e dtangle/pooled_counts',note_text,fixed=TRUE))
    stop('Entrada 09 NAO confirmou dtangle/08e pooled_counts. Conferir configuracao antes de executar.')
  prov <- 'CONFIRMED_09_NOTES_08e_pooled_counts'
} else {
  prov <- 'UNVERIFIED_NOTES_FILE_ABSENT'
  warning('09_LEIA_ANTES_DE_INTERPRETAR.txt ausente. Prova de origem do dtangle deve ser conferida ANTES da interpretacao.')
}
D <- fread(IN,na.strings=c('NA','','NaN'))
MD <- fread(META,na.strings=c('NA','','NaN'))
needed <- c('dataset','method','cell_type','sample_id','group','proportion')
if(!all(needed%in%names(D)))stop('Entrada 09 incompleta')
if(!all(c('sample_id','group')%in%names(MD)))stop('Metadata PD incompleto')
if(anyDuplicated(MD,by='sample_id'))stop('PD metadata IDs duplicados')
MD <- MD[,.(sample_id=as.character(sample_id),group=as.character(group))]
MD <- MD[order(sample_id)]
if(nrow(MD)!=73L || sum(MD$group=='PD')!=29L || sum(MD$group=='Control')!=44L)
  stop('PD: grupos/n esperados (29 PD, 44 Control) nao coincidem com metadata')
D <- D[dataset=='GSE68719']
if(anyDuplicated(D,by=c('method','sample_id','cell_type')))stop('09 PD possui duplicatas')
CT <- c('Ex','In','Ast','Oli','Opc','Mic','End','Per')
MT <- c('dtangle','BisqueRNA','MuSiC2')
if(!setequal(unique(D$method),MT)||!setequal(unique(D$cell_type),CT))stop('Metodos/tipos inesperados no PD')
chk <- D[,.(n=.N,tot=sum(proportion)),by=.(method,sample_id)]
if(any(chk$n!=8L)|any(abs(chk$tot-1)>1e-6))stop('Fractions PD incompletas/soma != 1')
if(any(!is.finite(D$proportion))||any(D$proportion<0|D$proportion>1))stop('Fractions fora de [0,1]')
# OLS de proporcao em pontos percentuais com HC3 (apenas group/case), sem covariaveis.
hc3 <- function(y, X) {
  n <- length(y);p <- ncol(X)
  if(nrow(X)!=n||qr(X)$rank!=p)stop('Matriz de desenho invalida')
  Q <- chol2inv(chol(crossprod(X)))
  beta <- as.numeric(Q%*%crossprod(X,y)); e <- as.numeric(y-X%*%beta)
  h <- rowSums((X%*%Q)*X)
  if(any(h>=1-1e-8))stop('HC3 indefinido')
  V <- Q %*% crossprod(X, X*as.numeric((e/(1-h))^2)) %*% Q
  b <- unname(beta[2]); s <- unname(sqrt(V[2,2])); df <- n-p
  if(!is.finite(s)||s<=0)stop('Erro padrao HC3 invalido')
  q <- qt(0.975,df)
  c(effect=b,ci_low=b-q*s,ci_high=b+q*s,robust_se=s,
    p_hc3=2*pt(-abs(b/s),df),df_resid=df)
}
rows <- list()
for(mt_name in MT){
  wide <- dcast(D[method==mt_name],sample_id+group~cell_type,value.var='proportion')
  if(!setequal(wide$sample_id,MD$sample_id))stop(mt_name,': amostras incompletas')
  wide <- wide[match(MD$sample_id,wide$sample_id)]
  if(!identical(as.character(wide$group),MD$group))stop(mt_name,': grupo conflitante')
  M <- as.matrix(wide[,..CT]); storage.mode(M)<-'double'
  if(anyNA(M)||max(abs(rowSums(M)-1))>1e-6)stop('Fractions inconsistentes')
  X <- cbind('(Intercept)'=1,case=as.integer(MD$group=='PD'))
  sets <- c(as.list(setNames(CT,CT)),list('Ex+In'=c('Ex','In'),
          'Oli+Opc'=c('Oli','Opc'),'End+Per'=c('End','Per')))
  for(nm in names(sets)){
    members <- sets[[nm]]
    frac <- rowSums(M[,members,drop=FALSE])
    zero_fraction <- mean(frac<=1e-12)
    sd_pp <- sd(100*frac)
    collapsed_members <- CT[colMeans(M<=1e-12)>=0.90 & CT%in%members]
    status <- if(zero_fraction>=0.90 || sd_pp<0.20) 'SKIP_COLLAPSED_TYPE' else 'PASS_UNADJUSTED'
    est <- c(effect=NA_real_,ci_low=NA_real_,ci_high=NA_real_,robust_se=NA_real_,p_hc3=NA_real_,df_resid=NA_real_)
    if(status=='PASS_UNADJUSTED'){
      tmp <- tryCatch(hc3(100*frac,X),error=function(e)e)
      if(inherits(tmp,'error'))status<-paste0('FAIL_FIT:',conditionMessage(tmp)) else est<-tmp
    }
    rows[[length(rows)+1L]] <- data.table(dataset='GSE68719',disease='PD',
      method=mt_name,cell_type=nm,analysis_family=if(length(members)==1L)'single_8types' else 'aggregate_3groups',
      analysis_class='UNADJUSTED',n_used=nrow(wide),n_control=sum(MD$group=='Control'),n_case=sum(MD$group=='PD'),
      effect=unname(est['effect']),ci_low=unname(est['ci_low']),ci_high=unname(est['ci_high']),
      robust_se=unname(est['robust_se']),p_hc3=unname(est['p_hc3']),df_resid=unname(est['df_resid']),
      mean_zero_fraction=zero_fraction,sd_pp=sd_pp,
      collapsed_constituents=paste(collapsed_members,collapse=';'),
      status=status,covariates='',reference_provenance=prov,
      note='Contraste bruto caso-controle; NAO controla idade, sexo, RIN, PMI ou batch')
  }
}
res <- rbindlist(rows,fill=TRUE)
res[,fdr_bh_family:=p.adjust(p_hc3,'BH'),by=.(method,analysis_family)]
if(nrow(res[status=='PASS_UNADJUSTED' & (!is.finite(effect)|!is.finite(p_hc3))]))
  stop('Valores faltantes em modelos marcados PASS_UNADJUSTED')
WR(res,'11u_PD_unadjusted_effects.tsv')
WR(data.table(dataset='GSE68719',status='EXPLORATORY_UNADJUSTED_ONLY',
  n_case=29L,n_control=44L,shared_controls_with_HD=44L,reference_provenance=prov,
  warning='NAO combinar p-valores como independentes com HD e NAO rotular como ajustado'),
  '11u_PD_status.tsv')
writeLines(c('PD sem ajuste: OLS(HC3) aplicado a diferenca percentual caso-controle.',
  'FDR BH por metodo e familia (oito tipos / tres agrupamentos).',
  'NAO constitui inferencia ajustada. Confundimento por idade/sexo/RIN/PMI pode persistir.',
  'HD e PD compartilham 44 controles; estimativas nao sao independentes.',
  'Manter regras de cautela para MuSiC2 (In/Opc) e tipos de baixa abundancia.'),
  file.path(OUT,'11u_INTERPRETACAO.txt'))
cat('11u concluido; contraste PD SEM AJUSTE em ',OUT,'\n',sep='')
