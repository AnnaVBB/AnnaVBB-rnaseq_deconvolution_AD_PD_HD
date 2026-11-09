#!/usr/bin/env Rscript
# =============================================================================
# 13b | INTEGRACAO DESCRITIVA DE EFEITOS CELULARES AJUSTADOS
# AD principal GSE125583; AD validacao ROSMAP; HD GSE64810; PD GSE68719.
# NUNCA meta-analisa ou combina valores-p: HD/PD compartilham 44 controles.
# 11b opcional (se rodado, incorpora compartimentos com seus IC HC3 validos).
# =============================================================================
options(stringsAsFactors=FALSE)
if(!requireNamespace('data.table',quietly=TRUE))stop('Precisa data.table')
suppressPackageStartupMessages(library(data.table))
ROOT<-normalizePath(Sys.getenv('DECONV_ROOT','.'),mustWork=TRUE)
P<-function(...)file.path(ROOT,...)
INPUT<-P('results','deconvolution','Integrated_Comparison','11_adjusted','11_adjusted_pp_HC3.tsv')
GROUPED<-P('results','deconvolution','Integrated_Comparison','11b_grouped','11b_grouped_adjusted_pp_HC3.tsv')
OUT<-P('results','deconvolution','Integrated_Comparison','13b_integracao_AD_main_ROSMAP_validacao')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
WR<-function(x,nm)fwrite(x,file.path(OUT,nm),sep='\t',quote=FALSE,na='NA')
if(!file.exists(INPUT))stop('Falta 11 ajustado: ',INPUT)
E<-fread(INPUT,na.strings=c('NA','','NaN'))
need<-c('dataset','method','cell_type','scale','status','effect','ci_low','ci_high','p_hc3','fdr_bh_8types','covariates')
if(!all(need%in%names(E)))stop('Entrada incompleta')
if(!setequal(unique(E$dataset),c('GSE125583','ROSMAP','GSE64810','GSE68719')))stop('Faltam coortes do 11')
if(anyDuplicated(E,by=c('dataset','method','cell_type','scale')))stop('Chaves 11 duplicadas')
E<-E[scale=='percentage' & status=='PASS',.(dataset,method,cell_type,
  effect,ci_low,ci_high,p_hc3,fdr=fdr_bh_8types,covariates)]
if(any(!is.finite(E$effect)|!is.finite(E$ci_low)|!is.finite(E$ci_high)|!is.finite(E$p_hc3)))stop('PASS invalido')
E[,`:=`(AD_role=fifelse(dataset=='GSE125583','AD_MAIN',fifelse(dataset=='ROSMAP','AD_VALIDATION','OTHER')),
 direction=ifelse(effect<0,'DOWN',ifelse(effect>0,'UP','ZERO')),
 fdr_significant=(!is.na(fdr)&fdr<.05),
 rare_type_caution=cell_type%in%c('Mic','End','Per'))]
WR(E,'13b_individual_effects.tsv')
CT<-c('Ex','In','Ast','Oli','Opc','Mic','End','Per'); METHODS<-c('dtangle','BisqueRNA','MuSiC2')
make_scenario<-function(ad) {
  rows<-list()
  for(mt in METHODS)for(ct in CT){
    opt<-c(AD=ad,HD='GSE64810',PD='GSE68719')
    vals<-lapply(opt,function(ds)E[dataset==ds & method==mt & cell_type==ct])
    complete<-all(vapply(vals,nrow,integer(1))==1L)
    effects<-vapply(vals,function(z)if(nrow(z)==1L)z$effect[1L] else NA_real_,numeric(1))
    sig<-vapply(vals,function(z)if(nrow(z)==1L)z$fdr_significant[1L] else NA,logical(1))
    direction_same<-if(complete) length(unique(sign(effects)))==1L else NA
    rows[[length(rows)+1L]]<-data.table(AD_dataset=ad,method=mt,cell_type=ct,
      AD_effect_pp=effects['AD'],HD_effect_pp=effects['HD'],PD_effect_pp=effects['PD'],
      AD_fdr05=sig['AD'],HD_fdr05=sig['HD'],PD_fdr05=sig['PD'],
      n_diseases_available=sum(is.finite(effects)),same_direction=direction_same,
      status=if(!complete)'INCOMPLETE' else if(direction_same)'SAME_DIRECTION_DESCRIPTIVE' else 'DIVERGENT_DIRECTION',
      warning=if(ct%in%c('Mic','End','Per'))'Low abundance; exploratory' else if(mt=='MuSiC2'&&ct%in%c('In','Opc'))
        'MuSiC2 structurally collapsed' else '',
      inference='DESCRIPTIVE_ONLY: shared HD-PD controls; distinct covariates')
  }
  rbindlist(rows,fill=TRUE)
}
MAIN<-make_scenario('GSE125583'); VALID<-make_scenario('ROSMAP')
WR(MAIN,'13b_AD_GSE125583_main_direction.tsv')
WR(VALID,'13b_AD_ROSMAP_validation_direction.tsv')
# ROSMAP valida nao passa a ser o efeito principal. Comparacao de sinal por tipo e metodo.
RV<-merge(E[dataset=='GSE125583',.(method,cell_type,main_effect_pp=effect,main_fdr=fdr)],
 E[dataset=='ROSMAP',.(method,cell_type,validation_effect_pp=effect,validation_fdr=fdr)],
 by=c('method','cell_type'),all=TRUE)
RV[,same_sign:=ifelse(is.finite(main_effect_pp)&is.finite(validation_effect_pp),
                     sign(main_effect_pp)==sign(validation_effect_pp),NA)]
WR(RV,'13b_AD_main_vs_ROSMAP_validation.tsv')
if(file.exists(GROUPED)) {
 G<-fread(GROUPED,na.strings=c('NA','','NaN'))
 req<-c('dataset','method','group_type','scale','status','effect','ci_low','ci_high','p_hc3','fdr_bh_3groups')
 if(!all(req%in%names(G)))stop('11b incompleto')
 G<-G[scale=='percentage'&status=='PASS',.(dataset,method,group_type,effect,
              ci_low,ci_high,p_hc3,fdr=fdr_bh_3groups)]
 if(nrow(G)){
  WR(G,'13b_grouped_effects_by_cohort.tsv')
  grp<-list()
  for(mt in METHODS)for(gg in unique(G$group_type)) {
    ds<-c('GSE125583','GSE64810','GSE68719','ROSMAP')
    ef<-vapply(ds,function(d){z<-G[dataset==d&method==mt&group_type==gg];if(nrow(z)==1L)z$effect else NA_real_},numeric(1))
    grp[[length(grp)+1L]]<-data.table(method=mt,group_type=gg,AD_main_pp=ef[1],HD_pp=ef[2],PD_pp=ef[3],AD_validation_pp=ef[4],
      same_direction_main=if(all(is.finite(ef[1:3])))length(unique(sign(ef[1:3])))==1L else NA,
      descriptive_only=TRUE)
  }
  WR(rbindlist(grp),'13b_grouped_direction.tsv')
 }
}
writeLines(c('AD principal: GSE125583, AD validacao: ROSMAP, HD: GSE64810, PD: GSE68719.',
 'A robustez AD deve ser avaliada por direcao/magnitude/IC, nao so FDR.',
 'Compartilhamento de 44 controles HD-PD => nao combinar efeitos como independentes.',
 'Modelos tem covariaveis distintas; nenhum p-valor combinado/efeito pan-doenca.',
 'Concordancia de sinal nao prova acuracia, perda celular absoluta ou mecanismo compartilhado.',
 'dtangle principal e configuracao 08e; 08g e sensibilidade metodologica separada.'),file.path(OUT,'13b_LEIA_ANTES.txt'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'13b_sessionInfo.txt'))
cat('13b CONCLUIDO: ',OUT,'\n',sep='')
