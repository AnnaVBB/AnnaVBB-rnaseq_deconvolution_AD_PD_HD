#!/usr/bin/env Rscript
# ==============================================================================
# 13 | INTEGRACAO TRANSDIAGNOSTICA DESCRITIVA (SEM META-ANALISE INDEVIDA)
# Consome efeitos 11 ajustados e, se necessario, PD 11u NAO AJUSTADO.
# - Apresenta ROSMAP e GSE125583 como CENARIOS AD separados (nao escolhe AD).
# - HD/PD partilham 44 controles => nao trata os efeitos como independentes.
# - Ajustes diferentes/PD nao ajustado => nenhuma inferencia pan-doenca combinada.
# - Sem dados de todas as doencas, identifica apenas o que esta pendente.
# Executar: Rscript scripts/phase02_deconvolution/13_integracao_transdiagnostica_descritiva.R
# ==============================================================================
options(stringsAsFactors=FALSE)
if(!requireNamespace('data.table',quietly=TRUE))stop('Pacote data.table ausente')
suppressPackageStartupMessages(library(data.table))
ROOT <- normalizePath(Sys.getenv('DECONV_ROOT','.'),mustWork=TRUE)
P <- function(...) file.path(ROOT,...)
IN11 <- P('results','deconvolution','Integrated_Comparison','11_adjusted','11_adjusted_pp_HC3.tsv')
IN11U <- P('results','deconvolution','Integrated_Comparison','11_unadjusted_PD','11u_PD_unadjusted_effects.tsv')
OUT <- P('results','deconvolution','Integrated_Comparison','13_transdiagnostico_descritivo')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
WR <- function(x,filename)fwrite(x,file.path(OUT,filename),sep='\t',na='NA',quote=FALSE)
CT <- c('Ex','In','Ast','Oli','Opc','Mic','End','Per')
METHODS <- c('dtangle','BisqueRNA','MuSiC2')
D <- list()
if(file.exists(IN11)) {
  x <- fread(IN11,na.strings=c('NA','','NaN'))
  need <- c('dataset','disease','method','cell_type','scale','status','effect','ci_low','ci_high','p_hc3','fdr_bh_8types','covariates')
  if(!all(need %in% names(x)))stop('Arquivo 11 ajustado sem colunas obrigatorias')
  x <- x[scale=='percentage' & status=='PASS' & cell_type%in%CT & method%in%METHODS]
  if(nrow(x))D[[length(D)+1L]] <- x[,.(dataset,disease,method,cell_type,
    effect,ci_low,ci_high,p_hc3,fdr=fdr_bh_8types,
    adjustment='ADJUSTED',covariates,source='11_adjusted_pp_HC3')]
}
if(file.exists(IN11U)) {
  x <- fread(IN11U,na.strings=c('NA','','NaN'))
  need <- c('dataset','disease','method','cell_type','analysis_family','analysis_class',
            'status','effect','ci_low','ci_high','p_hc3','fdr_bh_family')
  if(!all(need%in%names(x)))stop('Arquivo PD 11u sem colunas obrigatorias')
  x <- x[analysis_family=='single_8types' & status=='PASS_UNADJUSTED' &
         analysis_class=='UNADJUSTED' & cell_type%in%CT & method%in%METHODS]
  if(nrow(x))D[[length(D)+1L]] <- x[,.(dataset,disease,method,cell_type,
    effect,ci_low,ci_high,p_hc3,fdr=fdr_bh_family,
    adjustment='UNADJUSTED',covariates='',source='11u_PD_unadjusted_effects')]
}
if(!length(D)) {
  WR(data.table(status='PENDING_INPUT',note='Sem resultados 11 ajustados ou 11u PD; execute os scripts primeiro'),
     '13_status.tsv')
  cat('13 sem entradas; consulte 13_status.tsv\n')
  quit(save='no',status=0)
}
E <- rbindlist(D,fill=TRUE)
# Se PD estiver ajustado em 11 futuramente, usar somente esse resultado e nao duplicar o 11u.
if(any(E$dataset=='GSE68719' & E$adjustment=='ADJUSTED'))
 E <- E[dataset!='GSE68719' | adjustment=='ADJUSTED']
if(anyDuplicated(E,by=c('dataset','method','cell_type')))
 stop('Entradas duplicadas para a mesma coorte/metodo/tipo')
if(any(!is.finite(E$effect)|!is.finite(E$ci_low)|!is.finite(E$ci_high)|
       E$ci_low>E$effect|E$ci_high<E$effect))
 stop('Intervalos de confianca ou efeitos invalidos nas entradas')
E[,`:=`(direction=ifelse(effect>0,'increase',ifelse(effect<0,'decrease','zero')),
        low_abundance_caution=cell_type%in%c('Mic','End','Per'),
        inference_note='Efeito especifico da coorte; nao representa numero de celulas absoluto')]
WR(E,'13_effects_by_cohort_method.tsv')
# Cenarios de Alzheimer sempre separados. Um AD nao e automaticamente a coorte principal.
ADs <- c('GSE125583','ROSMAP')
scenarios <- list()
for(mt in METHODS) for(ad in ADs) for(ct in CT) {
  opt <- c(ad,'GSE64810','GSE68719')
  block <- E[dataset%in%opt & method==mt & cell_type==ct]
  rowget <- function(ds,col) {
    xx <- block[dataset==ds]
    if(nrow(xx)!=1L)return(NA)
    xx[[col]][1L]
  }
  n_have <- sum(opt%in%block$dataset)
  effect <- vapply(opt,function(ds)as.numeric(rowget(ds,'effect')),numeric(1))
  names(effect) <- c('AD','HD','PD')
  adjusts <- vapply(opt,function(ds)as.character(rowget(ds,'adjustment')),character(1))
  names(adjusts) <- c('AD','HD','PD')
  if(n_have==3L){
    signs <- sign(effect)
    stable <- all(is.finite(signs)) && all(signs!=0) && length(unique(signs))==1L
    label <- if(stable)'SAME_DIRECTION_DESCRIPTIVE' else 'HETEROGENEOUS_DIRECTION'
    if(any(adjusts!='ADJUSTED')) label <- paste0(label,'_MIXED_ADJUSTMENT')
  } else {
    stable <- NA
    label <- 'INCOMPLETE_COHORT_EFFECTS'
  }
  scenarios[[length(scenarios)+1L]] <- data.table(
    AD_scenario=ad,method=mt,cell_type=ct,n_diseases_available=n_have,
    AD_effect_pp=effect['AD'],HD_effect_pp=effect['HD'],PD_effect_pp=effect['PD'],
    AD_adjustment=adjusts['AD'],HD_adjustment=adjusts['HD'],PD_adjustment=adjusts['PD'],
    same_direction=stable,status=label,
    note='COMPARACAO DESCRITIVA: HD/PD compartilham 44 controles; nao combinar p-valores ou assumir independencia')
}
S <- rbindlist(scenarios,fill=TRUE)
WR(S,'13_direction_by_AD_scenario.tsv')
WR(S[method=='dtangle'],'13_dtangle_primary_descriptive.tsv')
stat <- data.table(
  check=c('Inputs usados','Coortes presentes','Metodos','AD ainda nao eleita','Dependencia HD_PD','Inferencia pan-doenca'),
  value=c(paste(unique(E$source),collapse=';'),paste(unique(E$dataset),collapse=';'),
          paste(unique(E$method),collapse=';'),
          'GSE125583 e ROSMAP em cenarios separados',
          '44 controles comuns: nao independencia','NAO REALIZADA: comparacao descritiva apenas'))
WR(stat,'13_status.tsv')
# Figura auxiliar: efeitos absolutos (nao a composicao verdadeira); apenas se ggplot2 existir.
if(requireNamespace('ggplot2',quietly=TRUE)){
  gdat <- E[method=='dtangle']
  if(nrow(gdat)){
    gg <- ggplot2::ggplot(gdat,ggplot2::aes(x=effect,y=cell_type,color=dataset))+
      ggplot2::geom_vline(xintercept=0,linetype='dashed',color='grey65')+
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=ci_low,xmax=ci_high),height=.15,
                            position=ggplot2::position_dodge(width=.6))+
      ggplot2::geom_point(position=ggplot2::position_dodge(width=.6),size=2)+
      ggplot2::labs(title='Efeitos por coorte — dtangle 08e (descritivo)',
        subtitle='PD pode estar sem ajuste; HD e PD compartilham controles. NAO e meta-analise.',
        x='Diferenca caso - controle (pontos percentuais)',y='Tipo celular')+
      ggplot2::theme_minimal(base_size=11)
    ggplot2::ggsave(file.path(OUT,'13_dtangle_effects_by_cohort.png'),gg,
                    width=12,height=7,dpi=160,bg='white')
  }
}
writeLines(c('13 nao produz efeito pan-doenca combinado nem p-valor combinado.',
  'ROSMAP e GSE125583 permanecem cenarios AD distintos enquanto a selecao estiver pendente.',
  'O resultado PD com label UNADJUSTED nunca deve ser chamado de ajustado.',
  'HD e PD compartilham 44 controles; dependencia invalida combinacao ingenua por efeitos aleatorios.',
  'Concordancia de direcao e exploratoria e NAO demonstra processo biologico comum.',
  'Tipos raros (Mic, End, Per) exigem ressalvas; MuSiC2 In/Opc sem confiabilidade.',
  'Para inferencia transdiagnostica, planejar metodo que trate dependencia e heterogeneidade.'),
  file.path(OUT,'13_LEIA_ANTES.txt'))
cat('13 concluido: INTEGRACAO DESCRITIVA em ',OUT,'\n',sep='')
