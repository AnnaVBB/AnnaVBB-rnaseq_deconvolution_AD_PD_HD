#!/usr/bin/env Rscript
# ==============================================================================
# 10d | PREPARAR EXTENSOES CANDIDATAS PARA HD, PD E ROSMAP (REVISAO OBRIGATORIA)
# AD/PD/HD - Fase 02 | Nao modifica arquivos de analise nem aprova modelos.
# Executar, na raiz do repositorio:
#   Rscript scripts/phase02_deconvolution/10d_preparar_metadados_HD_PD_ROSMAP.R
# Pre-requisitos: 10b v2, 10c e metadata ROSMAP de fase01c, quando disponivel.
# ==============================================================================
options(stringsAsFactors = FALSE)
if (!requireNamespace('data.table', quietly = TRUE)) stop('Falta data.table')
suppressPackageStartupMessages(library(data.table))
ROOT <- normalizePath(Sys.getenv('DECONV_ROOT', '.'), mustWork = TRUE)
P <- function(...) file.path(ROOT, ...)
BASE <- P('data', 'processed', 'deconvolution')
IN10C <- P('results','deconvolution','Integrated_Comparison','10_metadata','10c_HD_PD')
OUT <- P('results','deconvolution','Integrated_Comparison','10_metadata','10d_preparacao')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, file) fwrite(x, file.path(OUT,file), sep='\t', na='NA', quote=FALSE)
not_blank <- function(z) !is.na(z) & nzchar(trimws(as.character(z)))
assert_unique <- function(d, cols, label) {
  if (!all(cols %in% names(d))) stop(label, ': faltam ',paste(setdiff(cols,names(d)),collapse=', '))
  if (any(!vapply(d[, ..cols], function(z) all(not_blank(z)), logical(1))))
    stop(label, ': chaves ausentes')
  if (anyDuplicated(d, by = cols)) stop(label, ': chaves duplicadas')
}
read_required <- function(path, cols, label) {
  if (!file.exists(path)) stop(label, ': arquivo ausente: ',path)
  d <- fread(path, na.strings=c('NA','NaN','','NULL'))
  assert_unique(d, 'sample_id', label)
  if (!all(cols %in% names(d))) stop(label, ': faltam colunas ',paste(setdiff(cols,names(d)), collapse=', '))
  d[, sample_id := trimws(as.character(sample_id))]
  d
}
read_meta <- function(ds, expected_n, case_label) {
  d <- read_required(file.path(BASE, ds, 'metadata_deconvolution.tsv'), c('sample_id','group'), ds)
  d[, group := trimws(as.character(group))]
  if(nrow(d)!=expected_n || !setequal(unique(d$group),c('Control',case_label)))
    stop(ds, ': tamanho ou grupos inesperados')
  d[, .(sample_id,group)]
}
match_exact <- function(meta, ext, label) {
  if(!setequal(meta$sample_id,ext$sample_id))
    stop(label, ': extensao nao corresponde exatamente ao universo de amostras da etapa 05')
  ext <- ext[match(meta$sample_id,ext$sample_id)]
  if(!identical(meta$sample_id,ext$sample_id)) stop(label, ': erro no alinhamento')
  ext
}
number_or_na <- function(v, label) {
  s <- trimws(as.character(v))
  blank <- is.na(v)|!nzchar(s)
  ok <- blank | grepl('^[0-9]+(\\.[0-9]+)?$',s)
  if(any(!ok)) stop(label, ': valores numericos nao interpretaveis: ',paste(head(unique(s[!ok]),5),collapse='; '))
  x <- suppressWarnings(as.numeric(s)); x[blank] <- NA_real_
  if(any(!is.finite(x[!is.na(x)])))stop(label,': infinito ou NaN')
  x
}
coverage <- list(); status <- list(); add_status <- function(ds, result, note) {
  status[[length(status)+1L]] <<- data.table(dataset=ds,status=result,detail=note)
  cat(ds,' | ',result,' | ',note,'\n',sep='')
}
add_coverage <- function(ds, meta, ext, vars) {
  z <- copy(meta)
  for(v in vars) z[, (v) := ext[[v]][match(sample_id,ext$sample_id)]]
  coverage[[length(coverage)+1L]] <<- rbindlist(lapply(vars, function(v)
    z[,.(n=.N,available=sum(!is.na(get(v)) & not_blank(get(v)))),by=group][,
      `:=`(dataset=ds,variable=v,missing=n-available)]))
  cc <- z[, .(group, complete=complete.cases(z[, ..vars]))]
  add_status(ds,'COBERTURA',paste(paste(cc[,.(n=.N,complete=sum(complete)),by=group][,
    paste(group,complete,'/',n)],collapse='; '),'casos completos para ',paste(vars,collapse='+'),sep=''))
}
HD <- read_meta('GSE64810',69L,'HD')
PD <- read_meta('GSE68719',73L,'PD')
RO <- read_meta('ROSMAP',144L,'AD')

# HD + PD: aceita apenas os candidatos produzidos pelo 10c, nunca transferencias por prefixo.
shared_file <- file.path(IN10C,'10c_PD_HD_shared_controls_crosswalk.tsv')
hd_file <- file.path(IN10C,'10c_GSE64810_covariates_REVIEW_FIRST.tsv')
pd_file <- file.path(IN10C,'10c_GSE68719_covariates_REVIEW_FIRST.tsv')
if(all(file.exists(c(shared_file,hd_file,pd_file)))) {
  cm <- fread(shared_file,na.strings=c('NA','','NaN','NULL'))
  if(!all(c('sample_id','status_match','age_conflict','sex_conflict') %in% names(cm)))
    stop('Crosswalk 10c incompleto')
  if(nrow(cm)!=44L || anyDuplicated(cm,by='sample_id') ||
     any(cm$status_match!='EXACT_CONTROL_IN_BOTH') ||
     any(cm$age_conflict %in% TRUE) || any(cm$sex_conflict %in% TRUE) ||
     !setequal(cm$sample_id,PD[group=='Control',sample_id]) ||
     !all(cm$sample_id %in% HD[group=='Control',sample_id]))
    stop('Controles compartilhados NAO confirmados pelo cruzamento exato')
  add_status('HD_PD','PASS_CROSSWALK','44 controles coincidentes, rotulados Control em ambas as coortes')
  hd_raw <- read_required(hd_file,c('rin','pmi'),'10c_HD')
  pd_raw <- read_required(pd_file,c('age','sex','rin','pmi'),'10c_PD')
  hd_raw <- match_exact(HD,hd_raw,'HD')
  pd_raw <- match_exact(PD,pd_raw,'PD')
  HD_ext <- hd_raw[, .(sample_id,
                        rin=number_or_na(rin,'HD rin'),
                        pmi=number_or_na(pmi,'HD pmi'))]
  PD_ext <- pd_raw[, .(sample_id,
                        age=number_or_na(age,'PD age'),
                        sex=as.character(sex),
                        rin=number_or_na(rin,'PD rin'),
                        pmi=number_or_na(pmi,'PD pmi'))]
  PD_ext[!not_blank(sex),sex:=NA_character_]
  if(any(!is.na(PD_ext$sex) & !PD_ext$sex %in% c('M','F')))
    stop('PD: codificacao de sexo inesperada')
  write_tsv(HD_ext,'10d_GSE64810_covariates_REVIEW_FIRST.tsv')
  write_tsv(PD_ext,'10d_GSE68719_covariates_REVIEW_FIRST.tsv')
  add_coverage('GSE64810',HD,HD_ext,c('rin','pmi'))
  add_coverage('GSE68719',PD,PD_ext,c('age','sex','rin','pmi'))
  add_status('GSE64810','CANDIDATE_READY','Revisar unidades de RIN/PMI, arquivos originais e perdas antes de incorporar')
  add_status('GSE68719','EXPLORATORY_ONLY','Variaveis disponíveis em grupos diferentes; nao aprovar modelo ajustado')
} else {
  add_status('HD_PD','PENDING_10C',paste('Execute 10c; ausentes:',
     paste(basename(c(shared_file,hd_file,pd_file)[!file.exists(c(shared_file,hd_file,pd_file))]),collapse=', ')))
}

# ROSMAP: reutilizar a preparacao ORIGINAL da fase 01c e conferir identidade e diagnostico.
ro_candidates <- c(
 P('results','candidate_ROSMAP','metadata','metadata_ROSMAP.tsv'),
 P('results','candidate_ROSMAP','metadata','ROSMAP_metadata_main.tsv'),
 P('results','candidate_ROSMAP','qc','metadata_ROSMAP_main.tsv'),
 P('results','candidate_ROSMAP','qc','ROSMAP_metadata_main.tsv')
)
found <- ro_candidates[file.exists(ro_candidates)]
if(length(found)) {
  src <- found[1L]
  m <- read_required(src,c('group','msex','RIN','pmi','sequencingBatch'),'ROSMAP fase01c')
  if(!setequal(m$sample_id,RO$sample_id))
    stop('ROSMAP: IDs da fase01c nao coincidem com as 144 amostras da etapa 05')
  m <- m[match(RO$sample_id,m$sample_id)]
  if(!identical(as.character(m$group),as.character(RO$group)))
    stop('ROSMAP: divergencia de grupos entre fase01c e etapa 05')
  sexcoded <- trimws(as.character(m$msex))
  sexcoded[!not_blank(sexcoded)] <- NA_character_
  if(any(!is.na(sexcoded) & !sexcoded %in% c('0','1')))
    stop('ROSMAP: msex contem valores inesperados; rever codificacao sem presumir o significado de 0/1')
  batch <- trimws(as.character(m$sequencingBatch))
  batch[!not_blank(batch)] <- NA_character_
  RO_ext <- data.table(sample_id=RO$sample_id,
                       rin=number_or_na(m$RIN,'ROSMAP RIN'),
                       pmi=number_or_na(m$pmi,'ROSMAP PMI'),
                       msex=sexcoded,
                       sequencingBatch=batch)
  write_tsv(RO_ext,'10d_ROSMAP_covariates_REVIEW_FIRST.tsv')
  add_coverage('ROSMAP',RO,RO_ext,c('rin','pmi','msex','sequencingBatch'))
  add_status('ROSMAP','CANDIDATE_READY',paste('Fonte:',src,'; msex 0/1 preservado; idade 90+ nao inferida'))
} else {
  add_status('ROSMAP','SOURCE_NOT_FOUND','Rode/recupere fase01c/01_prepare_metadata.R e verifique results/candidate_ROSMAP/metadata/metadata_ROSMAP.tsv')
}

if(length(coverage)) write_tsv(rbindlist(coverage,fill=TRUE),'10d_covariate_coverage_by_group.tsv')
write_tsv(rbindlist(status,fill=TRUE),'10d_status.tsv')
# Proposta: nao altera o plano APROVADO. Colunas e texto podem ser refinados na revisao.
planfile <- file.path(BASE,'10_approved_model_plan.tsv')
if(file.exists(planfile)) {
  p <- fread(planfile,na.strings=c('NA','NaN'))
  if(nrow(p)!=4L||anyDuplicated(p,by='dataset'))stop('Plano anterior invalido')
  if('GSE64810'%in%p$dataset){
    p[dataset=='GSE64810',`:=`(approved='NO',numeric_covariates='rin,pmi',categorical_covariates='',
      justification='PROPOSTA: controlar qualidade RNA e intervalo post-mortem; revisar unidades, cobertura e distribuicao',
      provenance='GSE64810 GEO SOFT: extraido 10b e cruzado 10c; revisar antes da aprovacao')]
  }
  if('ROSMAP'%in%p$dataset){
    p[dataset=='ROSMAP',`:=`(approved='NO',numeric_covariates='rin,pmi',categorical_covariates='msex,sequencingBatch',
      justification='PROPOSTA: covariaveis tecnicas e sexo codificado; verificar confundimento batch x grupo e n final',
      provenance='ROSMAP metadata da fase01c; pareamento exato com etapa05 e regras de acesso')]
  }
  if('GSE68719'%in%p$dataset){
    p[dataset=='GSE68719',`:=`(approved='NO',numeric_covariates='',categorical_covariates='',
      justification='NAO APROVAR: covariaveis existentes nao cobrem simultaneamente casos e controles',
      provenance='GSE68719 e 44 controles compartilhados com GSE64810')]
  }
  write_tsv(p,'10d_model_plan_PROPOSTA_NAO_APROVADA.tsv')
}
writeLines(c(
 'NENHUM arquivo em data/processed foi modificado.',
 'Revisar 10d_status e 10d_covariate_coverage_by_group.',
 'HD: apos revisar RIN/PMI e unidades, copiar somente o candidato HD para data/processed/deconvolution/metadata_extension/GSE64810_covariates.tsv.',
 'ROSMAP: apos verificar fonte, acesso, diagnostico, batch e perdas, copiar candidato ROSMAP para metadata_extension/ROSMAP_covariates.tsv.',
 'PD: candidato e exploratorio; NAO copiar como se houvesse covariaveis em ambos os grupos.',
 'Executar 10_auditar_metadata_deconvolucao.R novamente apos incorporar as extensoes revisadas.',
 'Somente apos auditoria, ajustar manualmente 10_approved_model_plan.tsv para YES nas coortes liberadas.',
 'Nao usar o plano PROPOSTA_NAO_APROVADA como plano aprovado sem conferir os resultados.'
 ),file.path(OUT,'10d_LEIA_ANTES.txt'))
cat('\n10d concluido. Arquivos para revisao: ',OUT,'\n',sep='')
