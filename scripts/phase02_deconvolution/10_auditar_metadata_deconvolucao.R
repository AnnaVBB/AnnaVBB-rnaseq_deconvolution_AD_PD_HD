#!/usr/bin/env Rscript
# ==============================================================================
# 10 | AUDITORIA DE METADADOS E PLANEJAMENTO DE COVARIAVEIS
# TCC AD / PD / HD — NAO testa associacoes; NAO infere covariaveis ausentes.
#
# Entrada obrigatoria, por coorte:
#   data/processed/deconvolution/<COORTE>/metadata_deconvolution.tsv
# Extensao opcional, com provenance documentada:
#   data/processed/deconvolution/metadata_extension/<COORTE>_covariates.tsv
#   Chave sample_id obrigatoria; demais colunas devem ser NOVAS.
# Saidas: results/deconvolution/Integrated_Comparison/10_metadata/
#
# Uso: Rscript scripts/phase02_deconvolution/10_auditar_metadata_deconvolucao.R
# ==============================================================================
options(stringsAsFactors = FALSE)
if (!requireNamespace('data.table', quietly=TRUE)) stop('Instale data.table no renv, sem atualizar lock automaticamente.')
suppressPackageStartupMessages(library(data.table))
ROOT <- normalizePath(Sys.getenv('DECONV_ROOT', '.'), mustWork=TRUE)
P <- function(...) file.path(ROOT, ...)
OUT <- P('results','deconvolution','Integrated_Comparison','10_metadata')
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
DATASETS <- c('ROSMAP','GSE125583','GSE64810','GSE68719')
DISEASE <- c(ROSMAP='AD',GSE125583='AD',GSE64810='HD',GSE68719='PD')
EXPECTED_N <- c(ROSMAP=144L,GSE125583=289L,GSE64810=69L,GSE68719=73L)
CANDIDATES <- c('age','age_at_death','sex','gender','pmi','postmortem_interval',
 'rin','rna_integrity_number','batch','sequencing_batch','seq_batch',
 'region','brain_region','tissue_subtype','platform','study','donor_id','subject_id')
WR <- function(x,name) { fwrite(x,file.path(OUT,name),sep='\t',na='NA',quote=FALSE); invisible(x) }
check_keys <- function(x,fields,label) {
 if (!all(fields %in% names(x))) stop(label, ': missing columns ',paste(setdiff(fields,names(x)),collapse=','))
 for (k in fields) if (anyNA(x[[k]]) || any(trimws(as.character(x[[k]]))=='')) stop(label,': missing key ',k)
 if (anyDuplicated(x,by=fields)) stop(label,': duplicated keys ',paste(fields,collapse=','))
}
DICT <- list(); PROFILE <- list(); COV <- list(); COVER <- list(); STATUS <- list(); SHARED <- list()
for (ds in DATASETS) {
 cat('\n10 | ',ds,'\n',sep='')
 base_path <- P('data','processed','deconvolution',ds,'metadata_deconvolution.tsv')
 if (!file.exists(base_path)) stop('Missing metadata 05: ',base_path)
 x <- fread(base_path,na.strings=c('','NA','NaN','NULL'))
 check_keys(x,c('sample_id','group'),base_path)
 x[,sample_id:=trimws(as.character(sample_id))]
 x[,group:=trimws(as.character(group))]
 if (!setequal(unique(x$group),c('Control',unname(DISEASE[[ds]])))) stop(ds,': unexpected groups')
 if (nrow(x)!=EXPECTED_N[[ds]]) stop(ds,': expected ',EXPECTED_N[[ds]],' rows, got ',nrow(x))
 ext_path <- P('data','processed','deconvolution','metadata_extension',paste0(ds,'_covariates.tsv'))
 ext_used <- file.exists(ext_path)
 if (ext_used) {
   ext <- fread(ext_path,na.strings=c('','NA','NaN','NULL'))
   check_keys(ext,'sample_id',ext_path)
   ext[,sample_id:=trimws(as.character(sample_id))]
   extras <- setdiff(ext$sample_id,x$sample_id)
   if (length(extras)) stop(ds,': extension has sample IDs not present in 05: ',paste(head(extras,5),collapse=','))
   collisions <- intersect(setdiff(names(ext),'sample_id'),names(x))
   if (length(collisions)) stop(ds,': extension duplicates existing columns; rename or reconcile manually: ',paste(collisions,collapse=','))
   x <- merge(x,ext,by='sample_id',all.x=TRUE,sort=FALSE)
   check_keys(x,'sample_id',paste0(ds,' merged'))
 }
 x[, dataset := ds]
 x[, disease := unname(DISEASE[[ds]])]
 setcolorder(x,c('dataset','disease','sample_id','group',setdiff(names(x),c('dataset','disease','sample_id','group'))))
 WR(x,paste0('10_metadata_',ds,'_audited.tsv'))
 COVER[[ds]] <- data.table(dataset=ds,n_total=nrow(x),n_control=sum(x$group=='Control'),
    n_case=sum(x$group==DISEASE[[ds]]),n_fields=ncol(x),extension_used=ext_used,
    extension_file=if(ext_used) ext_path else NA_character_)
 cols <- setdiff(names(x),c('dataset','disease','sample_id','group'))
 for (nm in cols) {
   v <- x[[nm]]
   nonmissing <- !is.na(v) & trimws(as.character(v))!=''
   nu <- uniqueN(as.character(v[nonmissing])); lev <- if(nu<=10)paste(sort(unique(as.character(v[nonmissing]))),collapse=';') else NA_character_
   cand <- tolower(nm) %in% CANDIDATES || grepl('age|sex|gender|rin|batch|pmi|region|quality|qc|library|donor|subject',tolower(nm))
   DICT[[length(DICT)+1L]] <- data.table(dataset=ds,column=nm,R_type=class(v)[1],
     n_observed=sum(nonmissing),n_missing=sum(!nonmissing),n_distinct=nu,
     complete_pct=round(100*mean(nonmissing),2),candidate_for_review=cand,
     levels_if_10_or_fewer=lev)
   PROFILE[[length(PROFILE)+1L]] <- x[,.(n=.N,n_observed=sum(!is.na(get(nm)) & trimws(as.character(get(nm)))!='')),
                                      by=.(group)][,`:=`(dataset=ds,column=nm)]
 }
 STATUS[[ds]] <- data.table(dataset=ds,status='PASS_METADATA_AUDITED',
   note='Inventory only; NOT an adjusted case/control model')
 cat(sprintf('  n=%d, fields=%d, extension=%s\n',nrow(x),ncol(x),ext_used))
}
WR(rbindlist(COVER,fill=TRUE),'10_dataset_coverage.tsv')
WR(if(length(DICT))rbindlist(DICT,fill=TRUE) else data.table(dataset=character(),column=character()),
   '10_covariate_dictionary.tsv')
WR(if(length(PROFILE))rbindlist(PROFILE,fill=TRUE) else data.table(),
   '10_covariate_coverage_by_group.tsv')
WR(rbindlist(STATUS),'10_status.tsv')
# MODELO DE CONFIGURACAO: NAO preencher automaticamente termos com base nos resultados.
# O pesquisador e orientador devem deliberar o conjunto de covariaveis POR COORTE.
TEMPLATE <- data.table(dataset=DATASETS,approved='NO',
  numeric_covariates='',categorical_covariates='',
  cluster_column='',
  justification='PENDING: review technical/clinical availability and causal role',
  provenance='PENDING: original cohort metadata / publication',
  notes='Do not use severity/outcomes as automatic baseline confounders')
# Evitar sobrescrever um plano que o usuario possa ter preenchido.
planpath <- P('data','processed','deconvolution','10_approved_model_plan.tsv')
if (!file.exists(planpath)) WR(TEMPLATE,'10_MODEL_PLAN_TEMPLATE.tsv')
shared_info <- data.table(note=c('HD and PD reuse 44 controls: NEVER treat as independent external validation',
  'The metadata 05 cannot establish brain regions if region is missing',
  'No fields were silently imputed, converted, or pulled from the web',
  'Fill an approved model plan only after reviewing columns and missingness'))
WR(shared_info,'10_INTERPRETATION_NOTES.tsv')
writeLines(capture.output(sessionInfo()),file.path(OUT,'10_sessionInfo.txt'))
cat('\n10 | CONCLUIDO (auditoria; sem modelos). Saidas: ',OUT,'\n',sep='')
