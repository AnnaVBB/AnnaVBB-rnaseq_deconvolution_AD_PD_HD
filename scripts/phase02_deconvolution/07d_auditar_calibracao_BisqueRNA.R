#!/usr/bin/env Rscript
# ==============================================================================
# 07d — AUDITORIA DE CALIBRACAO E GENERALIZACAO DO BISQUERNA
# Projeto: deconvolucao AD / PD / HD | referencia Mathys 2019 pathology-free
# ==============================================================================
# OBJETIVOS (tres experimentos SEM alterar a etapa 07):
#   A. Alterar a composicao de NUCLEOS da referencia (bulk fixo).
#      Permite testar ancoragem a frequencias de referencia, com replicas.
#   B. Reutilizar as MESMAS 44 amostras-controle, com counts E libraries
#      identicos, incluindo casos HD, PD, ambos ou nenhum.
#      Permite isolar a dependencia do contexto bulk.
#   C. Validacao sintetica leave-one-donor-out (LODO): misturas de composicao
#      conhecida geradas por um doador excluido do treinamento; marcadores
#      recalculados APENAS nos doadores de treino.
#
# ALERTAS CIENTIFICOS:
#   * Saidas de A e B sao experimentos de sensibilidade, NAO precisao real.
#   * Misturas C sao da mesma plataforma snRNA-seq, nao bulk real conhecido.
#   * Subamostragem em A muda distribuicao de doadores por tipo, alem da
#     proporcao agregada; replicas e uma base tamanho-igual sao essenciais.
#   * Os valores 'verdadeiros' em C de fracao de NUCLEOS e de contribuicao
#     de RNA sao diferentes. Ambos sao avaliados, SEM equipara-los.
#   * Nao se escolhe uma configuracao por plausibilidade ou p-valor.
#   * RDS originais de 06, 07, 08 NUNCA sao escritos.
#
# Versao corrigida: 07d C usa duas chamadas consecutivas a data.table :=
# para calcular erros e erros absolutos, corrigindo objeto nao encontrado.
# Requisitos: R + BisqueRNA + Biobase + Matrix + data.table + limSolve.
# Nenhuma instalacao, nenhuma alteracao do renv.lock, sem testes inferenciais.
#
# EXEMPLOS (executar na raiz do repositorio):
#   Rscript scripts/phase02_deconvolution/07d_auditar_calibracao_BisqueRNA.R
#   BISQUE07D_PARTS=A Rscript scripts/phase02_deconvolution/07d_auditar_calibracao_BisqueRNA.R
#   BISQUE07D_PARTS=B Rscript scripts/phase02_deconvolution/07d_auditar_calibracao_BisqueRNA.R
#   BISQUE07D_PARTS=C Rscript scripts/phase02_deconvolution/07d_auditar_calibracao_BisqueRNA.R
#   BISQUE07D_A_REPS=3 BISQUE07D_C_HOLDOUTS=4 Rscript ...
#
# Saida: results/deconvolution/BisqueRNA_Audit/07d/{A,B,C,plots}
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. CONFIGURACAO
# ------------------------------------------------------------------------------
TYPES <- c('Ex','In','Ast','Oli','Opc','Mic','End','Per')
DATASETS <- c('ROSMAP','GSE125583','GSE64810','GSE68719')
REF_BASE <- 'data/reference/mathys2019/phase02_prepared/reference'
REF_COUNTS_PATH <- file.path(REF_BASE,'mathys2019_pathology_free_geneid_counts_sparse.rds')
REF_META_PATH <- file.path(REF_BASE,'mathys2019_pathology_free_cell_metadata_aligned.tsv')
BULK_ROOT <- 'data/processed/deconvolution'
ORIG07_ROOT <- 'results/deconvolution/BisqueRNA'
OUT <- 'results/deconvolution/BisqueRNA_Audit/07d'

raw_parts <- trimws(Sys.getenv('BISQUE07D_PARTS',unset='A,B,C'))
PARTS <- unique(trimws(strsplit(gsub(';',',',raw_parts,fixed=TRUE),',',fixed=TRUE)[[1]]))
PARTS <- toupper(PARTS[nzchar(PARTS)])
if (!length(PARTS) || any(!PARTS %in% c('A','B','C'))) {
  stop('BISQUE07D_PARTS deve listar A, B e/ou C separados por virgula.')
}
A_DATASET <- trimws(Sys.getenv('BISQUE07D_A_DATASET',unset='GSE64810'))
if (!A_DATASET %in% DATASETS) stop('BISQUE07D_A_DATASET invalido.')
A_N <- as.integer(Sys.getenv('BISQUE07D_A_NUCLEI',unset='6000'))
A_REPS <- as.integer(Sys.getenv('BISQUE07D_A_REPS',unset='2'))
C_HOLDOUTS <- as.integer(Sys.getenv('BISQUE07D_C_HOLDOUTS',unset='4'))
REUSE <- identical(Sys.getenv('BISQUE07D_REUSE',unset='1'),'1')
if (!is.finite(A_N) || A_N<1000L || A_N>20000L ||
    !is.finite(A_REPS) || A_REPS<1L || A_REPS>10L ||
    !is.finite(C_HOLDOUTS) || C_HOLDOUTS<1L || C_HOLDOUTS>12L) {
  stop('Parametros A_N, A_REPS ou C_HOLDOUTS fora dos limites permitidos.')
}
SEED <- 20261022L
TOL <- 1e-6
PACKAGE_NAMES <- c('BisqueRNA','Biobase','Matrix','data.table','limSolve')
missing <- PACKAGE_NAMES[!vapply(PACKAGE_NAMES,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop('Pacotes ausentes: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Biobase)
})
for (sub in c('A','B','C','plots')) dir.create(file.path(OUT,sub),recursive=TRUE,showWarnings=FALSE)

log_msg <- function(...) cat(sprintf(...),'\n',sep='')
write_tab <- function(x,p) {
  dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
  data.table::fwrite(x,p,sep='\t',na='NA',quote=FALSE)
}
read_req <- function(path) {
  if (!file.exists(path)) stop('Arquivo requerido ausente: ',path)
  path
}
check_matrix <- function(x,label) {
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop(label,': IDs de linhas/colunas ausentes ou duplicados')
  }
  if (anyNA(x) || any(!is.finite(x)) || min(x)<0) stop(label,': contagens invalidas')
  invisible(TRUE)
}
canon_meta <- function(meta,ids,label) {
  if (!all(c('TAG','cell_type','donor_id') %in% names(meta))) stop('Metadata sem TAG/donor_id/cell_type')
  if (anyDuplicated(meta$TAG)) stop('TAG duplicado na referencia')
  ix <- match(ids,as.character(meta$TAG))
  if (anyNA(ix)) stop('Nucleos sem metadados em ',label)
  ans <- data.table::copy(meta[ix])
  if (!identical(as.character(ans$TAG),ids)) stop('Metadata nao alinhado: ',label)
  ans[,`:=`(cell_type=as.character(cell_type),donor_id=as.character(donor_id),TAG=as.character(TAG))]
  if (anyNA(ans$cell_type) || anyNA(ans$donor_id) ||
      !all(ans$cell_type %in% TYPES)) stop('Metadata celular invalido')
  ans
}

# ------------------------------------------------------------------------------
# 1. ENTRADAS E METADADOS FIXOS
# ------------------------------------------------------------------------------
log_msg('07d | Iniciando auditorias [%s]',paste(PARTS,collapse=','))
ref_all <- readRDS(read_req(REF_COUNTS_PATH))
if (!inherits(ref_all,'dgCMatrix')) ref_all <- as(ref_all,'dgCMatrix')
check_matrix(ref_all,'Mathys')
ref_meta <- canon_meta(fread(read_req(REF_META_PATH)),colnames(ref_all),'Mathys')
if (ncol(ref_all)!=35140L || nrow(ref_all)!=17394L ||
    uniqueN(ref_meta$donor_id)!=24L || !setequal(unique(ref_meta$cell_type),TYPES)) {
  stop('Dimensoes/donors/tipos da referencia diferentes da fase 02.')
}
ref_lib <- as.numeric(Matrix::colSums(ref_all))
if (any(ref_lib<=0)) stop('Library size zero no snRNA-seq')
ref_donor_ids <- sort(unique(ref_meta$donor_id))

load_bulk <- function(ds) {
  root <- file.path(BULK_ROOT,ds)
  p_counts <- read_req(file.path(root,'bulk_counts_top300.rds'))
  p_meta <- read_req(file.path(root,'metadata_deconvolution.tsv'))
  p_lib <- read_req(file.path(root,'library_sizes.tsv'))
  x <- as.matrix(readRDS(p_counts)); storage.mode(x) <- 'double'
  check_matrix(x,paste('bulk',ds))
  m <- fread(p_meta)
  l <- fread(p_lib)
  if (!all(c('sample_id','group') %in% names(m)) ||
      !all(c('sample_id','library_size_full_input') %in% names(l))) {
    stop(ds,': metadados/library_sizes incompletos')
  }
  if (anyDuplicated(m$sample_id) || anyDuplicated(l$sample_id)) stop(ds,': IDs bulk duplicados')
  im <- match(colnames(x),as.character(m$sample_id))
  il <- match(colnames(x),as.character(l$sample_id))
  if (anyNA(im) || anyNA(il)) stop(ds,': metadata ou library size sem todas amostras')
  m <- m[im]; libs <- as.numeric(l$library_size_full_input[il])
  if (any(!is.finite(libs)) || any(libs<=0) ||
      any(as.numeric(colSums(x))>libs+1e-5)) stop(ds,': library sizes invalidas')
  list(counts=x,meta=m,lib=libs,dataset=ds)
}

# ------------------------------------------------------------------------------
# 2. ADAPTADOR BISQUERNA — reproduz normalizacao do script 07
# ------------------------------------------------------------------------------
# Bisque old.cpm=FALSE calcula CPM antes de selecionar markers. Para evitar
# uma matriz densa 17394 x 35140, usamos 1 linha CONTABIL de complemento
# ate o tamanho TOTAL da biblioteca. A linha nunca entra em genes.used.
# Verificar as mesmas invariantes da etapa 07.
fit_bisque <- function(ref_counts,ref_metadata,ref_libraries,
                       bulk_counts,bulk_libraries,marker_genes,run_name,
                       folder,context) {
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  genes <- intersect(as.character(marker_genes),intersect(rownames(ref_counts),rownames(bulk_counts)))
  genes <- genes[Matrix::rowSums(ref_counts[genes,,drop=FALSE])>0 &
                   rowSums(bulk_counts[genes,,drop=FALSE])>0]
  if (length(genes)<500L) stop(run_name,': menos de 500 genes para ajuste')
  if (length(ref_libraries)!=ncol(ref_counts) ||
      length(bulk_libraries)!=ncol(bulk_counts) ||
      nrow(ref_metadata)!=ncol(ref_counts)) stop(run_name,': dimensoes invalidas')
  if (!identical(ref_metadata$TAG,colnames(ref_counts))) stop('TAG/reference desalinhados')
  if (any(!ref_metadata$cell_type %in% TYPES)) stop('Tipo celular desconhecido')
  if (uniqueN(ref_metadata$donor_id)<3L) stop('Bisque requer pelo menos 3 doadores de referencia')
  sc <- as.matrix(ref_counts[genes,,drop=FALSE]); storage.mode(sc) <- 'double'
  bu <- as.matrix(bulk_counts[genes,,drop=FALSE]); storage.mode(bu) <- 'double'
  sc_rest <- ref_libraries-colSums(sc)
  bu_rest <- bulk_libraries-colSums(bu)
  if (min(sc_rest) < -1e-4 || min(bu_rest) < -1e-4) {
    stop(run_name,': denominadores completos menores do que os genes selecionados')
  }
  sc_rest[abs(sc_rest)<1e-7] <- 0
  bu_rest[abs(bu_rest)<1e-7] <- 0
  sc <- rbind(sc,'__LIBRARY_REMAINDER__'=sc_rest)
  bu <- rbind(bu,'__LIBRARY_REMAINDER__'=bu_rest)
  if (max(abs(colSums(sc)-ref_libraries))>1e-4 ||
      max(abs(colSums(bu)-bulk_libraries))>1e-4) {
    stop(run_name,': conservacao das bibliotecas falhou')
  }
  sc_pd <- data.frame(cellType=factor(ref_metadata$cell_type,levels=TYPES),
                      SubjectName=factor(ref_metadata$donor_id),
                      row.names=colnames(sc),check.names=FALSE)
  bl_pd <- data.frame(sample_id=colnames(bu),row.names=colnames(bu),check.names=FALSE)
  sc_eset <- Biobase::ExpressionSet(
    assayData=sc,phenoData=Biobase::AnnotatedDataFrame(sc_pd))
  bl_eset <- Biobase::ExpressionSet(
    assayData=bu,phenoData=Biobase::AnnotatedDataFrame(bl_pd))
  rm(sc,bu); invisible(gc(FALSE))
  time_start <- proc.time()[['elapsed']]
  fit <- BisqueRNA::ReferenceBasedDecomposition(
    bulk.eset=bl_eset,sc.eset=sc_eset,markers=genes,
    cell.types='cellType',subject.names='SubjectName',
    use.overlap=FALSE,old.cpm=FALSE,verbose=FALSE)
  elapsed <- proc.time()[['elapsed']]-time_start
  p <- as.matrix(fit$bulk.props)
  if (is.null(rownames(p)) || is.null(colnames(p)) ||
      !setequal(rownames(p),TYPES) ||
      !setequal(colnames(p),colnames(bulk_counts))) stop(run_name,': IDs na saida incorretos')
  p <- p[TYPES,colnames(bulk_counts),drop=FALSE]
  if (any(!is.finite(p)) || min(p) < -1e-8 || max(p)>1+1e-8 ||
      max(abs(colSums(p)-1))>TOL) stop(run_name,': proporcoes invalidas')
  if ('__LIBRARY_REMAINDER__' %in% as.character(fit$genes.used)) {
    stop(run_name,': linha tecnica entrou na estimacao!')
  }
  if (!is.null(fit$sc.props)) {
    q <- as.matrix(fit$sc.props)
    if (!setequal(rownames(q),TYPES)) stop('sc.props nao preservou os oito tipos')
    q <- q[TYPES,,drop=FALSE]
  } else q <- matrix(NA_real_,8,0,dimnames=list(TYPES,NULL))
  # Armazenar resultado auditavel antes de qualquer grafico posterior.
  saveRDS(list(proportions=p,sc_props=q,genes=fit$genes.used,
               elapsed=elapsed,version=as.character(packageVersion('BisqueRNA')),
               run_name=run_name,context=context),
          file.path(folder,'fit_audit.rds'),compress='xz')
  wide <- data.table(sample_id=colnames(p)); wide <- cbind(wide,as.data.table(t(p)))
  write_tab(wide,file.path(folder,'proportions_wide.tsv'))
  data.table(run=run_name,genes_used=length(fit$genes.used),
             n_bulk=ncol(p),n_ref=ncol(ref_counts),
             n_ref_donors=uniqueN(ref_metadata$donor_id),
             elapsed_seconds=elapsed,max_sum_error=max(abs(colSums(p)-1)),
             bisque_version=as.character(packageVersion('BisqueRNA')))
}

# Salvar partes ja concluidas e continuar outras execucoes sem rerodagem.
quick_hash <- function(v) {
  p <- tempfile('bisque07d_hash_')
  on.exit(unlink(p))
  writeLines(as.character(v),p,useBytes=TRUE)
  unname(tools::md5sum(p)[[1]])
}
source_signature <- function(paths) {
  info <- file.info(paths)
  paste(paths,info$size,format(info$mtime,'%Y-%m-%dT%H:%M:%S'),collapse=';')
}
SOURCE_STAMP <- source_signature(c(REF_COUNTS_PATH,REF_META_PATH))
run_fitted <- function(key,args,folder) {
  # Checkpoint nao depende apenas do nome: identifica colunas, nucleos,
  # marcadores, libraries e arquivos-base. Para refazer, USE REUSE=0.
  args$context <- paste0(args$context,'|',SOURCE_STAMP,'|',
    quick_hash(c(colnames(args$ref_counts),colnames(args$bulk_counts),
                 rownames(args$bulk_counts),args$marker_genes,
                 format(args$ref_libraries,digits=15),
                 format(args$bulk_libraries,digits=15),
                 format(sum(args$bulk_counts),digits=17),
                 format(sum(args$ref_counts),digits=17))))
  cache <- file.path(folder,'fit_audit.rds')
  if (REUSE && file.exists(cache)) {
    ans <- tryCatch(readRDS(cache),error=function(e) NULL)
    # Reuse somente se chave e parametros principais tiverem identidade.
    if (!is.null(ans) && identical(ans$run_name,key) &&
        identical(ans$version,as.character(packageVersion('BisqueRNA'))) &&
        identical(ans$context,args$context) &&
        identical(colnames(ans$proportions),colnames(args$bulk_counts)) &&
        identical(rownames(ans$proportions),TYPES) &&
        isTRUE(all.equal(unname(colSums(ans$proportions)),rep(1,ncol(ans$proportions)),tolerance=TOL))) {
      log_msg('    [checkpoint] %s',key)
      return(list(proportions=ans$proportions,sc_props=ans$sc_props,
                  genes_used=length(ans$genes),elapsed_seconds=0,
                  cached=TRUE))
    }
  }
  log_msg('    [ajuste] %s',key)
  tmp <- do.call(fit_bisque,c(args,list(run_name=key,folder=folder)))
  ans <- readRDS(cache)
  list(proportions=ans$proportions,sc_props=ans$sc_props,
       genes_used=tmp$genes_used,elapsed_seconds=tmp$elapsed_seconds,cached=FALSE)
}

statuses <- list()
record_status <- function(part,run,state,note='') {
  statuses[[length(statuses)+1L]] <<- data.table(
    part=part,run=run,status=state,note=as.character(note))
  write_tab(rbindlist(statuses,fill=TRUE),file.path(OUT,'07d_status.tsv'))
}
safe_call <- function(part,run,fun) {
  tryCatch({
    out <- fun()
    record_status(part,run,'PASS')
    out
  },error=function(e) {
    record_status(part,run,'FAIL',conditionMessage(e))
    message('ERRO ',part,'/',run,': ',conditionMessage(e))
    NULL
  })
}

# ------------------------------------------------------------------------------
# 3. TESTE A: ALTERAR PROPORCAO DE NUCLEOS DA REFERENCIA, BULK FIXO
# ------------------------------------------------------------------------------
if ('A' %in% PARTS) {
  log_msg('\n=== 07d A: composicao da referencia | bulk=%s ===',A_DATASET)
  bulk <- load_bulk(A_DATASET)
  target0 <- as.numeric(table(factor(ref_meta$cell_type,levels=TYPES)))/nrow(ref_meta)
  names(target0) <- TYPES
  shift_to <- function(prop,type,target) {
    out <- prop
    change <- target-out[[type]]
    out[[type]] <- target
    out[['Ex']] <- out[['Ex']]-change
    if (any(out<=0) || abs(sum(out)-1)>1e-8) stop('Cenario A invalido')
    out
  }
  scenario_props <- list(
    ref_like=target0,
    in_low=shift_to(target0,'In',max(.025,target0[['In']]-.08)),
    in_high=shift_to(target0,'In',target0[['In']]+.08),
    opc_low=shift_to(target0,'Opc',max(.008,target0[['Opc']]-.025)),
    opc_high=shift_to(target0,'Opc',target0[['Opc']]+.06),
    vascular_high={
      v <- shift_to(target0,'End',.008)
      shift_to(v,'Per',.010)
    }
  )
  prop_targets <- rbindlist(lapply(names(scenario_props),function(nm) {
    p <- scenario_props[[nm]]
    data.table(scenario=nm,cell_type=TYPES,target_fraction=as.numeric(p))
  }))
  write_tab(prop_targets,file.path(OUT,'A','A_reference_targets.tsv'))
  ref_idx_by_type <- lapply(TYPES,function(ct) which(ref_meta$cell_type==ct))
  names(ref_idx_by_type) <- TYPES
  results_A <- list(); attr_A <- list()
  for (rep_id in seq_len(A_REPS)) {
    set.seed(SEED+rep_id)
    perm <- lapply(ref_idx_by_type,sample)
    for (scenario in names(scenario_props)) {
      counts_target <- as.integer(round(scenario_props[[scenario]]*A_N))
      names(counts_target) <- TYPES
      counts_target[['Ex']] <- A_N - sum(counts_target[TYPES!='Ex'])
      if (any(counts_target<1L) || any(counts_target>lengths(ref_idx_by_type))) {
        stop('A: impossivel amostrar sem reposicao em ',scenario,
             '. Reduza BISQUE07D_A_NUCLEI.')
      }
      selected <- sort(unlist(lapply(TYPES,function(ct) head(perm[[ct]],counts_target[[ct]])),
                              use.names=FALSE))
      m <- ref_meta[selected]
      x <- ref_all[,selected,drop=FALSE]
      ll <- ref_lib[selected]
      run <- paste0('A_',A_DATASET,'_',scenario,'_r',rep_id)
      folder <- file.path(OUT,'A',run)
      args <- list(ref_counts=x,ref_metadata=m,ref_libraries=ll,
                   bulk_counts=bulk$counts,bulk_libraries=bulk$lib,
                   marker_genes=rownames(bulk$counts),
                   context=paste0(run,'|',A_N,'|',SEED+rep_id,'|nuclear_sample_order_v1'))
      ans <- safe_call('A',run,function() run_fitted(run,args,folder))
      if (is.null(ans)) next
      p <- ans$proportions
      prior <- if(ncol(ans$sc_props)) rowMeans(ans$sc_props) else rep(NA_real_,8)
      results_A[[run]] <- data.table(
        dataset=A_DATASET,replicate=rep_id,scenario=scenario,cell_type=TYPES,
        n_reference_nuclei=length(selected),
        fraction_reference_pool=as.numeric(counts_target[TYPES]/A_N),
        fraction_reference_donor_mean=as.numeric(prior),
        mean_bulk_estimate=as.numeric(rowMeans(p)),
        sd_bulk_estimate=as.numeric(apply(p,1,sd)),
        genes_used=ans$genes_used)
      attr_A[[run]] <- data.table(replicate=rep_id,scenario=scenario,
                                 cell_type=TYPES,n_nuclei=as.integer(counts_target[TYPES]))
      log_msg('  %-18s r%d: In=%.2f%% Opc=%.2f%% End+Per=%.2f%%',scenario,rep_id,
              100*mean(p['In',]),100*mean(p['Opc',]),
              100*mean(colSums(p[c('End','Per'),,drop=FALSE])))
      rm(x,ans); invisible(gc(FALSE))
    }
  }
  if (length(results_A)) {
    aa <- rbindlist(results_A)
    aa_base <- aa[scenario=='ref_like',.(replicate,cell_type,
                                      baseline_mean=mean_bulk_estimate)]
    aa <- merge(aa,aa_base,by=c('replicate','cell_type'),all.x=TRUE,sort=FALSE)
    aa[,delta_from_ref_like_pp:=100*(mean_bulk_estimate-baseline_mean)]
    setorder(aa,replicate,scenario,cell_type)
    write_tab(aa,file.path(OUT,'A','07d_A_individual_fits.tsv'))
    write_tab(rbindlist(attr_A),file.path(OUT,'A','07d_A_reference_sample_sizes.tsv'))
    sum_a <- aa[,.(n_replicates=.N,
                    reference_fraction_mean=mean(fraction_reference_pool),
                    estimated_fraction_mean=mean(mean_bulk_estimate),
                    estimated_fraction_sd_across_reps=sd(mean_bulk_estimate),
                    delta_from_ref_like_mean_pp=mean(delta_from_ref_like_pp),
                    delta_from_ref_like_sd_pp=sd(delta_from_ref_like_pp)),
                by=.(dataset,scenario,cell_type)]
    write_tab(sum_a,file.path(OUT,'A','07d_A_composition_dependence.tsv'))
  }
}

# ------------------------------------------------------------------------------
# 4. TESTE B: MESMOS CONTROLES / CONTAGENS IDENTICAS / CONTEXTOS DIFERENTES
# ------------------------------------------------------------------------------
if ('B' %in% PARTS) {
  log_msg('\n=== 07d B: contexto bulk (44 controles identicos) ===')
  bh <- load_bulk('GSE64810')
  bp <- load_bulk('GSE68719')
  common_genes <- intersect(rownames(bh$counts),rownames(bp$counts))
  hd_ctrl <- as.character(bh$meta$sample_id[bh$meta$group=='Control'])
  pd_ctrl <- as.character(bp$meta$sample_id[bp$meta$group=='Control'])
  ctrl_ids <- intersect(hd_ctrl,pd_ctrl)
  if (length(ctrl_ids)!=44L) stop('B: esperados exatamente 44 controles HD/PD compartilhados')
  hd_case <- as.character(bh$meta$sample_id[bh$meta$group=='HD'])
  pd_case <- as.character(bp$meta$sample_id[bp$meta$group=='PD'])
  if (length(hd_case)!=20L || length(pd_case)!=29L ||
      anyDuplicated(c(ctrl_ids,hd_case,pd_case))) stop('B: casos/controles invalidos')
  # Bulk-control CANONICO: SEMPRE usa HD, tambem no contexto PD.
  # Dessa forma counts e library_size sao exatamente os mesmos nos 4 ajustes.
  x_ctrl <- bh$counts[common_genes,ctrl_ids,drop=FALSE]
  l_ctrl <- bh$lib[match(ctrl_ids,colnames(bh$counts))]
  x_hd <- bh$counts[common_genes,hd_case,drop=FALSE]
  l_hd <- bh$lib[match(hd_case,colnames(bh$counts))]
  x_pd <- bp$counts[common_genes,pd_case,drop=FALSE]
  l_pd <- bp$lib[match(pd_case,colnames(bp$counts))]
  contexts <- list(
    controls_only=list(x=x_ctrl,lib=l_ctrl),
    controls_plus_HD=list(x=cbind(x_ctrl,x_hd),lib=c(l_ctrl,l_hd)),
    controls_plus_PD=list(x=cbind(x_ctrl,x_pd),lib=c(l_ctrl,l_pd)),
    controls_plus_both=list(x=cbind(x_ctrl,x_hd,x_pd),lib=c(l_ctrl,l_hd,l_pd))
  )
  audit_b <- rbindlist(lapply(names(contexts),function(nm) {
    v <- contexts[[nm]]
    data.table(context=nm,n_samples=ncol(v$x),n_controls=length(ctrl_ids),
               n_genes=nrow(v$x),reference_source='Mathys_24')
  }))
  if (!all(vapply(contexts,function(v) identical(v$x[,ctrl_ids,drop=FALSE],x_ctrl) &&
                  identical(v$lib[match(ctrl_ids,colnames(v$x))],l_ctrl),logical(1)))) {
    stop('B: matriz canonica dos controles diverge entre contextos')
  }
  write_tab(audit_b,file.path(OUT,'B','07d_B_context_input_audit.tsv'))
  write_tab(data.table(sample_id=ctrl_ids,canonical_source='GSE64810',
                       library_size=l_ctrl),file.path(OUT,'B','07d_B_shared_controls.tsv'))
  Bprops <- list()
  B_summary <- list()
  for (nm in names(contexts)) {
    run <- paste0('B_',nm)
    v <- contexts[[nm]]
    folder <- file.path(OUT,'B',run)
    args <- list(ref_counts=ref_all,ref_metadata=ref_meta,
                 ref_libraries=ref_lib,bulk_counts=v$x,
                 bulk_libraries=v$lib,marker_genes=common_genes,
                 context=paste0('B|',nm,'|ctrl_HD_canonical|',length(common_genes)))
    ans <- safe_call('B',run,function() run_fitted(run,args,folder))
    if (is.null(ans)) next
    cc <- ans$proportions[,ctrl_ids,drop=FALSE]
    Bprops[[nm]] <- cc
    B_summary[[nm]] <- data.table(context=nm,cell_type=TYPES,
                                  mean_controls=as.numeric(rowMeans(cc)),
                                  sd_controls=as.numeric(apply(cc,1,sd)),
                                  n_controls=ncol(cc))
    rm(ans); invisible(gc(FALSE))
  }
  if (length(B_summary)) {
    write_tab(rbindlist(B_summary),file.path(OUT,'B','07d_B_control_means.tsv'))
    out_pair <- list()
    combos <- if(length(Bprops)>1) combn(names(Bprops),2,simplify=FALSE) else list()
    for (combo in combos) {
      pa <- Bprops[[combo[[1]]]]
      pb <- Bprops[[combo[[2]]]]
      dd <- 100*abs(pa-pb)
      out_pair[[length(out_pair)+1L]] <- data.table(
        context_a=combo[[1]],context_b=combo[[2]],cell_type=TYPES,
        mean_absolute_difference_pp=as.numeric(rowMeans(dd)),
        median_absolute_difference_pp=as.numeric(apply(dd,1,median)),
        max_absolute_difference_pp=as.numeric(apply(dd,1,max)),
        mean_signed_difference_pp=as.numeric(rowMeans(100*(pb-pa))))
    }
    if (length(out_pair)) {
      write_tab(rbindlist(out_pair),file.path(OUT,'B','07d_B_paired_context_differences.tsv'))
    }
    # Reproducao do 07 original HD, se checkpoints originais estao presentes.
    orig_path <- file.path(ORIG07_ROOT,'GSE64810','main24','top300','bisque_result.rds')
    if ('controls_plus_HD' %in% names(Bprops) && file.exists(orig_path)) {
      orig <- readRDS(orig_path)$bulk.props
      if (!is.null(orig) && all(TYPES %in% rownames(orig)) &&
          all(ctrl_ids %in% colnames(orig))) {
        dif <- max(abs(orig[TYPES,ctrl_ids,drop=FALSE]-
                       Bprops[['controls_plus_HD']][TYPES,ctrl_ids,drop=FALSE]))
        write_tab(data.table(comparison='B_HD_vs_07_original',
                             max_abs_diff=dif,reproduced=(dif<1e-6)),
                  file.path(OUT,'B','07d_B_reproduction_07.tsv'))
      }
    }
  }
}

# ------------------------------------------------------------------------------
# 5. TESTE C: VALIDACAO SINTETICA LODO COM MARCADORES TRAIN-ONLY
# ------------------------------------------------------------------------------
COMPOSITIONS <- rbind(
  equal=c(Ex=.125,In=.125,Ast=.125,Oli=.125,Opc=.125,Mic=.125,End=.125,Per=.125),
  cortical_like=c(Ex=.38,In=.17,Ast=.12,Oli=.20,Opc=.05,Mic=.04,End=.02,Per=.02),
  inhibitory_enriched=c(Ex=.22,In=.30,Ast=.12,Oli=.20,Opc=.06,Mic=.05,End=.03,Per=.02),
  opc_enriched=c(Ex=.20,In=.12,Ast=.13,Oli=.22,Opc=.20,Mic=.07,End=.03,Per=.03),
  vascular_enriched=c(Ex=.23,In=.12,Ast=.15,Oli=.18,Opc=.06,Mic=.06,End=.10,Per=.10),
  low_in=c(Ex=.45,In=.01,Ast=.13,Oli=.25,Opc=.06,Mic=.06,End=.02,Per=.02),
  low_opc=c(Ex=.38,In=.16,Ast=.14,Oli=.22,Opc=.01,Mic=.05,End=.02,Per=.02),
  low_vascular=c(Ex=.43,In=.17,Ast=.14,Oli=.19,Opc=.04,Mic=.02,End=.005,Per=.005),
  low_all_rare=c(Ex=.51,In=.01,Ast=.16,Oli=.28,Opc=.01,Mic=.02,End=.005,Per=.005)
)
if (!identical(colnames(COMPOSITIONS),TYPES) ||
    any(COMPOSITIONS<=0) ||
    max(abs(rowSums(COMPOSITIONS)-1))>1e-10) stop('Composicoes C invalidas')
VIRTUAL_N <- 1000
VCOUNTS <- COMPOSITIONS * VIRTUAL_N
if (any(abs(VCOUNTS-round(VCOUNTS))>1e-8)) stop('Nucleos virtuais nao inteiros')

if ('C' %in% PARTS) {
  log_msg('\n=== 07d C: validacao sintetica LODO (%d doadores) ===',C_HOLDOUTS)
  support <- ref_meta[,.N,by=.(donor_id,cell_type)]
  support_wide <- dcast(support,donor_id~cell_type,value.var='N',fill=0L)
  eligible <- support_wide[
    rowSums(as.matrix(support_wide[,..TYPES])>0L)==length(TYPES)]
  if (nrow(eligible)<C_HOLDOUTS) stop('Menos doadores elegiveis que solicitado')
  # Mesma semente e regra da auditoria 08d, para comparabilidade exploratoria.
  set.seed(62026L)
  chosen <- eligible[sample(seq_len(nrow(eligible)),C_HOLDOUTS,replace=FALSE)]
  chosen[,holdout:=sprintf('holdout_%02d',seq_len(.N))]
  write_tab(chosen[,c('holdout',TYPES),with=FALSE],
            file.path(OUT,'C','07d_C_holdout_support.tsv'))

  # Pseudobulks por doador/tipo de TODOS OS GENES exclusivamente para
  # RESELECAO dos marcadores nos doadores restantes em cada fold.
  groups <- unique(ref_meta[,.(donor_id,cell_type)])
  setorder(groups,donor_id,cell_type)
  gid <- match(paste(ref_meta$donor_id,ref_meta$cell_type,sep='\r'),
               paste(groups$donor_id,groups$cell_type,sep='\r'))
  if (anyNA(gid)) stop('Grupos doador/tipo desalinhados')
  G <- Matrix::sparseMatrix(i=seq_len(ncol(ref_all)),j=gid,x=rep(1,ncol(ref_all)),
                             dims=c(ncol(ref_all),nrow(groups)))
  pb_counts <- ref_all %*% G
  pb_lib <- as.numeric(Matrix::colSums(pb_counts))
  if (any(pb_lib<=0)) stop('Pseudobulk de treino com biblioteca vazia')
  rm(G); invisible(gc(FALSE))

  train_only_markers <- function(donor) {
    kk <- which(groups$donor_id!=donor)
    pb <- as.matrix(pb_counts[,kk,drop=FALSE]); storage.mode(pb) <- 'double'
    pcpm <- sweep(pb,2,pb_lib[kk],'/')*1e6
    avg <- matrix(0,nrow(pcpm),length(TYPES),
                  dimnames=list(rownames(ref_all),TYPES))
    det <- avg
    for (ct in TYPES) {
      ii <- which(groups$cell_type[kk]==ct)
      if (!length(ii)) stop('Fold sem pseudobulk ',ct)
      avg[,ct] <- rowMeans(pcpm[,ii,drop=FALSE])
      det[,ct] <- rowMeans(pcpm[,ii,drop=FALSE]>1)
    }
    den <- pmax(rowSums(avg),1e-12)
    chosen_genes <- list()
    for (ct in TYPES) {
      sp <- avg[,ct]/den
      score <- sp*log1p(avg[,ct])*det[,ct]
      valid <- which(is.finite(score) & score>0)
      if (length(valid)<300L) stop('Menos de 300 candidatos treinados: ',ct)
      od <- order(-score[valid],-sp[valid],-avg[valid,ct],rownames(ref_all)[valid])
      chosen_genes[[ct]] <- rownames(ref_all)[valid[od[seq_len(300L)]]]
    }
    unique(unlist(chosen_genes,use.names=FALSE))
  }

  # Geracao de mistura in silico a partir da EXPRESSAO MEDIA POR NUCLEO
  # de cada tipo no doador de teste. O RNA completo do nucleo e incluido.
  make_mixture <- function(donor) {
    mm <- matrix(0,nrow(ref_all),length(TYPES),
                 dimnames=list(rownames(ref_all),TYPES))
    for (ct in TYPES) {
      jj <- which(ref_meta$donor_id==donor & ref_meta$cell_type==ct)
      if (!length(jj)) stop('Doador de teste sem tipo ',ct)
      mm[,ct] <- as.numeric(Matrix::rowMeans(ref_all[,jj,drop=FALSE]))
    }
    mix <- round(mm %*% t(VCOUNTS))
    colnames(mix) <- paste0('synthetic_',rownames(COMPOSITIONS))
    rownames(mix) <- rownames(ref_all)
    lib <- as.numeric(colSums(mix))
    if (any(lib<=0)) stop('Mistura simulada sem RNA')
    # Contribuicao de RNA baseada no total esperado de reads por tipo.
    size_per_type <- colSums(mm)
    weighted <- sweep(t(VCOUNTS),1,size_per_type,'*')
    rna_prop <- sweep(weighted,2,colSums(weighted),'/')
    list(counts=mix,lib=lib,
         nuclei=t(COMPOSITIONS),rna=rna_prop)
  }

  Cpred <- list(); Cmetrics <- list(); Csupport <- list()
  for (fold in seq_len(nrow(chosen))) {
    donor <- as.character(chosen$donor_id[[fold]])
    label <- as.character(chosen$holdout[[fold]])
    train_idx <- which(ref_meta$donor_id!=donor)
    if (any(ref_meta$donor_id[train_idx]==donor)) stop('LODO: vazamento de doador')
    markers <- train_only_markers(donor)
    write_tab(data.table(holdout=label,gene_id=markers),
              file.path(OUT,'C',paste0(label,'_train_markers.tsv')))
    mix <- make_mixture(donor)
    run <- paste0('C_',label)
    folder <- file.path(OUT,'C',run)
    args <- list(ref_counts=ref_all[,train_idx,drop=FALSE],
                 ref_metadata=ref_meta[train_idx],
                 ref_libraries=ref_lib[train_idx],
                 bulk_counts=mix$counts,bulk_libraries=mix$lib,
                 marker_genes=markers,
                 context=paste0(run,'|',length(markers),'|train_only|',paste(rownames(COMPOSITIONS),collapse=',')))
    ans <- safe_call('C',run,function() run_fitted(run,args,folder))
    if (is.null(ans)) next
    p <- ans$proportions
    if (!identical(colnames(p),colnames(mix$counts))) stop('C: sinteticos desalinhados')
    long <- rbindlist(lapply(seq_len(ncol(p)),function(i) {
      data.table(holdout=label,scenario=rownames(COMPOSITIONS)[[i]],
                 cell_type=TYPES,truth_nucleus=as.numeric(mix$nuclei[,i]),
                 truth_rna=as.numeric(mix$rna[,i]),estimated=as.numeric(p[,i]))
    }))
    # data.table avalia o lado direito de := antes de criar as colunas.
    # Portanto, calcular os erros primeiro, em uma chamada separada,
    # e os valores absolutos somente depois (corrige object not found).
    long[,`:=`(error_nucleus_pp=100*(estimated-truth_nucleus),
               error_rna_pp=100*(estimated-truth_rna))]
    long[,`:=`(abs_error_nucleus_pp=abs(error_nucleus_pp),
               abs_error_rna_pp=abs(error_rna_pp))]
    Cpred[[label]] <- long
    Cmetrics[[label]] <- data.table(
      holdout=label,n_scenarios=ncol(p),n_predictions=nrow(long),
      mae_nucleus_pp=mean(long$abs_error_nucleus_pp),
      rmse_nucleus_pp=sqrt(mean(long$error_nucleus_pp^2)),
      mae_rna_pp=mean(long$abs_error_rna_pp),
      rmse_rna_pp=sqrt(mean(long$error_rna_pp^2)),
      n_genes=ans$genes_used,
      max_sample_sum_error=max(abs(colSums(p)-1)))
    Csupport[[label]] <- data.table(holdout=label,
                                   train_donors=uniqueN(ref_meta$donor_id[train_idx]),
                                   train_nuclei=length(train_idx),
                                   test_nuclei=sum(ref_meta$donor_id==donor))
    log_msg('  %-12s MAE nucleo %.2f pp | MAE RNA %.2f pp',label,
            Cmetrics[[label]]$mae_nucleus_pp,Cmetrics[[label]]$mae_rna_pp)
    rm(args,ans,mix); invisible(gc(FALSE))
  }
  if (length(Cpred)) {
    pp <- rbindlist(Cpred)
    write_tab(pp,file.path(OUT,'C','07d_C_predictions.tsv'))
    write_tab(rbindlist(Cmetrics),file.path(OUT,'C','07d_C_metrics_by_holdout.tsv'))
    write_tab(rbindlist(Csupport),file.path(OUT,'C','07d_C_support.tsv'))
    bytype <- pp[,.(n_predictions=.N,
                    mae_nucleus_pp=mean(abs_error_nucleus_pp),
                    bias_nucleus_pp=mean(error_nucleus_pp),
                    mae_rna_pp=mean(abs_error_rna_pp),
                    bias_rna_pp=mean(error_rna_pp),
                    fraction_zero_for_positive_truth=mean(estimated<=1e-10 & truth_nucleus>0)),
                 by=cell_type]
    write_tab(bytype,file.path(OUT,'C','07d_C_metrics_by_type.tsv'))
    byscenario <- pp[,.(n_predictions=.N,
                        mae_nucleus_pp=mean(abs_error_nucleus_pp),
                        mae_rna_pp=mean(abs_error_rna_pp)),by=scenario]
    write_tab(byscenario,file.path(OUT,'C','07d_C_metrics_by_scenario.tsv'))
    global <- data.table(n_holdouts=uniqueN(pp$holdout),n_scenarios=uniqueN(pp$scenario),
                         n_predictions=nrow(pp),
                         mae_nucleus_pp=mean(pp$abs_error_nucleus_pp),
                         rmse_nucleus_pp=sqrt(mean(pp$error_nucleus_pp^2)),
                         mae_rna_pp=mean(pp$abs_error_rna_pp),
                         rmse_rna_pp=sqrt(mean(pp$error_rna_pp^2)))
    write_tab(global,file.path(OUT,'C','07d_C_metrics_global.tsv'))
    log_msg(' C global: MAE nucleos %.3f pp; MAE RNA %.3f pp',
            global$mae_nucleus_pp,global$mae_rna_pp)
  }
}

# ------------------------------------------------------------------------------
# 6. VISUALIZACOES DESCRITIVAS DOS EXPERIMENTOS (sem recalc. de proporcoes)
# ------------------------------------------------------------------------------
# Usar arquivos TSV como fonte. Erro em figura nao invalida estimativas salvas.
if (requireNamespace('ggplot2',quietly=TRUE)) {
  safe_plot <- function(name,fun) {
    tryCatch({
      p <- fun()
      ggplot2::ggsave(file.path(OUT,'plots',name),plot=p,
                      width=12,height=7,dpi=180,bg='white',limitsize=FALSE)
      record_status('PLOT',name,'PASS')
    },error=function(e) {
      record_status('PLOT',name,'WARN',conditionMessage(e))
    })
  }
  theme_clean <- ggplot2::theme_minimal(base_size=12) +
    ggplot2::theme(plot.title=ggplot2::element_text(face='bold'),
                   axis.text.x=ggplot2::element_text(angle=30,hjust=1),
                   plot.margin=ggplot2::margin(12,18,18,12))
  ap <- file.path(OUT,'A','07d_A_individual_fits.tsv')
  if ('A' %in% PARTS && file.exists(ap)) safe_plot('07d_A_ancoragem_referencia.png',function() {
    x <- fread(ap)[cell_type %in% c('In','Opc','End','Per')]
    x[,`:=`(reference_percent=100*fraction_reference_pool,
             bulk_percent=100*mean_bulk_estimate)]
    ggplot2::ggplot(x,ggplot2::aes(x=reference_percent,y=bulk_percent,
                                   color=scenario,shape=factor(replicate))) +
      ggplot2::geom_point(size=3,alpha=.8) +
      ggplot2::facet_wrap(~cell_type,scales='free',ncol=2) +
      ggplot2::labs(title='Bisque: dependencia da composicao da referencia',
                    subtitle='Mesmo bulk, composicao nuclear da referencia alterada; cada ponto e uma replica',
                    x='Frequencia de nucleos na referencia (%)',
                    y='Proporcao media estimada no bulk (%)',shape='Replica',color='Cenario')+
      theme_clean
  })
  bp <- file.path(OUT,'B','07d_B_paired_context_differences.tsv')
  if ('B' %in% PARTS && file.exists(bp)) safe_plot('07d_B_contexto_bulk.png',function() {
    x <- fread(bp)
    x[,comparison:=paste(context_a,'vs',context_b)]
    ggplot2::ggplot(x,ggplot2::aes(x=cell_type,y=mean_absolute_difference_pp,
                                   fill=comparison))+
      ggplot2::geom_col(position=ggplot2::position_dodge(width=.85),width=.8)+
      ggplot2::labs(title='Bisque: o mesmo controle em diferentes contextos bulk',
                    subtitle='Mesmos counts e tamanhos de biblioteca nos controles',
                    x='Tipo celular',y='Diferenca absoluta media (pontos percentuais)',
                    fill='Comparacao')+
      theme_clean
  })
  cp <- file.path(OUT,'C','07d_C_metrics_by_type.tsv')
  if ('C' %in% PARTS && file.exists(cp)) safe_plot('07d_C_erro_por_tipo.png',function() {
    x <- fread(cp)
    xl <- data.table::melt(x,id.vars='cell_type',
              measure.vars=c('mae_nucleus_pp','mae_rna_pp'),
              variable.name='estimando',value.name='mae_pp')
    xl[,estimando:=factor(estimando,
             levels=c('mae_nucleus_pp','mae_rna_pp'),
             labels=c('Fracao de nucleos','Contribuicao de RNA'))]
    ggplot2::ggplot(xl,ggplot2::aes(x=cell_type,y=mae_pp,fill=estimando))+
      ggplot2::geom_col(position='dodge',width=.7)+
      ggplot2::labs(title='Bisque: validacao sintetica, erro por tipo',
                    subtitle='Quatro doadores excluidos individualmente; estimandos diferentes',
                    x='Tipo celular',y='Erro absoluto medio (pp)',fill='Verdade')+
      theme_clean
  })
  predp <- file.path(OUT,'C','07d_C_predictions.tsv')
  if ('C' %in% PARTS && file.exists(predp)) safe_plot('07d_C_verdade_vs_estimativa_In_Opc.png',function() {
    x <- fread(predp)[cell_type %in% c('In','Opc')]
    x[,`:=`(truth_pp=100*truth_nucleus,estimate_pp=100*estimated)]
    ggplot2::ggplot(x,ggplot2::aes(x=truth_pp,y=estimate_pp,color=holdout))+
      ggplot2::geom_abline(slope=1,intercept=0,linetype=2,color='gray45')+
      ggplot2::geom_point(size=2,alpha=.75)+
      ggplot2::facet_wrap(~cell_type,scales='free')+
      ggplot2::labs(title='Bisque: recuperacao de In e Opc em misturas sinteticas',
                    subtitle='Linha pontilhada indica estimativa igual ao valor conhecido',
                    x='Fracao nuclear definida na mistura (%)',
                    y='Proporcao estimada (%)',color='Doador excluido')+
      theme_clean
  })
} else {
  record_status('PLOT','ggplot2','WARN','ggplot2 nao instalado; TSVs preservados')
}

# ------------------------------------------------------------------------------
# 7. FINALIZACAO, INTEGRIDADE E REGISTRO DE VERSOES
# ------------------------------------------------------------------------------
if (length(statuses)) {
  ss <- rbindlist(statuses,fill=TRUE)
  write_tab(ss,file.path(OUT,'07d_status.tsv'))
  log_msg('\n=== 07d STATUS ===')
  print(ss[,.(n_runs=.N),by=.(part,status)])
} else {
  stop('Nenhuma auditoria foi executada')
}
writeLines(capture.output(sessionInfo()),file.path(OUT,'07d_sessionInfo.txt'))
write_tab(data.table(
  property=c('parts','A_dataset','A_nuclei','A_reps','C_holdouts','reuse',
             'bisque_version','old_cpm','use_overlap','reference_nuclei','reference_donors'),
  value=c(paste(PARTS,collapse=','),A_DATASET,as.character(A_N),as.character(A_REPS),
          as.character(C_HOLDOUTS),as.character(REUSE),
          as.character(packageVersion('BisqueRNA')),'FALSE','FALSE',
          as.character(ncol(ref_all)),as.character(uniqueN(ref_meta$donor_id)))),
  file.path(OUT,'07d_run_parameters.tsv'))
log_msg('Arquivos em: %s',OUT)
if (any(vapply(statuses,function(s) any(s$status=='FAIL'),logical(1)))) {
  stop('07d terminou com falhas; verificar 07d_status.tsv (outputs anteriores preservados)')
}
log_msg('07d concluido sem falhas computacionais. Ainda NAO constitui validacao externa.')
