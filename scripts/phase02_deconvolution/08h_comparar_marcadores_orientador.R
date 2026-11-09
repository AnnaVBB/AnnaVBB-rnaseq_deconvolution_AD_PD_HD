#!/usr/bin/env Rscript
# =============================================================================
# 08h v3 | Comparar os marcadores do orientador com nosso script 04 e dtangle 08e/08g
# Projeto: deconvolucao AD-PD-HD | Mathys 2019
#
# INSTALACAO: copiar os arquivos deste pacote para a raiz do repositorio.
# EXECUCAO:  Rscript scripts/phase02_deconvolution/08h_comparar_marcadores_orientador.R
# O script NAO modifica 04, 08e, 08g, 09, 11, nem proporcoes do bMIND.
#
# Comparacoes:
#  A) Professor dtangle_top100_automatico vs nossos 08e (2 referencias, por coorte)
#     e 08g (pseudobulk por doador, por coorte).
#  B) Professor top300_especificidade vs nosso 04 top300 por marker_score.
#  C) Opcional: professor top300 vs ranking local alternativo por especificidade
#     em CPM *linear*, gerado apenas para diagnostico. Nao equivale ao ranking
#     do professor que usa log2CPM nem corrige diferencas de universo genico.
#  D) Comparacao interna entre as duas listas do professor (criterios distintos).
#
# LIMITACAO FUNDAMENTAL: as listas fornecidas NAO incluem a totalidade do
# universo de ~9300 genes testaveis do orientador. Portanto, nao e possivel
# reproduzir o mesmo universo genico a partir desta planilha isoladamente.
# Jaccard considera as LISTAS observadas, e Spearman so os IDs compartilhados.
# =============================================================================
options(stringsAsFactors = FALSE)
if (!requireNamespace('data.table', quietly=TRUE)) stop('Instale data.table no ambiente R atual')
suppressPackageStartupMessages(library(data.table))
ROOT <- normalizePath(Sys.getenv('DECONV_ROOT', '.'), mustWork = TRUE)
P <- function(...) file.path(ROOT, ...)
CELLS <- c('Ex','In','Ast','Oli','Opc','Mic','End','Per')
TEACHER <- P('data','reference','mathys2019','orientador_markers')
OUT <- P('results','deconvolution','dtangle','08h_comparacao_marcadores_orientador')
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
W <- function(x, fname) { fwrite(x, file.path(OUT, fname), sep='\t', na='NA', quote=FALSE); invisible(x) }
STATUS <- list(); ADD_STATUS <- function(stage,status,detail) {
  STATUS[[length(STATUS)+1L]] <<- data.table(stage=stage,status=status,detail=detail)
  cat(sprintf('08h | %s | %s | %s\n',stage,status,detail))
}

# Todos os identificadores sao Entrez; comparar como texto, nunca por symbol.
check_ranked <- function(x,label,k,group_cols='cell_type') {
  need <- c(group_cols,'cell_type','gene_id','rank')
  miss <- setdiff(need,names(x))
  if (length(miss)) stop(label,': colunas ausentes: ',paste(miss,collapse=','))
  x[, `:=`(gene_id=trimws(as.character(gene_id)),cell_type=trimws(as.character(cell_type)),rank=as.integer(rank))]
  if (anyNA(x$gene_id)||any(!grepl('^[0-9]+$',x$gene_id))||anyNA(x$rank)||
      any(!x$cell_type %in% CELLS)) stop(label,': Entrez, ranking ou tipo celular invalido')
  if (anyDuplicated(x,by=c(group_cols,'gene_id')) || anyDuplicated(x,by=c(group_cols,'rank'))) {
    stop(label,': GeneID ou rank duplicado em um mesmo tipo/grupo')
  }
  # Conferir numeros e cobertura dos ranks, sem depender de identical(),
  # que exige tambem identidade estrita de atributos e tipos R.
  tests <- x[, .(
    n = .N,
    n_ranks = uniqueN(rank),
    min_rank = min(rank),
    max_rank = max(rank),
    faltantes = paste(setdiff(seq_len(k), rank), collapse = ','),
    fora_do_intervalo = paste(setdiff(rank, seq_len(k)), collapse = ',')
  ), by = group_cols]
  invalid <- tests[
    n != k | n_ranks != k | min_rank != 1L | max_rank != k |
      nzchar(faltantes) | nzchar(fora_do_intervalo)
  ]
  if (nrow(invalid)) {
    stop(label, ': rankings invalidos (esperado 1..', k, ' por grupo).\n',
         paste(capture.output(print(invalid)), collapse='\n'))
  }
  x
}

read_teacher <- function(fname,k,kind) {
  path <- file.path(TEACHER,fname)
  if (!file.exists(path)) stop('Arquivo do orientador ausente: ',path)
  x <- fread(path, colClasses=list(character=c('entrez','symbol')))
  if (!all(c('cell_type','rank','entrez','symbol') %in% names(x))) stop('Colunas invalidas em: ',path)
  x[,`:=`(gene_id=as.character(entrez), rank=as.integer(rank),symbol=as.character(symbol))]
  x <- x[,.(cell_type,rank,gene_id,symbol)]
  x <- check_ranked(x,paste0('orientador/',kind),k)
  if (!setequal(unique(x$cell_type),CELLS) || nrow(x)!=8L*k) stop(kind,': esperados 8 tipos x ',k)
  x[,source:=kind]
  x
}

teacher100 <- read_teacher('dtangle_top100_automatico.tsv',100L,'orientador_dtangle_top100')
teacher300 <- read_teacher('top300_especificidade.tsv',300L,'orientador_especificidade_top300')
ADD_STATUS('PLANILHA','PASS',sprintf('top100: %d atribuicoes / %d IDs unicos; top300: %d atribuicoes / %d IDs unicos',
  nrow(teacher100),uniqueN(teacher100$gene_id),nrow(teacher300),uniqueN(teacher300$gene_id)))

# Confirma que a planilha nao associa simbolos conflitantes a um mesmo Entrez.
symbol_map <- unique(rbindlist(list(
  teacher100[,.(gene_id,symbol)],teacher300[,.(gene_id,symbol)]
)))
sym_conflict <- symbol_map[,.(n_symbols=uniqueN(symbol), symbols=paste(unique(symbol),collapse=';')),by=gene_id][n_symbols>1L]
W(sym_conflict,'08h_conflicts_entrez_symbol.tsv')
if (nrow(sym_conflict) > 0L) {
  ADD_STATUS('ANOTACAO','WARN',paste(nrow(sym_conflict),'Entrez com simbolos conflitantes'))
} else {
  ADD_STATUS('ANOTACAO','PASS','Sem simbolos conflitantes dentro da planilha')
}

# Auditar os Entrez GFAP=2670 e ADGRV1=84059 (na planilha do professor).
spot <- rbindlist(list(teacher100,teacher300),use.names=TRUE,fill=TRUE)[gene_id %in% c('2670','84059')]
W(spot[,.(source,cell_type,rank,gene_id,symbol)],'08h_spotcheck_GFAP_ADGRV1.tsv')

# Calcular comparacao de duas listas de ranked Entrez sem misturar tipos.
compare_ranked <- function(a,b,source_label,local_label,dataset=NA_character_,reference_mode=NA_character_) {
  metric_rows <- list(); gene_rows <- list()
  for (ct in CELLS) {
    a1 <- a[cell_type==ct,.(gene_id, rank_advisor=rank, symbol)]
    b1 <- b[cell_type==ct,.(gene_id, rank_local=rank)]
    if (!nrow(a1)||!nrow(b1)) next
    aa <- a1$gene_id; bb <- b1$gene_id
    overlap <- intersect(aa,bb); unions <- union(aa,bb)
    rho <- if (length(overlap)>=3L) suppressWarnings(
      cor(a1$rank_advisor[match(overlap,aa)],b1$rank_local[match(overlap,bb)],method='spearman')) else NA_real_
    # A correlacao so considera genes compartilhados: nao e uma medida de acuracia.
    metric_rows[[length(metric_rows)+1L]] <- data.table(
      comparison=source_label, local_source=local_label, dataset=dataset, reference_mode=reference_mode,
      cell_type=ct,n_advisor=length(aa),n_local=length(bb),n_common=length(overlap),
      n_union=length(unions),jaccard=length(overlap)/length(unions),
      spearman_rank_overlap=as.numeric(rho),
      n_top10_teacher_found=sum(a1[rank_advisor<=10,gene_id] %in% bb),
      n_teacher_missing=length(setdiff(aa,bb)), n_local_only=length(setdiff(bb,aa)))
    z <- merge(a1,b1,by='gene_id',all=TRUE,sort=FALSE)
    z[,`:=`(comparison=source_label,local_source=local_label,dataset=dataset,
            reference_mode=reference_mode,cell_type=ct,
            presence=fifelse(!is.na(rank_advisor)&!is.na(rank_local),'BOTH',
                    fifelse(!is.na(rank_advisor),'ADVISOR_ONLY','LOCAL_ONLY')))]
    gene_rows[[length(gene_rows)+1L]] <- z[,.(comparison,local_source,dataset,reference_mode,cell_type,
      gene_id,symbol,rank_advisor,rank_local,presence)]
  }
  list(summary=rbindlist(metric_rows,fill=TRUE), genes=rbindlist(gene_rows,fill=TRUE))
}
RESULT <- list(); DETAILS <- list()
ADD <- function(ret) {
  if (nrow(ret$summary)) RESULT[[length(RESULT)+1L]] <<- ret$summary
  if (nrow(ret$genes)) DETAILS[[length(DETAILS)+1L]] <<- ret$genes
}

# A) Comparacao entre listas internas do orientador (metodos diferentes).
ADD(compare_ranked(teacher100,teacher300,'orientador_dtangle_top100_vs_especificidade_top300',
                   'orientador_especificidade_top300'))

# B) Comparacao com 04, usando o conjunto PRINCIPAL real do nosso script.
our04 <- P('data','reference','mathys2019','phase02_prepared','markers','mathys2019_top300_markers.tsv')
if (file.exists(our04)) {
  a <- fread(our04)
  if (!all(c('cell_type','gene_id','marker_rank') %in% names(a))) stop('04 top300: colunas inesperadas')
  a <- a[,.(cell_type, gene_id=as.character(gene_id),rank=as.integer(marker_rank))]
  # Alguns runs do script 04 exportam apenas top100 apesar do nome do arquivo.
  # Nao confundir essa lista com top300 nem interromper a auditoria dtangle.
  n_by_cell <- a[,.(n=.N),by=cell_type]
  if (!setequal(n_by_cell$cell_type,CELLS)) {
    ADD_STATUS('04_TOP300','WARN','Arquivo 04 sem os oito tipos; comparacao 04 ignorada')
  } else if (all(n_by_cell$n == 300L)) {
    a <- check_ranked(a,'04_marker_score_top300',300L)
    ADD(compare_ranked(teacher300,a,'orientador_especificidade_top300_vs_04_marker_score_top300',
                       '04_marker_score_top300'))
    ADD_STATUS('04_TOP300','PASS','Comparacao contra top300 principal, ordenado por marker_score')
  } else if (all(n_by_cell$n == 100L)) {
    a <- check_ranked(a,'04_marker_score_top100',100L)
    ADD(compare_ranked(teacher300[rank<=100L],a,
                       'orientador_especificidade_top100_de_300_vs_04_marker_score_top100',
                       '04_marker_score_top100'))
    ADD_STATUS('04_TOP300','WARN','Arquivo 04 contem apenas top100 por tipo; comparacao exploratoria top100 feita; top300 completo pendente')
  } else {
    ADD_STATUS('04_TOP300','WARN',paste('Contagens por tipo diferentes de 100/300; comparacao 04 ignorada:',
               paste(paste(n_by_cell$cell_type,n_by_cell$n,sep='='),collapse='; ')))
  }
  if (all(n_by_cell$n %in% c(100L,300L)) &&
      (all(n_by_cell$n == 100L) || all(n_by_cell$n == 300L))) {
    ast_local <- a[cell_type=='Ast' & rank<=10]
    ast_prof100 <- teacher100[cell_type=='Ast' & rank<=10,.(cell_type,gene_id,rank)]
    ast_prof300 <- teacher300[cell_type=='Ast' & rank<=10,.(cell_type,gene_id,rank)]
    W(rbindlist(list(
      ast_local[,.(source='04_marker_score',cell_type,rank,gene_id)],
      ast_prof100[,.(source='orientador_dtangle_top100',cell_type,rank,gene_id)],
      ast_prof300[,.(source='orientador_especificidade_top300',cell_type,rank,gene_id)]
    )), '08h_top10_Ast_orientador_e_04.tsv')
  }
} else {
  ADD_STATUS('04_TOP300','PENDING','Falta: data/reference/mathys2019/phase02_prepared/markers/mathys2019_top300_markers.tsv')
}

# C) Sensibilidade metodologica opcional: especificidade local (CPM linear).
metricfile <- P('data','reference','mathys2019','phase02_prepared','markers','mathys2019_marker_metrics_long.tsv')
if (file.exists(metricfile)) {
  m <- fread(metricfile)
  if (all(c('cell_type','gene_id','specificity','mean_CPM') %in% names(m))) {
    m <- m[cell_type %in% CELLS & is.finite(specificity) & specificity>0 & is.finite(mean_CPM) & mean_CPM>0]
    m[,gene_id:=as.character(gene_id)]
    setorder(m,cell_type,-specificity,-mean_CPM,gene_id)
    m[,rank:=seq_len(.N),by=cell_type]
    m <- m[rank<=300,.(cell_type,gene_id,rank)]
    if (nrow(m)==2400L) {
      m <- check_ranked(m,'04_specificity_linear_cpm_top300',300L)
      ADD(compare_ranked(teacher300,m,'orientador_especificidade_log2CPM_vs_04_especificidade_CPM_linear',
                         '04_specificity_linear_CPM_exploratory'))
      ADD_STATUS('04_ESPECIFICIDADE','PASS','Sensibilidade exploratoria; metricas de especificidade em escalas DIFERENTES')
    } else {
      ADD_STATUS('04_ESPECIFICIDADE','SKIP','Menos de 300 candidatos em um ou mais tipos')
    }
  } else {
    ADD_STATUS('04_ESPECIFICIDADE','SKIP','Colunas de especificidade/CPM ausentes')
  }
} else {
  ADD_STATUS('04_ESPECIFICIDADE','SKIP','Arquivo de metricas do 04 ausente; sensibilidade opcional')
}

# D) Comparar dtangle automatico: cada coorte e configuracao separadamente.
read_local_dtangle <- function(path,label) {
  if (!file.exists(path)) { ADD_STATUS(label,'PENDING',paste('Arquivo ausente:',path)); return(NULL) }
  x <- fread(path)
  if (!all(c('dataset','reference_mode','cell_type','rank','gene_id') %in% names(x))) stop(label,': colunas inesperadas')
  x[,`:=`(dataset=as.character(dataset),reference_mode=as.character(reference_mode),
          cell_type=as.character(cell_type),gene_id=as.character(gene_id),rank=as.integer(rank))]
  x <- check_ranked(x,label,100L,group_cols=c('dataset','reference_mode','cell_type'))
  ADD_STATUS(label,'PASS',paste('Configuracoes encontradas:',uniqueN(x,by=c('dataset','reference_mode'))))
  x
}
for (item in list(
  list('08e',P('results','deconvolution','dtangle','08e_replicacao_orientador','08e_marker_assignments.tsv')),
  list('08g',P('results','deconvolution','dtangle','08g_media_por_doador','08g_marker_assignments.tsv'))
)) {
  label <- item[[1]]; local <- read_local_dtangle(item[[2]],label)
  if (is.null(local)) next
  configs <- unique(local[,.(dataset,reference_mode)])
  for (i in seq_len(nrow(configs))) {
    ds <- configs$dataset[i]; mode <- configs$reference_mode[i]
    b <- local[dataset==ds & reference_mode==mode,.(cell_type,gene_id,rank)]
    ADD(compare_ranked(teacher100,b,'orientador_dtangle_top100_vs_local_dtangle',
                       label,dataset=ds,reference_mode=mode))
  }
}

result <- rbindlist(RESULT,use.names=TRUE,fill=TRUE)
detail <- rbindlist(DETAILS,use.names=TRUE,fill=TRUE)
setorder(result,comparison,local_source,dataset,reference_mode,cell_type)
W(result,'08h_comparacao_jaccard_spearman.tsv')
W(detail,'08h_marcadores_gene_a_gene.tsv')
W(result[comparison=='orientador_dtangle_top100_vs_local_dtangle'], '08h_dtangle_top100_vs_08e_08g.tsv')
W(result[comparison %in% c('orientador_especificidade_top300_vs_04_marker_score_top300',
                            'orientador_especificidade_log2CPM_vs_04_especificidade_CPM_linear')],
  '08h_top300_vs_04.tsv')
W(result[comparison=='orientador_dtangle_top100_vs_especificidade_top300'],
  '08h_comparacao_interna_orientador.tsv')
W(result[comparison=='orientador_especificidade_top100_de_300_vs_04_marker_score_top100'],
  '08h_top100_especificidade_vs_04_top100_EXPLORATORIO.tsv')
W(rbindlist(STATUS,fill=TRUE),'08h_status.tsv')

# Gráfico de Jaccard, se o ambiente possuir ggplot2; ausencia nao impede analise.
if (requireNamespace('ggplot2',quietly=TRUE) && nrow(result)>0) {
  suppressPackageStartupMessages(library(ggplot2))
  graph_dt <- result[comparison=='orientador_dtangle_top100_vs_local_dtangle']
  if (nrow(graph_dt)) {
    graph_dt[,config:=paste(local_source,dataset,reference_mode,sep=' / ')]
    g <- ggplot(graph_dt,aes(x=cell_type,y=jaccard,fill=config))+
      geom_col(position='dodge')+coord_cartesian(ylim=c(0,1))+
      labs(title='08h | Jaccard: top100 do orientador vs dtangle local',
           subtitle='Universos genicos distintos; comparacao de listas observadas, nao acuracia',
           x='Tipo celular',y='Jaccard')+theme_minimal(base_size=10)+
      theme(legend.position='bottom')
    ggplot2::ggsave(file.path(OUT,'08h_jaccard_dtangle.png'),g,width=13,height=7,dpi=160,bg='white')
  }
}

cat('\n08h FINALIZADO. Resultados em: ',OUT,'\n',sep='')
cat('Compare 08h_status.tsv, 08h_dtangle_top100_vs_08e_08g.tsv e 08h_top300_vs_04.tsv.\n')
cat('IMPORTANTE: nao ha informacao suficiente para igualar o universo de ~9,3 mil genes do orientador.\n')
