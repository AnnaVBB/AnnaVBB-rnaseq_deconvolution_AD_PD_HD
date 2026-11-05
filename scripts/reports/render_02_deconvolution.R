#!/usr/bin/env Rscript
# =============================================================================
# Renderizacao do relatorio integral da Fase 02: 02_deconvolution.Rmd
# Uso (na raiz do repositorio):
#   Rscript scripts/reports/render_02_deconvolution.R
#   Rscript scripts/reports/render_02_deconvolution.R /caminho/do/repositorio
# Variaveis opcionais:
#   DECONV_ROOT=/caminho/repo
#   DECONV_FULL_TABLES=1  (inclui todas as tabelas adicionais - HTML grande)
#   DECONV_MAX_ROWS=60    (linhas por tabela; so preview, TSV integral preservado)
#   DECONV_MAX_FIGURES=80 (numero maximo de figuras por secao)
# Nao executa novamente a deconvolucao nem sobrescreve as saidas dos scripts.
# =============================================================================
args <- commandArgs(trailingOnly=TRUE)
root_arg <- if(length(args)) args[[1L]] else Sys.getenv('DECONV_ROOT', unset='.')
root <- normalizePath(root_arg, winslash='/', mustWork=TRUE)
if(!dir.exists(file.path(root,'scripts','phase02_deconvolution')))
  stop('A raiz deve conter scripts/phase02_deconvolution: ', root)
input <- file.path(root,'scripts','reports','02_deconvolution.Rmd')
if(!file.exists(input)) stop('Coloque o Rmd em scripts/reports/02_deconvolution.Rmd: ',input)
for(pkg in c('rmarkdown','knitr','data.table')) {
  if(!requireNamespace(pkg,quietly=TRUE)) stop('Instale o pacote R ausente: ',pkg)
}
rows <- suppressWarnings(as.integer(Sys.getenv('DECONV_MAX_ROWS','60')))
if(is.na(rows) || rows<1L) stop('DECONV_MAX_ROWS deve ser inteiro >= 1')
figs <- suppressWarnings(as.integer(Sys.getenv('DECONV_MAX_FIGURES','80')))
if(is.na(figs) || figs<1L) stop('DECONV_MAX_FIGURES deve ser inteiro >= 1')
all_tables <- tolower(Sys.getenv('DECONV_FULL_TABLES','0')) %in% c('1','yes','true','sim')
out <- file.path(root,'results','reports','deconvolution')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cat('Raiz: ',root,'\nEntrada: ',input,'\nDestino: ',out,'\n',sep='')
ans <- rmarkdown::render(
  input=input,
  output_format='html_document',
  output_file='02_deconvolution.html',
  output_dir=out,
  params=list(project_root=root,max_rows_table=rows,
              include_all_tables=all_tables,max_figures_per_stage=figs),
  knit_root_dir=root,
  envir=new.env(parent=globalenv()),
  quiet=FALSE
)
if(!file.exists(ans)) stop('Renderizacao nao gerou o HTML: ',ans)
cat('\nRELATORIO CONCLUIDO: ',normalizePath(ans,winslash='/'),
    '\nManifesto: ',file.path(out,'02_deconvolution_file_manifest.tsv'),'\n',sep='')
