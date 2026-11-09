#!/usr/bin/env Rscript
# Inventario exploratorio dos arquivos de figuras e locais de exportacao.
# Nao certifica eixos, cores, histogramas ou DPI sem inspecao grafica.
options(stringsAsFactors=FALSE)
if (!requireNamespace('data.table',quietly=TRUE)) stop('Necessario data.table')
library(data.table)
ROOT <- normalizePath(Sys.getenv('DECONV_ROOT','.'),mustWork=TRUE)
P <- function(...) file.path(ROOT,...)
OUT <- P('results','reports','plot_style_audit')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(x,f) fwrite(x,file.path(OUT,f),sep='\t',quote=FALSE,na='NA')
files <- list.files(P('scripts'),recursive=TRUE,full.names=TRUE,pattern='\\.(R|Rmd)$',ignore.case=TRUE)
patterns <- list(ggsave='ggsave\\s*\\(', png_device='\\bpng\\s*\\(', pdf_device='\\bpdf\\s*\\(',
 hist='geom_histogram\\s*\\(|\\bhist\\s*\\(',
 explicit_labels='labs\\s*\\(|xlab\\s*\\(|ylab\\s*\\(',
 theme='theme_.*\\s*\\(', source_shared_style='00_plot_style\\.R',
 cell_palette='scale_(color|fill)_manual\\s*\\(')
items <- lapply(files,function(f){
 lines <- tryCatch(readLines(f,warn=FALSE,encoding='UTF-8'),error=function(e)character())
 hits <- lapply(patterns,function(p) which(grepl(p,lines,perl=TRUE)))
 data.table(script=substring(f,nchar(ROOT)+2L),
            total_lines=length(lines),
            n_ggsave=length(hits$ggsave),n_png=length(hits$png_device),
            n_pdf=length(hits$pdf_device),n_histogram=length(hits$hist),
            n_label_calls=length(hits$explicit_labels),
            uses_shared_style=length(hits$source_shared_style)>0,
            review_required=length(hits$ggsave)+length(hits$png_device)+length(hits$pdf_device)+length(hits$hist)>0)
})
if(length(items)) write_tsv(rbindlist(items),'99_plot_script_inventory.tsv')
# Ler largura/altura do PNG usando IHDR (nao implica resolucao de impressao).
pngs <- if(dir.exists(P('results'))) list.files(P('results'),recursive=TRUE,full.names=TRUE,pattern='\\.png$',ignore.case=TRUE) else character()
png_info <- function(f){
 con <- file(f,'rb');on.exit(close(con),add=TRUE)
 raw <- readBin(con,what='raw',n=24L)
 if(length(raw)<24L)return(c(NA_real_,NA_real_))
 be <- function(v)sum(as.integer(v)*256^(3:0))
 c(be(raw[17:20]),be(raw[21:24]))
}
if(length(pngs)) {
 out <- rbindlist(lapply(pngs,function(f){d<-png_info(f)
   data.table(file=f,width_px=d[1],height_px=d[2],
       manual_review='Conferir eixo X/Y, legenda, paleta, contorno histograma e exportacao a 300 dpi')
 }))
 write_tsv(out,'99_png_dimensions.tsv')
}
writeLines(c('Este inventario e exploratorio.','A existencia de labels em codigo nao garante X e Y em todos os plots.',
 'Largura em pixels nao certifica DPI de publicacao.',
 'A integracao do 00_plot_style.R e o reprocessamento de cada script ainda sao necessarios.'),
 file.path(OUT,'99_LIMITATIONS.txt'))
cat('Inventario gerado: ',OUT,'\n',sep='')
