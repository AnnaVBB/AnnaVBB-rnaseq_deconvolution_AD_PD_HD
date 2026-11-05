#!/usr/bin/env Rscript
# =============================================================================
# 08e | REPRODUCAO DO DTANGLE NA CONFIGURACAO INFORMADA PELO ORIENTADOR
# Projeto: deconvolucao de bulk RNA-seq AD / HD / PD com Mathys 2019
# =============================================================================
# OBJETIVO:
#   Comparar com o resultado do orientador, isolando a diferenca de configuracao
#   do dtangle. Nao substitui as etapas 08, 08c ou 08d.
#
# PARAMETROS INFORMADOS PELO ORIENTADOR:
#   - Bulk e referencia em log2(CPM + 1)
#   - TODOS os genes GeneID compartilhados (nao o arquivo de top300)
#   - Um perfil medio por tipo celular (oito tipos)
#   - dtangle seleciona 100 marcadores por tipo, automaticamente
#
# PARAMETROS EXPLICITADOS NESTA REPRODUCAO:
#   - data_type = "rna-seq" (gamma calculado pelo pacote)
#   - marker_method = "ratio" (padrao do dtangle)
#   - markers = NULL; n_markers = 100L
#   - sem normalizacao quantilica (o orientador nao informou quantilica)
#   - tamanho da biblioteca bulk = library_size_full_input do passo 05
#
# AMBIGUIDADE QUE RESTA CONFIRMAR COM O ORIENTADOR:
#   "Media por tipo" pode significar:
#   (1) pooled_counts: somar contagens de todos os nucleos do tipo; CPM pela
#       biblioteca SOMADA do tipo; depois log2(CPM + 1).
#       Replica a representacao 'pooled' das etapas 08c/08d.
#   (2) mean_nucleus_cpm: normalizar cada nucleo em CPM usando sua biblioteca
#       completa, tirar a media por tipo e depois log2(CPM_medio + 1).
#   O script executa AMBAS para tornar a ambiguidade auditavel. A segunda e
#   uma sensibilidade, nao necessariamente igual a implementacao do orientador.
#
# SAIDAS: results/deconvolution/dtangle/08e_replicacao_orientador/
#   08e_run_summary.tsv, 08e_celltype_summary.tsv,
#   08e_proportions_long_all.tsv, 08e_marker_assignments.tsv,
#   08e_comparison_with_08c.tsv, 08e_status.tsv e plots/.
#
# USO (na raiz do repositorio, com renv ativado):
#   Rscript scripts/phase02_deconvolution/08e_reproduzir_dtangle_orientador.R
#   DTANGLE08E_DATASET=GSE64810 Rscript scripts/phase02_deconvolution/08e_reproduzir_dtangle_orientador.R
#   DTANGLE08E_MODE=pooled_counts Rscript scripts/phase02_deconvolution/08e_reproduzir_dtangle_orientador.R
#
# NAO atualiza renv.lock e NAO altera outputs 06/07/08/08c/08d.
# =============================================================================
options(stringsAsFactors = FALSE)

PKGS <- c("data.table", "Matrix", "dtangle")
missing_pkgs <- PKGS[!vapply(PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Pacotes ausentes: ", paste(missing_pkgs, collapse = ", "))
suppressPackageStartupMessages({library(data.table); library(Matrix)})

CELL_TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
ALL_DATASETS <- c("ROSMAP", "GSE125583", "GSE64810", "GSE68719")
EXPECTED_N <- c(ROSMAP=144L, GSE125583=289L, GSE64810=69L, GSE68719=73L)
REF_MODES <- c("pooled_counts", "mean_nucleus_cpm")
N_MARKERS <- 100L
MARKER_METHOD <- "ratio"
DATA_TYPE <- "rna-seq"
SEED <- 20261022L
SUM_TOL <- 1e-6
OUTPUT_ROOT <- "results/deconvolution/dtangle/08e_replicacao_orientador"
REF_ROOT <- "data/reference/mathys2019/phase02_prepared/reference"
BULK_ROOT <- "data/processed/deconvolution"
REF_MATRIX <- file.path(REF_ROOT, "mathys2019_pathology_free_geneid_counts_sparse.rds")
REF_META <- file.path(REF_ROOT, "mathys2019_pathology_free_cell_metadata_aligned.tsv")
OLD_08C <- "results/deconvolution/dtangle/08c_sensitivity/08c_proportions_long_all.tsv"

selected <- trimws(Sys.getenv("DTANGLE08E_DATASET", ""))
DATASETS <- if (nzchar(selected)) trimws(strsplit(selected, ",", fixed=TRUE)[[1]]) else ALL_DATASETS
if (anyNA(DATASETS) || !length(DATASETS) || any(!DATASETS %in% ALL_DATASETS)) {
  stop("DTANGLE08E_DATASET: selecione ROSMAP,GSE125583,GSE64810,GSE68719")
}
DATASETS <- unique(DATASETS)
selected_modes <- trimws(Sys.getenv("DTANGLE08E_MODE", ""))
MODES <- if (nzchar(selected_modes)) trimws(strsplit(selected_modes, ",", fixed=TRUE)[[1]]) else REF_MODES
if (anyNA(MODES) || !length(MODES) || any(!MODES %in% REF_MODES)) {
  stop("DTANGLE08E_MODE deve ser pooled_counts ou mean_nucleus_cpm")
}
MODES <- unique(MODES)
dir.create(OUTPUT_ROOT, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUTPUT_ROOT, "plots"), recursive=TRUE, showWarnings=FALSE)

write_tsv <- function(x, basename) {
  data.table::fwrite(x, file.path(OUTPUT_ROOT, basename), sep="\t", na="NA", quote=FALSE)
}
check_ids <- function(x, tag) {
  if (is.null(x) || !length(x) || anyNA(x) || any(!nzchar(as.character(x))) || anyDuplicated(x)) {
    stop(tag, ": IDs ausentes/vazios/duplicados")
  }
}
check_counts <- function(x, tag) {
  if (is.null(dim(x)) || is.null(rownames(x)) || is.null(colnames(x))) stop(tag, ": matriz sem dimnames")
  check_ids(rownames(x), paste(tag, "GeneIDs")); check_ids(colnames(x), paste(tag, "amostras"))
  if (inherits(x,"sparseMatrix")) vals <- x@x else vals <- as.numeric(x)
  if (any(!is.finite(vals)) || any(vals < 0) || any(abs(vals - round(vals)) > 1e-7)) {
    stop(tag, ": counts nao inteiros, negativos ou nao finitos")
  }
}
for (p in c(REF_MATRIX, REF_META)) if (!file.exists(p)) stop("Falta arquivo: ", p)

cat("08e | Carregando Mathys pathology-free main24...\n")
sc <- readRDS(REF_MATRIX)
if (!inherits(sc,"dgCMatrix")) sc <- methods::as(sc,"dgCMatrix")
check_counts(sc, "Mathys")
if (nrow(sc) != 17394L || ncol(sc) != 35140L) {
  stop("Referencia Mathys difere de 17.394 genes x 35.140 nucleos")
}
md <- data.table::fread(REF_META)
if (!all(c("TAG","donor_id","cell_type") %in% names(md))) stop("Metadata Mathys incompleto")
md[, `:=`(TAG=as.character(TAG), donor_id=as.character(donor_id), cell_type=as.character(cell_type))]
check_ids(md$TAG, "TAG metadata")
ii <- match(colnames(sc), md$TAG)
if (anyNA(ii)) stop("TAG da referencia nao localizado no metadata")
md <- md[ii]
if (!identical(md$TAG, colnames(sc)) || data.table::uniqueN(md$donor_id) != 24L ||
    !setequal(unique(md$cell_type), CELL_TYPES)) stop("Integridade da referencia Mathys main24 falhou")
idx <- match(md$cell_type, CELL_TYPES)
if (anyNA(idx)) stop("Tipo celular desconhecido na referencia")
type_sizes <- tabulate(idx, nbins=length(CELL_TYPES))
cat("Nucleos por tipo: ", paste(paste(CELL_TYPES,type_sizes,sep="="),collapse=", "), "\n", sep="")

# Calcular DOIS tipos diferentes de perfil medio, ainda na escala CPM linear.
# Computar a multiplicacao matricial com matrizes esparsas evita criar
# uma matriz densa genes x 35140 nucleos.
nuc_lib <- as.numeric(Matrix::colSums(sc))
if (any(!is.finite(nuc_lib)) || any(nuc_lib <= 0)) stop("Nucleo com library size <= 0")

# Variante 1 (08c pooled): contagens somadas por tipo, CPM pelo total do tipo.
A <- Matrix::sparseMatrix(i=seq_len(ncol(sc)), j=idx,
                          x=rep(1, ncol(sc)),
                          dims=c(ncol(sc),length(CELL_TYPES)))
pooled_counts <- sc %*% A
colnames(pooled_counts) <- CELL_TYPES
rownames(pooled_counts) <- rownames(sc)
pooled_lib <- as.numeric(Matrix::colSums(pooled_counts))
pooled_cpm <- sweep(as.matrix(pooled_counts), 2L, pooled_lib, "/") * 1e6

# Variante 2: CPM de CADA nucleo, media desses CPMs dentro do tipo.
# Essa media considera o mesmo peso para cada nucleo (nao para cada doador).
weights <- 1 / (nuc_lib * type_sizes[idx])
B <- Matrix::sparseMatrix(i=seq_len(ncol(sc)), j=idx, x=weights,
                          dims=c(ncol(sc),length(CELL_TYPES)))
mean_nucleus_cpm <- as.matrix(sc %*% B) * 1e6
colnames(mean_nucleus_cpm) <- CELL_TYPES
rownames(mean_nucleus_cpm) <- rownames(sc)
stopifnot(identical(rownames(pooled_cpm),rownames(mean_nucleus_cpm)))
if (any(!is.finite(pooled_cpm)) || any(!is.finite(mean_nucleus_cpm))) stop("Referencia CPM invalida")
ref_variants <- list(pooled_counts=pooled_cpm, mean_nucleus_cpm=mean_nucleus_cpm)
write_tsv(data.table(cell_type=CELL_TYPES, n_nuclei=as.integer(type_sizes),
                     pooled_library_counts=pooled_lib,
                     median_library_counts_per_nucleus=as.numeric(tapply(nuc_lib,idx,median))),
          "08e_reference_support.tsv")
rm(A,B,pooled_counts,pooled_cpm,mean_nucleus_cpm,sc)
invisible(gc(FALSE))

# Importante: queremos TODOS os genes compartilhados, nao apenas o top300.
# Os genes que sao zero no bulk inteiro ou na referencia inteira nao sao
# informativos para selecionar marcadores; remover de modo uniforme por coorte.
summary_rows <- list(); celltype_rows <- list(); pred_rows <- list()
marker_rows <- list(); compare_rows <- list(); status_rows <- list()

read_old <- if (file.exists(OLD_08C)) data.table::fread(OLD_08C) else NULL

for (ds in DATASETS) {
  cat("\n------------------------------\n08e | ", ds, "\n", sep="")
  path_counts <- file.path(BULK_ROOT, ds, "bulk_counts_reference_intersection.rds")
  path_meta <- file.path(BULK_ROOT, ds, "metadata_deconvolution.tsv")
  path_lib <- file.path(BULK_ROOT, ds, "library_sizes.tsv")
  tryCatch({
    missing <- c(path_counts,path_meta,path_lib)[!file.exists(c(path_counts,path_meta,path_lib))]
    if (length(missing)) stop("Arquivos da etapa 05 ausentes: ",paste(missing,collapse="; "))
    bulk <- readRDS(path_counts)
    check_counts(bulk, paste(ds,"bulk"))
    if (ncol(bulk) != EXPECTED_N[[ds]]) stop(ds, ": numero inesperado de amostras")
    bm <- data.table::fread(path_meta)
    ls <- data.table::fread(path_lib)
    if (!all(c("sample_id","group") %in% names(bm)) ||
        !all(c("sample_id","library_size_full_input") %in% names(ls))) {
      stop("metadata ou library_sizes.tsv invalido")
    }
    bm[, sample_id:=as.character(sample_id)]
    ls[, sample_id:=as.character(sample_id)]
    check_ids(bm$sample_id, "sample IDs metadata")
    check_ids(ls$sample_id, "sample IDs library")
    if (!setequal(bm$sample_id,colnames(bulk))) stop("metadata e counts divergentes")
    bm <- bm[match(colnames(bulk),bm$sample_id)]
    libs <- as.numeric(ls$library_size_full_input[match(colnames(bulk),ls$sample_id)])
    if (anyNA(libs) || any(!is.finite(libs)) || any(libs<=0)) stop("Library sizes invalidas")
    if (any(libs + 1e-6 < as.numeric(Matrix::colSums(bulk)))) {
      stop("Biblioteca completa inferior aos counts presentes")
    }
    genes <- rownames(bulk)[rownames(bulk) %in% rownames(ref_variants[[1]])]
    if (length(genes) < 1000L) stop("Menos de 1000 genes compartilhados")
    b <- bulk[genes, , drop=FALSE]
    # Remover genes sem contagens ao longo do bulk todo e referencia inteira.
    keep <- as.numeric(Matrix::rowSums(b)) > 0 &
            rowSums(ref_variants[[1]][genes, , drop=FALSE]) > 0
    genes <- genes[keep]
    if (length(genes) < 1000L) stop("Menos de 1000 genes apos filtro")
    b <- as.matrix(b[keep,,drop=FALSE])
    if (!identical(rownames(b),genes)) stop("Genes bulk desalinhados")
    # Calcular CPM com denominador do arquivo completo da etapa 05.
    # Se esse arquivo for o full-reference-intersection e nao o genoma
    # completo, library_sizes.tsv e a fonte auditada apropriada.
    bulk_log <- t(log2(sweep(b, 2L, libs, "/") * 1e6 + 1))
    if (!identical(rownames(bulk_log),bm$sample_id) ||
        !identical(colnames(bulk_log),genes) || any(!is.finite(bulk_log))) {
      stop("Matriz bulk log2CPM com ordem ou valores invalidos")
    }
    rm(b,bulk)
    pure <- setNames(lapply(seq_along(CELL_TYPES),function(k) as.integer(k)), CELL_TYPES)

    for (mode in MODES) {
      cat("  ", mode, ": ",length(genes)," genes; n_markers=",N_MARKERS,"\n", sep="")
      ref_log <- t(log2(ref_variants[[mode]][genes,,drop=FALSE] + 1))
      if (!identical(colnames(bulk_log),colnames(ref_log))) stop("Erro de GeneID/ordem")
      if (!identical(rownames(ref_log),CELL_TYPES)) stop("Erro de tipo celular/ordem")
      t0 <- proc.time()[["elapsed"]]
      set.seed(SEED)
      fit <- dtangle::dtangle(Y=bulk_log, references=ref_log,
                             pure_samples=pure,
                             markers=NULL, n_markers=N_MARKERS,
                             marker_method=MARKER_METHOD, data_type=DATA_TYPE)
      elapsed <- as.numeric(proc.time()[["elapsed"]]-t0)
      p <- as.matrix(fit$estimates)
      if (nrow(p) != length(bm$sample_id) || ncol(p) != length(CELL_TYPES)) {
        stop("Dimensoes inesperadas nas estimativas")
      }
      if (!is.null(rownames(p)) && !identical(rownames(p),bm$sample_id)) {
        if (!setequal(rownames(p),bm$sample_id)) stop("IDs de amostras divergentes")
        p <- p[bm$sample_id,,drop=FALSE]
      }
      if (!is.null(colnames(p)) && !identical(colnames(p),CELL_TYPES)) {
        if (!setequal(colnames(p),CELL_TYPES)) stop("Nomes dos tipos divergentes")
        p <- p[,CELL_TYPES,drop=FALSE]
      }
      rownames(p) <- bm$sample_id; colnames(p) <- CELL_TYPES
      if (any(!is.finite(p)) || any(p < -1e-8) ||
          max(abs(rowSums(p)-1)) > SUM_TOL) stop("Estimativas nao validas; sem correcao silenciosa")
      if (any(p < 0)) stop("Existem proporcoes negativas")
      marker_indices <- fit$markers
      if (length(marker_indices) != length(CELL_TYPES)) stop("Numero de listas de marcadores inesperado")
      for (i in seq_along(CELL_TYPES)) {
        mi <- as.integer(marker_indices[[i]])
        if (length(mi) != N_MARKERS || anyNA(mi) || any(mi<1L | mi>length(genes))) {
          stop("Marcadores automaticos nao atingiram 100 para ",CELL_TYPES[[i]])
        }
        marker_rows[[length(marker_rows)+1L]] <- data.table(
          dataset=ds, reference_mode=mode, cell_type=CELL_TYPES[[i]],
          rank=seq_along(mi), gene_id=genes[mi], gene_column_index=mi)
      }
      est <- data.table::as.data.table(p, keep.rownames="sample_id")
      est_long <- data.table::melt(est, id.vars="sample_id", variable.name="cell_type", value.name="proportion")
      est_long[, `:=`(dataset=ds, reference_mode=mode)]
      est_long[, group:=bm$group[match(sample_id,bm$sample_id)]]
      pred_rows[[length(pred_rows)+1L]] <- est_long[,.(dataset,reference_mode,sample_id,group,cell_type,proportion)]
      cell_sum <- est_long[,.(n_samples=.N, mean_proportion=mean(proportion),
                                median_proportion=median(proportion),
                                sd_pp=100*stats::sd(proportion),
                                fraction_exact_zero=mean(proportion==0),
                                range_pp=100*(max(proportion)-min(proportion))),
                           by=.(dataset,reference_mode,cell_type)]
      celltype_rows[[length(celltype_rows)+1L]] <- cell_sum
      summary_rows[[length(summary_rows)+1L]] <- data.table(
        dataset=ds, reference_mode=mode, n_samples=nrow(p), genes_used=length(genes),
        n_reference_profiles=8L, n_markers_per_type=N_MARKERS,
        data_type=DATA_TYPE, marker_method=MARKER_METHOD,
        gamma=as.numeric(fit$gamma),
        mean_In_percent=100*mean(p[,"In"]),
        mean_Opc_percent=100*mean(p[,"Opc"]),
        mean_vascular_percent=100*mean(p[,"End"]+p[,"Per"]),
        max_sum_error=max(abs(rowSums(p)-1)),
        elapsed_seconds=elapsed, dtangle_version=as.character(utils::packageVersion("dtangle")))
      cat(sprintf("      End+Per=%.2f%%; In=%.2f%%; Opc=%.2f%%\n",
                  100*mean(p[,"End"]+p[,"Per"]),100*mean(p[,"In"]),100*mean(p[,"Opc"])))
      # Comparar amostra A amostra com o 08c pooled_raw anterior, se presente.
      if (!is.null(read_old) && mode=="pooled_counts" &&
          all(c("dataset","condition","sample_id","cell_type","proportion") %in% names(read_old))) {
        prior <- read_old[read_old$dataset == ds & read_old$condition == "pooled_raw"]
        if (nrow(prior)) {
          paired <- merge(est_long[,.(sample_id,cell_type,new_estimate=proportion)],
                          prior[,.(sample_id,cell_type,old_08c_estimate=proportion)],
                          by=c("sample_id","cell_type"),all=TRUE)
          if (nrow(paired) != nrow(p)*length(CELL_TYPES) ||
              anyNA(paired$new_estimate) || anyNA(paired$old_08c_estimate)) {
            stop("Comparacao 08e/08c nao possui pares completos")
          }
          cmp <- paired[,.(mae_pp=100*mean(abs(new_estimate-old_08c_estimate)),
                           mean_new=mean(new_estimate),mean_old=mean(old_08c_estimate)),
                        by=cell_type]
          cmp[,dataset:=ds]
          compare_rows[[length(compare_rows)+1L]] <- cmp
        }
      }
      status_rows[[length(status_rows)+1L]] <- data.table(
        dataset=ds,reference_mode=mode,check="fit_integrity",status="PASS",
        info=paste(nrow(p),"amostras;",length(genes),"genes;",N_MARKERS,"marcadores/tipo"))
    }
    rm(bm,libs,bulk_log)
    invisible(gc(FALSE))
  }, error=function(e) {
    message("FALHA em ", ds, ": ",conditionMessage(e))
    status_rows[[length(status_rows)+1L]] <<- data.table(
      dataset=ds, reference_mode="ALL",check="execution",status="FAIL",info=conditionMessage(e))
  })
  # Checkpoint incremental apos cada coorte: preserva as saidas anteriores.
  if (length(pred_rows)) write_tsv(rbindlist(pred_rows,fill=TRUE),"08e_proportions_long_all.tsv")
  if (length(summary_rows)) write_tsv(rbindlist(summary_rows,fill=TRUE),"08e_run_summary.tsv")
  if (length(celltype_rows)) write_tsv(rbindlist(celltype_rows,fill=TRUE),"08e_celltype_summary.tsv")
  if (length(marker_rows)) write_tsv(rbindlist(marker_rows,fill=TRUE),"08e_marker_assignments.tsv")
  if (length(compare_rows)) write_tsv(rbindlist(compare_rows,fill=TRUE),"08e_comparison_with_08c.tsv")
  if (length(status_rows)) write_tsv(rbindlist(status_rows,fill=TRUE),"08e_status.tsv")
}

# Gráficos descritivos. A falha em PNG NAO invalida as estimativas.
if (length(summary_rows) && requireNamespace("ggplot2",quietly=TRUE)) {
  tryCatch({
    d <- rbindlist(summary_rows)
    d_long <- data.table::melt(d,
      id.vars=c("dataset","reference_mode"),
      measure.vars=c("mean_In_percent","mean_Opc_percent","mean_vascular_percent"),
      variable.name="population",value.name="mean_percent")
    d_long[,population:=factor(population,
      levels=c("mean_In_percent","mean_Opc_percent","mean_vascular_percent"),
      labels=c("In","Opc","End + Per"))]
    g <- ggplot2::ggplot(d_long,ggplot2::aes(x=dataset,y=mean_percent,fill=reference_mode))+
      ggplot2::geom_col(position=ggplot2::position_dodge(width=.75),width=.68)+
      ggplot2::facet_wrap(~population,scales="free_y")+
      ggplot2::theme_minimal(base_size=12)+
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,hjust=1),
                     plot.margin=ggplot2::margin(15,18,18,15))+
      ggplot2::labs(title="08e: dtangle com todos os genes e 100 marcadores/tipo",
                   subtitle="Duas definicoes auditaveis de perfil medio por tipo",
                   x="Coorte",y="Proporcao media estimada (%)",fill="Referencia")
    ggplot2::ggsave(file.path(OUTPUT_ROOT,"plots","08e_means_by_reference_mode.png"),
                    g,width=13,height=6.5,dpi=170)
  },error=function(e) warning("Falha apenas ao gerar figura: ",conditionMessage(e)))
}

capture.output(sessionInfo(),file=file.path(OUTPUT_ROOT,"08e_sessionInfo.txt"))
if (length(status_rows) && any(rbindlist(status_rows,fill=TRUE)$status=="FAIL")) {
  stop("08e terminou com falhas. Verifique 08e_status.tsv; resultados anteriores preservados.")
}
cat("\n=== 08e CONCLUIDO ===\nResultados: ", OUTPUT_ROOT,"\n",sep="")
if (length(summary_rows)) print(rbindlist(summary_rows)[,.(dataset,reference_mode,genes_used,mean_In_percent,mean_Opc_percent,mean_vascular_percent)])
