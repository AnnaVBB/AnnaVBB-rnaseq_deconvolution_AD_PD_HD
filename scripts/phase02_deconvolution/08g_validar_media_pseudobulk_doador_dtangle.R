#!/usr/bin/env Rscript
# =============================================================================
# 08g | AUDITORIA DA MEDIA DE CPM POR DOADOR PARA DTANGLE
# Projeto: AD / PD / HD, referencia Mathys 2019
# =============================================================================
# A configuracao 08e pooled_counts agrega todos os nucleos por tipo (sem
# equalizar o peso dos doadores). A 08e mean_nucleus_cpm tira a media
# de CPM por nucleo (tambem sem equalizar doadores).
#
# Nesta auditoria implementamos a definicao confirmada pelo orientador:
#   1) somar counts por doador x tipo => pseudobulk
#   2) dividir pela biblioteca de cada pseudobulk e multiplicar por 1e6
#   3) media aritmetica dos CPMs dos doadores OBSERVADOS por tipo
#   4) log2(CPM medio + 1) antes do dtangle
# Os tipos End/Per NAO sao imputados como zero em doadores sem nucleos.
#
# Mantem exatamente as demais configuracoes do 08e:
#   - todos os GeneIDs compartilhados entre o bulk e a referencia;
#   - CPM bulk com biblioteca completa, registrada na etapa 05;
#   - marker_method='ratio', n_markers=100, data_type='rna-seq';
#   - gamma estimado pelo dtangle; sem normalizacao quantilica.
#
# Execucao (a partir da raiz do projeto):
#   Rscript scripts/phase02_deconvolution/08g_validar_media_pseudobulk_doador_dtangle.R
# Uma coorte (opcional):
#   DTANGLE08G_DATASET=GSE64810 Rscript scripts/phase02_deconvolution/08g_validar_media_pseudobulk_doador_dtangle.R
#
# Saida isolada: results/deconvolution/dtangle/08g_media_por_doador/
# Nao substitui ou sobrescreve os resultados 08e/08f/09/12.
# =============================================================================
options(stringsAsFactors = FALSE)

PKGS <- c("data.table", "Matrix", "dtangle")
missing_pkgs <- PKGS[!vapply(PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Pacotes ausentes: ", paste(missing_pkgs, collapse = ", "))
suppressPackageStartupMessages({library(data.table); library(Matrix)})

CELL_TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
ALL_DATASETS <- c("ROSMAP", "GSE125583", "GSE64810", "GSE68719")
EXPECTED_N <- c(ROSMAP=144L, GSE125583=289L, GSE64810=69L, GSE68719=73L)
REF_MODES <- c("donor_pseudobulk_mean_cpm")
N_MARKERS <- 100L
MARKER_METHOD <- "ratio"
DATA_TYPE <- "rna-seq"
SEED <- 20261022L
SUM_TOL <- 1e-6
OUTPUT_ROOT <- "results/deconvolution/dtangle/08g_media_por_doador"
REF_ROOT <- "data/reference/mathys2019/phase02_prepared/reference"
BULK_ROOT <- "data/processed/deconvolution"
REF_MATRIX <- file.path(REF_ROOT, "mathys2019_pathology_free_geneid_counts_sparse.rds")
REF_META <- file.path(REF_ROOT, "mathys2019_pathology_free_cell_metadata_aligned.tsv")
OLD_08C <- "results/deconvolution/dtangle/08c_sensitivity/08c_proportions_long_all.tsv"
OLD_08E_P <- "results/deconvolution/dtangle/08e_replicacao_orientador/08e_proportions_long_all.tsv"
OLD_08E_M <- "results/deconvolution/dtangle/08e_replicacao_orientador/08e_marker_assignments.tsv"

selected <- trimws(Sys.getenv("DTANGLE08G_DATASET", ""))
DATASETS <- if (nzchar(selected)) trimws(strsplit(selected, ",", fixed=TRUE)[[1]]) else ALL_DATASETS
if (anyNA(DATASETS) || !length(DATASETS) || any(!DATASETS %in% ALL_DATASETS)) {
  stop("DTANGLE08G_DATASET: selecione ROSMAP,GSE125583,GSE64810,GSE68719")
}
DATASETS <- unique(DATASETS)
selected_modes <- trimws(Sys.getenv("DTANGLE08G_MODE", ""))
MODES <- if (nzchar(selected_modes)) trimws(strsplit(selected_modes, ",", fixed=TRUE)[[1]]) else REF_MODES
if (anyNA(MODES) || !length(MODES) || any(!MODES %in% REF_MODES)) {
  stop("DTANGLE08G_MODE deve ser donor_pseudobulk_mean_cpm")
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

cat("08g | Carregando Mathys pathology-free main24...\n")
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

# Construir pseudobulks doador x tipo sem materializar 35.140 colunas densas.
# Para End/Per, doadores sem aquele tipo sao ausentes, NAO perfis zero.
# Uma coluna da matriz P representa uma combinacao observada doador x tipo.
md[, pair_id := paste(donor_id, cell_type, sep="__08g__")]
pair_levels <- unique(md$pair_id)
pair_index <- match(md$pair_id, pair_levels)
pair_types <- md$cell_type[match(pair_levels, md$pair_id)]
pair_donors <- md$donor_id[match(pair_levels, md$pair_id)]
if (length(pair_levels) != 181L) {
  warning("Numero de pseudobulks doador x tipo difere dos 181 esperados: ", length(pair_levels))
}

P <- Matrix::sparseMatrix(
  i=seq_len(ncol(sc)), j=pair_index, x=rep(1, ncol(sc)),
  dims=c(ncol(sc),length(pair_levels))
)
pb_counts <- sc %*% P
colnames(pb_counts) <- pair_levels
rownames(pb_counts) <- rownames(sc)
pb_lib <- as.numeric(Matrix::colSums(pb_counts))
if (any(!is.finite(pb_lib)) || any(pb_lib <= 0)) {
  stop("Pseudobulk com biblioteca nula ou invalida")
}
pb_cpm <- sweep(as.matrix(pb_counts), 2L, pb_lib, "/") * 1e6
if (max(abs(colSums(pb_cpm) - 1e6)) > 1e-3) {
  stop("CPM dos pseudobulks nao soma 1e6: denominador incorreto")
}

# Cada pseudobulk de um tipo recebe peso 1/N_de_doadores_observados.
pair_type_index <- match(pair_types, CELL_TYPES)
if (anyNA(pair_type_index)) stop("Tipo inesperado nos pseudobulks")
n_donors_observed <- tabulate(pair_type_index, nbins=length(CELL_TYPES))
if (any(n_donors_observed < 2L)) stop("Tipo celular representado por menos de 2 doadores")
W <- Matrix::sparseMatrix(
  i=seq_along(pair_levels), j=pair_type_index,
  x=1/n_donors_observed[pair_type_index],
  dims=c(length(pair_levels), length(CELL_TYPES))
)
donor_mean_cpm <- as.matrix(pb_cpm %*% W)
colnames(donor_mean_cpm) <- CELL_TYPES
rownames(donor_mean_cpm) <- rownames(sc)
if (any(!is.finite(donor_mean_cpm)) || any(donor_mean_cpm < 0)) {
  stop("CPM medio por doador invalido")
}
if (max(abs(colSums(donor_mean_cpm) - 1e6)) > 1e-3) {
  stop("Perfil medio nao soma 1e6: pesos dos doadores incorretos")
}
stopifnot(length(unique(pair_donors))==24L)
ref_variants <- list(donor_pseudobulk_mean_cpm=donor_mean_cpm)

# Auditoria de pesos: 1/N de doadores por tipo, sem ponderacao por numero de nucleos.
write_tsv(data.table(
  cell_type=CELL_TYPES,
  n_observed_donors=as.integer(n_donors_observed),
  n_nuclei=as.integer(type_sizes),
  n_pseudobulks=as.integer(n_donors_observed),
  effective_weight_per_donor=1/n_donors_observed,
  comment="Media CPM dos pseudobulks de doadores observados; ausentes nao imputados"
), "08g_reference_support.tsv")
write_tsv(data.table(
  donor_id=pair_donors, cell_type=pair_types,
  library_size_pseudobulk=pb_lib,
  n_nuclei=as.integer(tabulate(pair_index,nbins=length(pair_levels)))
), "08g_pseudobulk_support.tsv")
rm(P,W,pb_counts,pb_cpm,donor_mean_cpm,sc)
invisible(gc(FALSE))

# Importante: queremos TODOS os genes compartilhados, nao apenas o top300.
# Os genes que sao zero no bulk inteiro ou na referencia inteira nao sao
# informativos para selecionar marcadores; remover de modo uniforme por coorte.
summary_rows <- list(); celltype_rows <- list(); pred_rows <- list()
marker_rows <- list(); compare_rows <- list(); status_rows <- list()

read_old <- if (file.exists(OLD_08C)) data.table::fread(OLD_08C) else NULL
read_08e_p <- if (file.exists(OLD_08E_P)) data.table::fread(OLD_08E_P) else NULL
read_08e_m <- if (file.exists(OLD_08E_M)) data.table::fread(OLD_08E_M) else NULL
vs08e_rows <- list()
marker_vs08e_rows <- list()


for (ds in DATASETS) {
  cat("\n------------------------------\n08g | ", ds, "\n", sep="")
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
      # Comparar os 100 marcadores por tipo com o 08e, se disponivel.
      # Jaccard = intersecao / uniao; Spearman no conjunto de genes comum.
      if (!is.null(read_08e_m) && all(c("dataset","reference_mode","cell_type","gene_id","rank") %in% names(read_08e_m))) {
        newly_selected <- data.table::rbindlist(marker_rows)[dataset==ds & reference_mode==mode]
        for (old_mode in c("pooled_counts","mean_nucleus_cpm")) {
          old_selected <- read_08e_m[dataset==ds & reference_mode==old_mode]
          for (ct in CELL_TYPES) {
            newct <- newly_selected[cell_type==ct]
            oldct <- old_selected[cell_type==ct]
            if (!nrow(oldct)) next
            common <- intersect(newct$gene_id, oldct$gene_id)
            j <- length(common)/length(union(newct$gene_id, oldct$gene_id))
            rnk <- if (length(common) >= 3L) suppressWarnings(stats::cor(
              newct$rank[match(common,newct$gene_id)],
              oldct$rank[match(common,oldct$gene_id)], method="spearman")) else NA_real_
            marker_vs08e_rows[[length(marker_vs08e_rows)+1L]] <- data.table(
              dataset=ds,reference_mode_new=mode,reference_mode_old=old_mode,
              cell_type=ct,n_common_markers=length(common),jaccard=j,
              spearman_rank_overlap=rnk)
          }
        }
      }
      est <- data.table::as.data.table(p, keep.rownames="sample_id")
      est_long <- data.table::melt(est, id.vars="sample_id", variable.name="cell_type", value.name="proportion")
      est_long[, `:=`(dataset=ds, reference_mode=mode)]
      est_long[, group:=bm$group[match(sample_id,bm$sample_id)]]
      pred_rows[[length(pred_rows)+1L]] <- est_long[,.(dataset,reference_mode,sample_id,group,cell_type,proportion)]
      # Comparacao amostra-a-amostra com as duas versoes anteriores de 08e.
      if (!is.null(read_08e_p) && all(c("dataset","reference_mode","sample_id","cell_type","proportion") %in% names(read_08e_p))) {
        for (old_mode in c("pooled_counts","mean_nucleus_cpm")) {
          oldpart <- read_08e_p[dataset==ds & reference_mode==old_mode,
                                .(sample_id,cell_type,old_proportion=proportion)]
          if (!nrow(oldpart)) next
          newpart <- est_long[,.(sample_id,cell_type,new_proportion=proportion)]
          matched <- merge(newpart,oldpart,by=c("sample_id","cell_type"),all=FALSE)
          if (nrow(matched)!=nrow(p)*length(CELL_TYPES)) {
            stop("Confronto com 08e incompleto: ",ds," / ",old_mode)
          }
          cmp <- matched[,.(n_samples=.N,
             mean_new_percent=100*mean(new_proportion),
             mean_old_percent=100*mean(old_proportion),
             mean_diff_pp=100*mean(new_proportion-old_proportion),
             mae_pp=100*mean(abs(new_proportion-old_proportion)),
             max_abs_diff_pp=100*max(abs(new_proportion-old_proportion)),
             spearman=suppressWarnings(stats::cor(new_proportion,old_proportion,method="spearman"))),
             by=cell_type]
          cmp[, `:=`(dataset=ds,new_mode=mode,old_mode=old_mode)]
          vs08e_rows[[length(vs08e_rows)+1L]] <- cmp
        }
      }
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
      if (!is.null(read_old) && mode=="donor_pseudobulk_mean_cpm" &&
          all(c("dataset","condition","sample_id","cell_type","proportion") %in% names(read_old))) {
        prior <- read_old[read_old$dataset == ds & read_old$condition == "pooled_raw"]
        if (nrow(prior)) {
          paired <- merge(est_long[,.(sample_id,cell_type,new_estimate=proportion)],
                          prior[,.(sample_id,cell_type,old_08c_estimate=proportion)],
                          by=c("sample_id","cell_type"),all=TRUE)
          if (nrow(paired) != nrow(p)*length(CELL_TYPES) ||
              anyNA(paired$new_estimate) || anyNA(paired$old_08c_estimate)) {
            stop("Comparacao 08g/08c nao possui pares completos")
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
  if (length(pred_rows)) write_tsv(rbindlist(pred_rows,fill=TRUE),"08g_proportions_long_all.tsv")
  if (length(summary_rows)) write_tsv(rbindlist(summary_rows,fill=TRUE),"08g_run_summary.tsv")
  if (length(celltype_rows)) write_tsv(rbindlist(celltype_rows,fill=TRUE),"08g_celltype_summary.tsv")
  if (length(marker_rows)) write_tsv(rbindlist(marker_rows,fill=TRUE),"08g_marker_assignments.tsv")
  if (length(compare_rows)) write_tsv(rbindlist(compare_rows,fill=TRUE),"08g_comparison_with_08c.tsv")
  if (length(status_rows)) write_tsv(rbindlist(status_rows,fill=TRUE),"08g_status.tsv")
  if (length(vs08e_rows)) write_tsv(rbindlist(vs08e_rows,fill=TRUE),"08g_comparison_with_08e.tsv")
  if (length(marker_vs08e_rows)) write_tsv(rbindlist(marker_vs08e_rows,fill=TRUE),"08g_marker_comparison_with_08e.tsv")
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
      ggplot2::labs(title="08g: dtangle com media de CPM entre pseudobulks de doadores",
                   subtitle="Pseudobulk por doador, CPM, media entre doadores, log2(CPM+1)",
                   x="Coorte",y="Proporcao media estimada (%)",fill="Referencia")
    ggplot2::ggsave(file.path(OUTPUT_ROOT,"plots","08g_means_by_reference_mode.png"),
                    g,width=13,height=6.5,dpi=170)
  },error=function(e) warning("Falha apenas ao gerar figura: ",conditionMessage(e)))
}

capture.output(sessionInfo(),file=file.path(OUTPUT_ROOT,"08g_sessionInfo.txt"))
if (length(status_rows) && any(rbindlist(status_rows,fill=TRUE)$status=="FAIL")) {
  stop("08g terminou com falhas. Verifique 08e_status.tsv; resultados anteriores preservados.")
}
cat("\n=== 08g CONCLUIDO ===\nResultados: ", OUTPUT_ROOT,"\n",sep="")
if (length(summary_rows)) print(rbindlist(summary_rows)[,.(dataset,reference_mode,genes_used,mean_In_percent,mean_Opc_percent,mean_vascular_percent)])
