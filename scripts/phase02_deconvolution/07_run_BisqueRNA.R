#!/usr/bin/env Rscript
# ==============================================================================
# 07_run_BisqueRNA.R — FASE 02 / TCC AD–HD–PD [v5: boxplots com margens + prior audit]
# ==============================================================================
# OBJETIVO
#   Executar a deconvolução de referência do BisqueRNA nas quatro coortes:
#     ROSMAP (AD principal), GSE125583 (AD validação),
#     GSE64810 (HD), GSE68719 (PD).
#   Referência Mathys 2019 pathology-free: 24 doadores, 35.140 núcleos,
#   oito tipos: Ex, In, Ast, Oli, Opc, Mic, End, Per.
#   Marcadores primários: união dos top300 construídos no script 04.
#
# DIFERENÇAS EM RELAÇÃO AO MuSiC2 (ETAPA 06)
#   - BisqueRNA::ReferenceBasedDecomposition(), não MarkerBasedDecomposition().
#   - Recebe Biobase::ExpressionSet com counts não transformados.
#   - A saída bulk.props contém proporções: linhas = tipos; colunas = amostras.
#   - Em ausência de pares bulk/snRNA-seq adequados à calibração supervisionada,
#     use.overlap = FALSE para TODAS as coortes.
#   - old.cpm = FALSE: o Bisque normaliza as bibliotecas ANTES de subsetting.
#
# ATENÇÃO: NORMALIZAÇÃO CORRETA E USO DE MEMÓRIA
#   O BisqueRNA, em sua implementação de referência, guarda um ExpressionSet
#   DENSO e implementa CountsToCPM() dividindo cada coluna por sua soma.
#   Alimentar somente top300 faria CPM depender do total de genes selecionados.
#
#   Para preservar as library sizes originais sem criar um objeto denso de
#   17.394 genes x 35.140 núcleos, cada matriz entregue ao Bisque tem:
#       a) counts REAIS de genes top300 disponíveis; e
#       b) UMA linha técnica '__LIBRARY_REMAINDER__', calculada como
#          (library size da matriz completa) - (soma dos counts top300).
#
#   Com old.cpm=FALSE, a API primeiro calcula CPM com o denominador correto
#   e DEPOIS remove a linha técnica, porque 'markers' contém APENAS GeneIDs.
#   Logo, essa linha NUNCA participa de sc.ref, transformação nem regressão.
#   É um recurso explícito de CONTABILIDADE, e NÃO gene biológico ou imputação.
#   Há verificações matemáticas abaixo, e sua presença nos genes finais gera
#   falha imediata. A implementação foi consultada no repositório CRAN/BisqueRNA.
#
#   Dependência importante: a ordem normalização -> subsetting é garantida na
#   API documentada quando old.cpm=FALSE. Registrar versão instalada.
#
# CONTROLE DE QUALIDADE
#   - Integridade de todos os IDs e alinhamento bulk/reference/meta.
#   - Counts inteiros, finitos, não negativos, library sizes positivas.
#   - Conservação das library sizes após incluir a linha técnica.
#   - Conservação dos valores CPM esperados (teste matemático).
#   - Oito tipos, proporções finitas, não negativas e soma ~1 por amostra.
#   - Fração de zeros / quase zeros; SD e amplitude entre amostras.
#   - Degenerescência: pouca variação + composição próxima de 12,5% cada.
#   - Degenerescência gera WARN, NÃO um PASS biológico e NÃO é corrigida.
#   - Nenhuma análise inferencial casos vs controles nesta etapa.
#   - NENHUM colapso de tipos antes da análise de sensibilidade pré-definida.
#
# OUTPUTS
#   results/deconvolution/BisqueRNA/<dataset>/main24/top300/
#      proportions_wide.tsv, proportions_long.tsv,
#      proportion_summary_by_group.tsv, celltype_variability.tsv,
#      sample_diagnostics.tsv, input_audit.tsv, reference_support.tsv,
#      library_size_audit.tsv, parameters.tsv, genes_used.tsv,
#      bisque_result.rds, run_summary.tsv, status.tsv, 07_plot_status.tsv, plots/*.png
#
#   results/deconvolution/BisqueRNA/
#      07_run_plan.tsv, 07_run_summary.tsv, 07_proportions_long_all.tsv,
#      07_proportion_summary_all.tsv, 07_celltype_variability_all.tsv,
#      07_sample_diagnostics_all.tsv, 07_status.tsv,
#      07_output_manifest.tsv, 07_vs_MuSiC2.tsv (se disponível),
#      07_api_audit.tsv, 07_sessionInfo.txt, 07_plot_status_all.tsv, plots/*.png
#
# EXECUÇÃO (da raiz do repositório)
#   Rscript scripts/phase02_deconvolution/07_run_BisqueRNA.R
#
# INSTALAÇÃO (somente SE necessário, preferencialmente no seu ambiente renv):
#   install.packages('BisqueRNA')
#   BiocManager::install('Biobase')
#   install.packages(c('data.table','Matrix','ggplot2'))
#
# NÃO instala pacotes, não altera o renv.lock, não sobrescreve a etapa 06.
# ==============================================================================

# ==============================================================================
# 0. CONFIGURAÇÃO REPRODUTÍVEL
# ==============================================================================
MARKER_SETS_TO_RUN <- c("top300")   # top100/top500 apenas após QC principal
RUN_DATASETS <- c("ROSMAP", "GSE125583", "GSE64810", "GSE68719")
# Para diagnosticar UMA coorte sem alterar o script, use BISQUE_DATASET.
# Exemplos: BISQUE_DATASET=GSE64810 Rscript ... (Bash/WSL)
#           $env:BISQUE_DATASET="GSE64810"; Rscript ... (PowerShell)
RUN_DATASET_FILTER <- trimws(Sys.getenv("BISQUE_DATASET",unset=""))
if (nzchar(RUN_DATASET_FILTER)) {
  selected <- trimws(strsplit(RUN_DATASET_FILTER,",",fixed=TRUE)[[1]])
  if (!length(selected) || any(!selected %in% RUN_DATASETS)) {
    stop("BISQUE_DATASET inválido. Valores aceitos: ROSMAP, GSE125583, GSE64810, GSE68719.")
  }
  RUN_DATASETS <- unique(selected)
}
RUN_ROSMAP_21_DONOR_SENSITIVITY <- FALSE # apenas se decidido posteriormente

EXPECTED_CELL_TYPES <- c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")
EXPECTED_MAIN_NUCLEI <- 35140L
EXPECTED_MAIN_DONORS <- 24L
EXPECTED_REFERENCE_GENES <- 17394L
EXPECTED_MARKER_UNIONS <- c(top100=793L, top300=2242L, top500=3495L)

# Alertas operacionais, definidos antes de inspecionar os resultados do 07.
ZERO_TOL <- 1e-10
NEAR_ZERO_TOL <- 1e-4             # <0,01% de proporção
PROP_SUM_TOL <- 1e-3
NEGATIVE_TOL <- 1e-8
LOW_SD_PP <- 0.20                 # desvio-padrão inferior a 0,20 pp
LOW_RANGE_PP <- 0.50              # amplitude inferior a 0,50 pp
NEAR_UNIFORM_DISTANCE_PP <- 2.0  # RMSE médio em relação a 12,5%, em pp
UNIFORM_CELL_TOL_PP <- 2.0       # |proporção - 12,5%| <= 2 pp
MAX_NUCLEUS_MATRIX_GIB <- 2.0    # interromper se entrada densa for excessiva

BASE_SEED <- 20261021L

# v2: guardar o resultado do solver imediatamente, para não repetir um cálculo
# bem-sucedido caso haja erro de pós-processamento/figura. A reutilização
# exige que os tamanhos, nomes de genes/amostras e arquivos-fonte coincidam.
REUSE_VERIFIED_CHECKPOINTS <- TRUE

BULK_ROOT <- file.path("data", "processed", "deconvolution")
REF_ROOT <- file.path("data", "reference", "mathys2019", "phase02_prepared", "reference")
MARKER_ROOT <- file.path("data", "reference", "mathys2019", "phase02_prepared", "markers")
RESULT_ROOT <- file.path("results", "deconvolution", "BisqueRNA")

REF_COUNTS_FILE <- file.path(REF_ROOT, "mathys2019_pathology_free_geneid_counts_sparse.rds")
REF_META_FILE <- file.path(REF_ROOT, "mathys2019_pathology_free_cell_metadata_aligned.tsv")
STATUS_05_FILE <- file.path("results", "deconvolution", "bulk_preparation", "05_status.tsv")
MUSIC2_LONG_FILE <- file.path("results", "deconvolution", "MuSiC2", "06_proportions_long_all.tsv")

REGISTRY <- data.frame(
  dataset=c("ROSMAP", "GSE125583", "GSE64810", "GSE68719"),
  disease=c("AD", "AD", "HD", "PD"),
  cohort_role=c("principal", "validation", "principal", "principal"),
  stringsAsFactors=FALSE
)

PALETTE_CELL <- c(Ex="#3B5BDB", In="#7950F2", Ast="#0CA678",
                  Oli="#F08C00", Opc="#FCC419", Mic="#C2255C",
                  End="#15AABF", Per="#868E96")
PALETTE_GROUP <- c(Control="#94A3B8", AD="#7B2CBF",
                   HD="#277DA1", PD="#E75480")

# ==============================================================================
# 1. DEPENDÊNCIAS E FUNÇÕES BÁSICAS
# ==============================================================================
REQUIRED_PACKAGES <- c("BisqueRNA", "Biobase", "data.table", "Matrix", "ggplot2", "limSolve")
missing_packages <- REQUIRED_PACKAGES[!vapply(REQUIRED_PACKAGES, requireNamespace,
                                                logical(1), quietly=TRUE)]
if (length(missing_packages)) {
  stop("Pacotes ausentes: ", paste(missing_packages, collapse=", "),
       ". Instale-os no ambiente apropriado, sem atualizar o renv automaticamente.")
}
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Biobase)
  library(ggplot2)
})
dir.create(RESULT_ROOT, recursive=TRUE, showWarnings=FALSE)

say <- function(...) cat(sprintf(...), "\n", sep="")

# Registro de estágio por execução, útil caso o tryCatch suprima o traceback.
mark_stage <- function(run_dir, stage) {
  dir.create(run_dir,recursive=TRUE,showWarnings=FALSE)
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stage, sep="\t")
  cat(line,"\n",file=file.path(run_dir,"07_execution_stages.tsv"),append=TRUE)
  say("[07 diagnóstico] %s",stage)
}

file_signature <- function(path) {
  inf <- file.info(path)
  if (nrow(inf)!=1L || is.na(inf$size)) fail("Arquivo de entrada inacessível: ",path)
  list(path=normalizePath(path,winslash="/",mustWork=TRUE),
       size=as.numeric(inf$size),
       modified=as.numeric(inf$mtime))
}

fit_fingerprint <- function(inputs) {
  list(dataset=inputs$dataset,marker_set=inputs$marker_set,
       sample_ids=as.character(inputs$meta$sample_id),
       gene_ids=as.character(inputs$genes),
       cell_types=EXPECTED_CELL_TYPES,reference_nuclei=ncol(ref_counts),
       bisque_version=as.character(utils::packageVersion("BisqueRNA")),
       source_files=lapply(c(REF_COUNTS_FILE,REF_META_FILE,inputs$count_file,
                             inputs$metadata_file,inputs$lib_file,
                             inputs$marker_file),file_signature))
}
fail <- function(...) stop(paste0(...), call.=FALSE)
clean_id <- function(x) trimws(as.character(x))
write_tsv <- function(dt, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  data.table::fwrite(dt, path, sep="\t", na="NA", quote=FALSE)
}
make_status <- function(ds, ref, set, check, status, observed, expected="", note="") {
  data.table(dataset=ds, reference_version=ref, marker_set=set,
             check=check, status=status, observed=as.character(observed),
             expected=as.character(expected), note=as.character(note))
}
check_ids <- function(v, label) {
  if (is.null(v) || !length(v) || anyNA(v) || any(!nzchar(clean_id(v))) ||
      anyDuplicated(v)) fail(label, ": identificadores ausentes/duplicados.")
  invisible(TRUE)
}
check_counts <- function(x, label) {
  if (is.null(rownames(x)) || is.null(colnames(x))) fail(label, ": dimnames ausentes.")
  check_ids(rownames(x), paste0(label," genes"))
  check_ids(colnames(x), paste0(label," amostras/núcleos"))
  vals <- if (inherits(x,"sparseMatrix")) x@x else as.numeric(x)
  if (anyNA(vals) || any(!is.finite(vals)) || any(vals < 0) ||
      any(abs(vals-round(vals)) > 1e-8)) {
    fail(label, ": counts ausentes, fracionários, não finitos ou negativos.")
  }
  invisible(TRUE)
}

# Verifica que old.cpm=FALSE normaliza pelo total, inclusive o resíduo.
# Esta identidade é usada tanto como unidade matemática quanto nos dados reais.
verify_remainder_identity <- function(selected, total, label) {
  if (length(total) != ncol(selected)) fail(label, ": library size desalinhada.")
  if (any(!is.finite(total)) || any(total <= 0)) fail(label, ": library size inválida.")
  selected_sum <- colSums(selected)
  remain <- as.numeric(total - selected_sum)
  if (any(remain < -1e-6)) fail(label, ": soma dos marcadores excede library size.")
  remain[abs(remain) < 1e-8] <- 0
  recon <- selected_sum + remain
  if (max(abs(recon - total)) > 1e-5) fail(label, ": falha conservação de library size.")
  list(remainder=remain,
       selected_sum=as.numeric(selected_sum),
       maximum_abs_library_error=max(abs(recon-total)),
       min_remainder=min(remain),
       max_remainder=max(remain))
}

# Mini-teste para evitar eventual alteração inadvertida do cálculo de CPM.
mini <- matrix(c(2,3,4,5,6,7), nrow=2,
               dimnames=list(c("101","102"), c("s1","s2","s3")))
mini_total <- c(20,30,40)
mini_qc <- verify_remainder_identity(mini, mini_total, "Mini-teste CPM")
mini_full <- rbind(mini, "__LIBRARY_REMAINDER__"=mini_qc$remainder)
mini_cpm <- sweep(mini_full, 2, colSums(mini_full), "/")*1e6
mini_expected <- sweep(mini, 2, mini_total, "/")*1e6
if (!isTRUE(all.equal(unname(mini_cpm[1:2,,drop=FALSE]),
                      unname(mini_expected), tolerance=1e-10))) {
  fail("Falha no teste unitário: CPM não corresponde ao denominador completo.")
}
# Testar também o próprio normalizador instalado do BisqueRNA, não apenas
# nossa álgebra, para detectar uma mudança futura da implementação.
if (!exists("CountsToCPM", envir=asNamespace("BisqueRNA"), inherits=FALSE)) {
  fail("A versão instalada do BisqueRNA não possui CountsToCPM; API interna mudou.")
}
mini_eset <- Biobase::ExpressionSet(assayData=mini_full)
bisque_cpm_fn <- get("CountsToCPM", envir=asNamespace("BisqueRNA"))
mini_real_cpm <- Biobase::exprs(bisque_cpm_fn(mini_eset))[1:2,,drop=FALSE]
if (!isTRUE(all.equal(unname(mini_real_cpm),unname(mini_expected),
                      tolerance=1e-10))) {
  fail("A normalização CPM do BisqueRNA instalado falhou no teste de biblioteca.")
}
rm(mini,mini_full,mini_cpm,mini_expected,mini_eset,mini_real_cpm)

# ==============================================================================
# 2. CHECAGENS DE API E ARQUIVOS DA FASE 02
# ==============================================================================
fun <- get("ReferenceBasedDecomposition", envir=asNamespace("BisqueRNA"))
expected_formals <- c("bulk.eset", "sc.eset", "markers", "cell.types",
                      "subject.names", "use.overlap", "verbose", "old.cpm")
missing_formals <- setdiff(expected_formals, names(formals(fun)))
if (length(missing_formals)) fail("API BisqueRNA incompatível. Argumentos ausentes: ",
                                 paste(missing_formals, collapse=", "))
api_audit <- data.table(
  package="BisqueRNA", version=as.character(packageVersion("BisqueRNA")),
  entrypoint="ReferenceBasedDecomposition",
  use_overlap=FALSE, old_cpm=FALSE,
  marker_selection="genes topK explícitos + linha técnica de library size",
  output="bulk.props (tipos x amostras)",
  all_formals=paste(names(formals(fun)),collapse=";"))
write_tsv(api_audit, file.path(RESULT_ROOT,"07_api_audit.tsv"))

if (!file.exists(STATUS_05_FILE)) fail("Execute a etapa 05: ", STATUS_05_FILE)
status05 <- fread(STATUS_05_FILE)
if (any(status05$status=="FAIL",na.rm=TRUE)) fail("A etapa 05 possui FAIL; revise antes do 07.")
if (!file.exists(REF_COUNTS_FILE) || !file.exists(REF_META_FILE)) {
  fail("A referência Mathys da etapa 03 não foi encontrada em ", REF_ROOT)
}
if (length(setdiff(RUN_DATASETS, REGISTRY$dataset))) fail("Dataset desconhecido em RUN_DATASETS.")
if (length(setdiff(MARKER_SETS_TO_RUN, names(EXPECTED_MARKER_UNIONS)))) {
  fail("MARKER_SETS_TO_RUN deve conter apenas top100/top300/top500.")
}
if (RUN_ROSMAP_21_DONOR_SENSITIVITY) {
  fail("Sensibilidade ROSMAP 21 doadores ainda não implementada neste 07; ",
       "não ativar antes da decisão metodológica e da auditoria dessa referência.")
}

say("============================================================")
say("07_run_BisqueRNA.R — Referência Mathys pathology-free / 8 tipos")
say("============================================================")
say("BisqueRNA %s | old.cpm=FALSE | use.overlap=FALSE",
    as.character(packageVersion("BisqueRNA")))

# Metadados e counts nucleares originais: NÃO usar os 181 pseudobulks como células.
# Na transformação sem pares, Bisque calcula sc.props contando núcleos por
# tipo e doador, e precisa dessas observações individuais.
ref_counts <- readRDS(REF_COUNTS_FILE)
if (!inherits(ref_counts,"dgCMatrix")) ref_counts <- as(ref_counts,"dgCMatrix")
check_counts(ref_counts,"Mathys reference")
if (nrow(ref_counts)!=EXPECTED_REFERENCE_GENES ||
    ncol(ref_counts)!=EXPECTED_MAIN_NUCLEI) {
  fail("Dimensões Mathys inesperadas: ",nrow(ref_counts)," x ",ncol(ref_counts))
}
ref_meta <- fread(REF_META_FILE)
required_meta <- c("TAG","donor_id","cell_type")
if (length(setdiff(required_meta,names(ref_meta)))) fail("Metadata Mathys sem TAG/donor_id/cell_type.")
ref_meta[, TAG:=clean_id(TAG)]
ref_meta[, donor_id:=clean_id(donor_id)]
ref_meta[, cell_type:=clean_id(cell_type)]
check_ids(ref_meta$TAG,"TAG da referência")
idx_meta <- match(colnames(ref_counts),ref_meta$TAG)
if (anyNA(idx_meta)) fail("Existem núcleos Mathys sem metadados.")
ref_meta <- ref_meta[idx_meta]
if (!identical(colnames(ref_counts),ref_meta$TAG)) fail("Referência fora de ordem.")
if (uniqueN(ref_meta$donor_id)!=EXPECTED_MAIN_DONORS ||
    !setequal(unique(ref_meta$cell_type), EXPECTED_CELL_TYPES)) {
  fail("Números de doadores ou tipos celulares diferentes do esperado.")
}
ref_libraries <- as.numeric(Matrix::colSums(ref_counts))
if (any(ref_libraries<=0)) fail("Existem núcleos com library size zero na referência.")

reference_support <- ref_meta[,.(n_nuclei=.N,n_donors=uniqueN(donor_id)),by=cell_type]
reference_support[, cell_type:=factor(cell_type,levels=EXPECTED_CELL_TYPES)]
setorder(reference_support,cell_type)
reference_support[,cell_type:=as.character(cell_type)]

# ==============================================================================
# 3. LEITURA CONTROLADA DE UMA EXECUÇÃO
# ==============================================================================
read_run_inputs <- function(dataset, marker_set) {
  ds_row <- REGISTRY[REGISTRY$dataset==dataset,,drop=FALSE]
  disease <- ds_row$disease[[1]]
  role <- ds_row$cohort_role[[1]]

  bulk_file <- file.path(BULK_ROOT,dataset,paste0("bulk_counts_",marker_set,".rds"))
  bulk_meta_file <- file.path(BULK_ROOT,dataset,"metadata_deconvolution.tsv")
  lib_file <- file.path(BULK_ROOT,dataset,"library_sizes.tsv")
  marker_file <- file.path(MARKER_ROOT,paste0("mathys2019_",marker_set,"_markers.tsv"))
  needed <- c(bulk_file,bulk_meta_file,lib_file,marker_file)
  if (any(!file.exists(needed))) fail("Arquivos ausentes para ",dataset,": ",
                                    paste(needed[!file.exists(needed)],collapse="; "))

  bulk <- as.matrix(readRDS(bulk_file))
  storage.mode(bulk) <- "double"
  check_counts(bulk,paste0(dataset," bulk ",marker_set))

  meta <- fread(bulk_meta_file)
  for (x in c("sample_id","group","disease","cohort_role")) {
    if (!x %in% names(meta)) fail(dataset,": coluna ausente no metadata: ",x)
  }
  meta[,sample_id:=clean_id(sample_id)]
  check_ids(meta$sample_id,paste0(dataset," metadata"))
  if (!setequal(meta$sample_id,colnames(bulk))) fail(dataset,": bulk/metadata discordantes.")
  meta <- meta[match(colnames(bulk),meta$sample_id)]
  if (!identical(meta$sample_id,colnames(bulk))) fail(dataset,": falha na ordenação de amostras.")
  if (!setequal(unique(meta$group),c("Control",disease))) {
    fail(dataset,": grupos esperados Control e ",disease," não encontrados.")
  }
  if (any(meta$disease!=disease) || any(meta$cohort_role!=role)) {
    fail(dataset,": doença/papel diferentes do registro do estudo.")
  }

  libs <- fread(lib_file)
  if (!all(c("sample_id","library_size_full_input") %in% names(libs))) {
    fail(dataset,": arquivo library_sizes.tsv sem colunas esperadas.")
  }
  libs[,sample_id:=clean_id(sample_id)]
  check_ids(libs$sample_id,paste0(dataset," library sizes"))
  lib_idx <- match(meta$sample_id,libs$sample_id)
  if (anyNA(lib_idx)) fail(dataset,": library sizes não cobrem todas as amostras.")
  bulk_libraries <- as.numeric(libs$library_size_full_input[lib_idx])
  if (any(!is.finite(bulk_libraries)) || any(bulk_libraries<=0)) {
    fail(dataset,": library sizes inválidas.")
  }

  marks <- fread(marker_file)
  if (!all(c("gene_id","cell_type","marker_rank") %in% names(marks))) {
    fail("Marcadores ",marker_set,": colunas insuficientes.")
  }
  marks[,gene_id:=clean_id(gene_id)]
  if (!setequal(unique(marks$cell_type),EXPECTED_CELL_TYPES)) fail("Marcadores sem oito tipos.")
  k <- as.integer(sub("top","",marker_set))
  if (any(marks[, .N, by=cell_type]$N != k)) fail("Marcadores por tipo não somam topK.")
  marker_union <- unique(marks$gene_id)
  if (length(marker_union) != EXPECTED_MARKER_UNIONS[[marker_set]]) {
    fail("União ",marker_set," difere da referência 04: ",length(marker_union))
  }
  if (!all(rownames(bulk) %in% marker_union)) fail(dataset,": bulk contém genes fora de ",marker_set)

  common <- rownames(bulk)[rownames(bulk) %in% rownames(ref_counts)]
  if (!length(common)) fail(dataset,": sem genes em comum com Mathys.")
  bulk <- bulk[common,,drop=FALSE]
  sc <- ref_counts[common,,drop=FALSE]
  # Sem criar zeros para genes ausentes: remover apenas linhas totalmente zeradas.
  present <- rowSums(bulk)>0 & Matrix::rowSums(sc)>0
  used <- common[present]
  if (length(used)<500L) fail(dataset,": menos de 500 genes com sinal em ambas matrizes.")
  bulk <- bulk[used,,drop=FALSE]
  sc <- sc[used,,drop=FALSE]

  list(dataset=dataset,disease=disease,role=role,marker_set=marker_set,
       bulk=bulk,sc=sc,meta=meta,bulk_lib=bulk_libraries,
       genes_initial=common,genes=used,genes_removed=common[!present],
       marker_union=marker_union, count_file=bulk_file,
       metadata_file=bulk_meta_file,lib_file=lib_file,marker_file=marker_file)
}

# ==============================================================================
# 4. PREPARAR ExpressionSets COM CORREÇÃO EXPLÍCITA DA LIBRARY SIZE
# ==============================================================================
build_bisque_esets <- function(inputs) {
  # Duas cópias densas podem existir temporariamente no Bisque (biobase).
  gib <- (length(inputs$genes)+1) * ncol(ref_counts) * 8 / (1024^3)
  if (gib>MAX_NUCLEUS_MATRIX_GIB) {
    fail("Matriz nuclear densa estimada em ",round(gib,2)," GiB (limite ",
         MAX_NUCLEUS_MATRIX_GIB,"). Restrinja MARKER_SETS_TO_RUN; não faça ",
         "conversão densa de todos os genes sem planejar memória.")
  }

  sc_selected <- as.matrix(inputs$sc)
  storage.mode(sc_selected) <- "double"
  sc_remainders <- verify_remainder_identity(sc_selected,ref_libraries,
                                             paste0(inputs$dataset,"/sc"))
  bulk_remainders <- verify_remainder_identity(inputs$bulk,inputs$bulk_lib,
                                               paste0(inputs$dataset,"/bulk"))
  remainder_feature <- "__LIBRARY_REMAINDER__"
  if (remainder_feature %in% inputs$genes) fail("ID técnico colide com um gene real.")

  sc_with_remainder <- rbind(sc_selected,
                             "__LIBRARY_REMAINDER__"=sc_remainders$remainder)
  bulk_with_remainder <- rbind(inputs$bulk,
                               "__LIBRARY_REMAINDER__"=bulk_remainders$remainder)
  colnames(sc_with_remainder) <- colnames(ref_counts)
  colnames(bulk_with_remainder) <- colnames(inputs$bulk)
  if (max(abs(colSums(sc_with_remainder)-ref_libraries)) > 1e-5 ||
      max(abs(colSums(bulk_with_remainder)-inputs$bulk_lib)) > 1e-5) {
    fail(inputs$dataset,": não foi possível restaurar denominadores de biblioteca.")
  }
  # Check de CPM após preparar o mesmo tipo de entrada utilizado pelo pacote.
  # Testado em genes sentinela/amostras, sem densificar novamente toda matriz.
  gix <- seq_len(min(8L,length(inputs$genes)))
  six <- seq_len(min(5L,ncol(sc_with_remainder)))
  bix <- seq_len(min(5L,ncol(bulk_with_remainder)))
  actual_sc_cpm <- sweep(sc_with_remainder[gix,six,drop=FALSE],2,
                         colSums(sc_with_remainder[,six,drop=FALSE]),"/")*1e6
  expected_sc_cpm <- sweep(sc_selected[gix,six,drop=FALSE],2,
                           ref_libraries[six],"/")*1e6
  actual_bulk_cpm <- sweep(bulk_with_remainder[gix,bix,drop=FALSE],2,
                           colSums(bulk_with_remainder[,bix,drop=FALSE]),"/")*1e6
  expected_bulk_cpm <- sweep(inputs$bulk[gix,bix,drop=FALSE],2,
                             inputs$bulk_lib[bix],"/")*1e6
  if (max(abs(actual_sc_cpm-expected_sc_cpm))>1e-7 ||
      max(abs(actual_bulk_cpm-expected_bulk_cpm))>1e-7) {
    fail(inputs$dataset,": CPM de teste não usa library sizes completas.")
  }

  sc_pheno <- data.frame(cellType=factor(ref_meta$cell_type,levels=EXPECTED_CELL_TYPES),
                         SubjectName=factor(ref_meta$donor_id),
                         row.names=ref_meta$TAG, check.names=FALSE)
  bulk_pheno <- data.frame(group=inputs$meta$group,
                           row.names=inputs$meta$sample_id, check.names=FALSE)
  sc_eset <- Biobase::ExpressionSet(assayData=sc_with_remainder,
                                    phenoData=Biobase::AnnotatedDataFrame(sc_pheno))
  bulk_eset <- Biobase::ExpressionSet(assayData=bulk_with_remainder,
                                      phenoData=Biobase::AnnotatedDataFrame(bulk_pheno))
  if (!identical(Biobase::sampleNames(sc_eset),ref_meta$TAG) ||
      !identical(Biobase::sampleNames(bulk_eset),inputs$meta$sample_id)) {
    fail(inputs$dataset,": ExpressionSets desalinhados.")
  }

  library_audit <- data.table(
    dataset=inputs$dataset, marker_set=inputs$marker_set,
    source=c("Mathys nuclei","bulk samples"),
    n_profiles=c(ncol(sc_selected),ncol(inputs$bulk)),
    library_min=c(min(ref_libraries),min(inputs$bulk_lib)),
    library_median=c(median(ref_libraries),median(inputs$bulk_lib)),
    library_max=c(max(ref_libraries),max(inputs$bulk_lib)),
    remainder_min=c(sc_remainders$min_remainder,bulk_remainders$min_remainder),
    remainder_max=c(sc_remainders$max_remainder,bulk_remainders$max_remainder),
    maximum_abs_library_error=c(sc_remainders$maximum_abs_library_error,
                                bulk_remainders$maximum_abs_library_error),
    cpm_test_passed=TRUE)

  # Libera matrizes intermediárias quando a função retornar.
  list(sc_eset=sc_eset,bulk_eset=bulk_eset,library_audit=library_audit,
       dense_gib=as.numeric(gib))
}

# ==============================================================================
# 5. VALIDAR AS PROPORÇÕES E DETECTAR SAÍDAS DEGENERADAS
# ==============================================================================
extract_proportions <- function(result, inputs) {
  if (is.null(result$bulk.props)) fail(inputs$dataset,": Bisque sem bulk.props.")
  prop <- as.matrix(result$bulk.props)
  storage.mode(prop) <- "double"
  if (is.null(rownames(prop)) || is.null(colnames(prop))) fail("Bisque sem dimnames.")
  if (!setequal(rownames(prop),EXPECTED_CELL_TYPES)) {
    fail(inputs$dataset,": o Bisque não retornou os oito tipos esperados.")
  }
  if (!setequal(colnames(prop),inputs$meta$sample_id)) {
    fail(inputs$dataset,": Bisque retornou amostras diferentes do bulk.")
  }
  prop <- prop[EXPECTED_CELL_TYPES,inputs$meta$sample_id,drop=FALSE]
  if (anyNA(prop) || any(!is.finite(prop))) fail(inputs$dataset,": proporções NA/Inf.")
  if (min(prop) < -NEGATIVE_TOL) fail(inputs$dataset,": proporções negativas.")
  if (max(abs(colSums(prop)-1))>PROP_SUM_TOL) {
    fail(inputs$dataset,": soma das proporções não é ~1. Não renormalizar artificialmente.")
  }
  if (!is.null(result$genes.used) &&
      "__LIBRARY_REMAINDER__" %in% as.character(result$genes.used)) {
    fail(inputs$dataset,": linha técnica entrou na deconvolução! ABORTADO.")
  }
  if (!is.null(result$genes.used) &&
      any(!result$genes.used %in% inputs$genes)) {
    fail(inputs$dataset,": Bisque utilizou gene fora da interseção aprovada.")
  }
  # Retornar amostras x tipos, sem truncar ou renormalizar valores.
  t(prop)
}

summarize_bisque <- function(proportions, meta, dataset, disease, role, set) {
  wide <- data.table(sample_id=rownames(proportions), group=meta$group)
  wide <- cbind(wide,as.data.table(proportions))
  long <- data.table::melt(wide,id.vars=c("sample_id","group"),
                            variable.name="cell_type",value.name="proportion",
                            variable.factor=FALSE)
  long[,`:=`(dataset=dataset,disease=disease,cohort_role=role,
             reference_version="main24",marker_set=set)]
  setcolorder(long,c("dataset","disease","cohort_role","reference_version",
                     "marker_set","sample_id","group","cell_type","proportion"))

  summary <- long[,.(n_samples=.N,
                      mean_proportion=mean(proportion),
                      median_proportion=median(proportion),
                      sd_proportion=sd(proportion),
                      min_proportion=min(proportion),
                      max_proportion=max(proportion)),
                  by=.(dataset,disease,cohort_role,reference_version,marker_set,group,cell_type)]

  # Variabilidade individual em CONTROLES, CASOS e na COORTE COMPLETA.
  all_dt <- copy(long)
  all_dt[,group:="ALL"]
  long_for_qc <- rbindlist(list(long,all_dt),use.names=TRUE)
  variability <- long_for_qc[,.(
    n_samples=.N,
    mean_proportion=mean(proportion),
    median_proportion=median(proportion),
    sd_pp=if (.N>1) 100*sd(proportion) else NA_real_,
    range_pp=100*(max(proportion)-min(proportion)),
    n_zero=sum(proportion<=ZERO_TOL),
    fraction_zero=mean(proportion<=ZERO_TOL),
    fraction_near_zero=mean(proportion<NEAR_ZERO_TOL),
    n_distinct_6dec=uniqueN(round(proportion,6))),
    by=.(dataset,reference_version,marker_set,group,cell_type)]

  uniform <- 1/length(EXPECTED_CELL_TYPES)
  uniform_rmse <- sqrt(rowMeans((proportions-uniform)^2))*100
  all_close_uniform <- rowSums(abs(proportions-uniform)*100 <= UNIFORM_CELL_TOL_PP)==
    length(EXPECTED_CELL_TYPES)
  sample_diag <- data.table(
    dataset=dataset,reference_version="main24",marker_set=set,
    sample_id=rownames(proportions),group=meta$group,
    proportion_sum=rowSums(proportions),
    euclidean_distance_to_uniform_pp=sqrt(rowSums((proportions-uniform)^2))*100,
    rmse_to_uniform_pp=uniform_rmse,
    all_eight_within_2pp_uniform=all_close_uniform,
    n_exact_zeros=rowSums(proportions<=ZERO_TOL))

  all_var <- variability[group=="ALL"]
  is_near_constant <- all(all_var$sd_pp<LOW_SD_PP & all_var$range_pp<LOW_RANGE_PP)
  is_near_uniform <- median(uniform_rmse)<NEAR_UNIFORM_DISTANCE_PP &&
    mean(all_close_uniform)>=0.90
  is_uniform_degenerate <- is_near_constant && is_near_uniform
  is_constant_degenerate <- is_near_constant
  qc <- list(
    degenerate_uniform=is_uniform_degenerate,
    degenerate_constant=is_constant_degenerate,
    low_variance_celltypes=sum(all_var$sd_pp<LOW_SD_PP),
    n_unique_vectors=uniqueN(as.data.table(round(proportions,6))),
    fraction_samples_uniform=mean(all_close_uniform),
    median_uniform_rmse_pp=median(uniform_rmse),
    median_celltype_sd_pp=median(all_var$sd_pp))

  list(wide=wide,long=long,summary=summary,variability=variability,
       sample_diag=sample_diag, qc=qc)
}

# ==============================================================================
# 6. FIGURAS COLORIDAS — INTEGRAÇÃO À ETAPA 07 (v4: pilhas preservadas)
# ==============================================================================
base_theme <- function() theme_minimal(base_size=12)+
  theme(plot.title=element_text(face="bold",size=14),
        axis.title=element_text(size=11),
        panel.grid.minor=element_blank(),
        legend.position="bottom",strip.text=element_text(face="bold"))

# Salva uma figura e verifica que o PNG existe e não está vazio.
# Um erro isolado de arquivo gráfico NÃO invalida nem recalcule o ajuste Bisque.
save_plot_07 <- function(plot_object, path, width, height, dpi=170L) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  error_message <- tryCatch({
    ggplot2::ggsave(filename=path, plot=plot_object,
                    width=width, height=height, dpi=dpi,
                    bg="white")
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size<=0) {
      stop("Imagem não foi criada ou tem tamanho zero.")
    }
    NA_character_
  }, error=function(e) conditionMessage(e))
  data.table(file=path, status=if (is.na(error_message)) "PASS" else "WARN",
             error=error_message)
}

plot_one_dataset <- function(out, disease, dataset, plot_dir) {
  dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)
  d <- copy(out$long)
  d[,pct:=proportion*100]
  d[,cell_type:=factor(cell_type,levels=EXPECTED_CELL_TYPES)]
  # ATENÇÃO: "disease" também é o nome de uma coluna de d! Dentro de
  # data.table, usar levels=c("Control", disease) recuperaria a COLUNA inteira,
  # repetindo HD/PD/AD centenas de vezes e disparando:
  # "factor level [3] is duplicated".
  # Construir níveis FORA de [.data.table e atribuir com data.table::set.
  disease_label <- as.character(disease)[1L]
  group_levels <- c("Control", disease_label)
  if (anyDuplicated(group_levels) || !setequal(unique(as.character(d[["group"]])), group_levels)) {
    fail(dataset, ": grupos inesperados para os gráficos: ",
         paste(unique(as.character(d[["group"]])), collapse=", "))
  }
  data.table::set(d, j="group",
                  value=factor(as.character(d[["group"]]), levels=group_levels))

  # v5 — Mais espaço vertical, margens internas/externas e pontos sem recorte.
  # Escalas Y continuam livres: NÃO comparar alturas absolutas de facetas distintas.
  # "expand" adiciona respiro no topo/rodapé, mesmo para valores extremos.
  g1 <- ggplot(d,aes(x=group,y=pct,fill=group))+
    geom_boxplot(width=.60,outlier.shape=NA,alpha=.78,linewidth=.45)+
    geom_point(aes(color=group),
               position=position_jitter(width=.12,height=0,seed=BASE_SEED),
               alpha=.38,size=.80)+
    facet_wrap(~cell_type,scales="free_y",ncol=4)+
    scale_fill_manual(values=PALETTE_GROUP,drop=FALSE)+
    scale_color_manual(values=PALETTE_GROUP,drop=FALSE)+
    scale_x_discrete(expand=expansion(add=.42))+
    scale_y_continuous(expand=expansion(mult=c(.14,.23)))+
    coord_cartesian(clip="off")+
    labs(title=paste0("BisqueRNA • ",dataset," • proporções individuais"),
         subtitle="Controle × doença • cada ponto representa uma amostra • escalas Y independentes",
         x=NULL,y="Contribuição estimada (%)",fill="Grupo",color="Grupo")+
    base_theme()+
    theme(axis.text.x=element_text(angle=15,hjust=1,size=10),
          strip.text=element_text(face="bold",size=12),
          panel.spacing.x=grid::unit(1.4,"lines"),
          panel.spacing.y=grid::unit(1.5,"lines"),
          plot.margin=margin(t=24,r=32,b=28,l=32,unit="pt"),
          legend.margin=margin(t=10,unit="pt"),
          plot.subtitle=element_text(margin=margin(b=14,unit="pt")))
  per_plot_status <- list()
  per_plot_status[[length(per_plot_status)+1L]] <- save_plot_07(
    g1, file.path(plot_dir,"07_boxplots_by_celltype.png"),16.5,10.5,185)

  g2dt <- d[,.(mean_pct=mean(pct)),by=.(group,cell_type)]
  # v4 — Auditoria obrigatória antes de gerar o gráfico empilhado.
  # A soma das médias de oito proporções fechadas deve ser ~100% por grupo.
  # Nesta etapa, NÃO alterar as estimativas ou renormalizá-las para desenhar.
  g2tot <- g2dt[,.(sum_mean_pct=sum(mean_pct),n_types=.N),by=.(group)]
  if (any(g2tot$n_types != length(EXPECTED_CELL_TYPES)) ||
      any(!is.finite(g2tot$sum_mean_pct)) ||
      any(abs(g2tot$sum_mean_pct-100)>1e-6)) {
    fail(dataset, ": gráfico de composição — médias não somam 100%. ",
         paste(g2tot$group, round(g2tot$sum_mean_pct,6), collapse="; "))
  }
  write_tsv(g2tot, file.path(dirname(plot_dir),"07_plot_group_totals.tsv"))
  # IMPORTANTE: scale_y_continuous(limits=c(0,100)) aplica censura antes do
  # desenho. Por erros de arredondamento (~1e-13), a pilha pode superar
  # ligeiramente 100%; ggplot2 então pode eliminar INTEIRO um segmento (Ex),
  # causando barras artificiais de ~47–55%. O zoom com coord_cartesian NÃO
  # remove os dados antes de empilhar e conserva os oito segmentos.
  g2 <- ggplot(g2dt,aes(x=group,y=mean_pct,fill=cell_type))+
    geom_col(width=.62,color="white",linewidth=.25)+
    scale_fill_manual(values=PALETTE_CELL,drop=FALSE)+
    scale_y_continuous(breaks=seq(0,100,25),
                       expand=expansion(mult=c(0,0.02)))+
    coord_cartesian(ylim=c(0,100))+
    labs(title=paste0("BisqueRNA • ",dataset," • composição média"),
         subtitle="Médias descritivas: não implicam abundância absoluta de células",
         x=NULL,y="Contribuição média estimada (%)",fill="Tipo celular")+base_theme()
  per_plot_status[[length(per_plot_status)+1L]] <- save_plot_07(
    g2, file.path(plot_dir,"07_mean_composition_stacked.png"),10,6,170)

  q <- copy(out$variability[group=="ALL"])
  q[,cell_type:=factor(cell_type,levels=EXPECTED_CELL_TYPES)]
  g3 <- ggplot(q,aes(x=cell_type,y=sd_pp,fill=cell_type))+
    geom_col(show.legend=FALSE,width=.72)+
    geom_hline(yintercept=LOW_SD_PP,linetype="dashed",color="#C92A2A")+
    scale_fill_manual(values=PALETTE_CELL)+
    labs(title=paste0("BisqueRNA • ",dataset," • variabilidade entre amostras"),
         subtitle=paste0("Linha tracejada = alerta técnico SD < ",LOW_SD_PP," ponto percentual"),
         x="Tipo celular",y="Desvio-padrão (pontos percentuais)")+base_theme()
  per_plot_status[[length(per_plot_status)+1L]] <- save_plot_07(
    g3, file.path(plot_dir,"07_variability_sd.png"),10,5.2,170)

  g4 <- ggplot(q,aes(x=cell_type,y=100*fraction_zero,fill=cell_type))+
    geom_col(show.legend=FALSE,width=.72)+
    scale_fill_manual(values=PALETTE_CELL)+
    scale_y_continuous(limits=c(0,100))+
    labs(title=paste0("BisqueRNA • ",dataset," • estimativas iguais a zero"),
         subtitle="Zero estimado não significa ausência biológica",
         x="Tipo celular",y="Amostras com proporção ~0 (%)")+base_theme()
  per_plot_status[[length(per_plot_status)+1L]] <- save_plot_07(
    g4, file.path(plot_dir,"07_zero_fraction.png"),10,5.2,170)
  rbindlist(per_plot_status, use.names=TRUE, fill=TRUE)
}

# ==============================================================================
# 7. EXECUTAR UMA COORTE E SALVAR SAÍDAS INDIVIDUAIS
# ==============================================================================
run_one <- function(dataset, marker_set, seed) {
  run_dir <- file.path(RESULT_ROOT,dataset,"main24",marker_set)
  dir.create(run_dir,recursive=TRUE,showWarnings=FALSE)
  plot_dir <- file.path(run_dir,"plots")
  # O log é reiniciado em cada execução; mantém histórico de estágios recentes.
  writeLines("timestamp\tstage",file.path(run_dir,"07_execution_stages.tsv"))
  mark_stage(run_dir,"START")
  say("\n------------------------------------------------------------")
  say("Bisque: %s | %s | main24",dataset,marker_set)

  inputs <- read_run_inputs(dataset,marker_set)
  mark_stage(run_dir,"Entradas auditadas")
  prepared <- build_bisque_esets(inputs)
  mark_stage(run_dir,"ExpressionSets preparados")
  audit <- data.table(metric=c("bulk_samples","reference_nuclei","reference_donors",
                                "genes_before_allzero_filter","genes_input","genes_removed_allzero",
                                "technical_remainder_feature","estimated_dense_sc_gib",
                                "bulk_library_size_source","reference_library_size_source"),
                      value=as.character(c(ncol(inputs$bulk),ncol(inputs$sc),EXPECTED_MAIN_DONORS,
                               length(inputs$genes_initial),length(inputs$genes),
                               length(inputs$genes_removed),"__LIBRARY_REMAINDER__",
                               round(prepared$dense_gib,3),
                               "05/library_sizes.tsv full input",
                               "03/Mathys nucleus-level full GeneID counts")))
  write_tsv(audit,file.path(run_dir,"input_audit.tsv"))
  write_tsv(prepared$library_audit,file.path(run_dir,"library_size_audit.tsv"))
  write_tsv(reference_support,file.path(run_dir,"reference_support.tsv"))
  write_tsv(data.table(gene_id=inputs$genes),file.path(run_dir,"genes_input.tsv"))
  write_tsv(data.table(gene_id=inputs$genes_removed),file.path(run_dir,"genes_removed_allzero.tsv"))
  pars <- data.table(parameter=c("dataset","reference_version","marker_set","seed",
                                  "function","use.overlap","old.cpm","cell.types",
                                  "subject.names","technical_feature","reference_nuclei"),
                     value=c(dataset,"main24",marker_set,as.character(seed),
                             "ReferenceBasedDecomposition","FALSE","FALSE","cellType",
                             "SubjectName","__LIBRARY_REMAINDER__",
                             as.character(ncol(inputs$sc))))
  write_tsv(pars,file.path(run_dir,"parameters.tsv"))

  set.seed(seed)
  start <- proc.time()[["elapsed"]]
  say("Genes: %d | amostras: %d | núcleos: %d",length(inputs$genes),
      ncol(inputs$bulk),ncol(inputs$sc))
  checkpoint_file <- file.path(run_dir,"07_raw_fit_checkpoint.rds")
  fingerprint <- fit_fingerprint(inputs)
  reuse <- FALSE
  result <- NULL
  if (REUSE_VERIFIED_CHECKPOINTS && file.exists(checkpoint_file)) {
    mark_stage(run_dir,"Verificando checkpoint anterior")
    cp <- tryCatch(readRDS(checkpoint_file),error=function(e) NULL)
    if (is.list(cp) && identical(cp$fingerprint,fingerprint) &&
        is.list(cp$result) && !is.null(cp$result$bulk.props)) {
      result <- cp$result
      reuse <- TRUE
      mark_stage(run_dir,"Checkpoint compatível recuperado; sem novo ajuste")
    } else {
      mark_stage(run_dir,"Checkpoint inválido/desatualizado: refazer ajuste")
    }
    rm(cp)
  }
  if (!reuse) {
    mark_stage(run_dir,"ENTRANDO em BisqueRNA::ReferenceBasedDecomposition")
    result <- BisqueRNA::ReferenceBasedDecomposition(
      bulk.eset=prepared$bulk_eset,
      sc.eset=prepared$sc_eset,
      markers=inputs$genes,
      cell.types="cellType",
      subject.names="SubjectName",
      use.overlap=FALSE,
      verbose=TRUE,
      old.cpm=FALSE
    )
    mark_stage(run_dir,"SAIU do solver BisqueRNA com sucesso")
    # Salvar ANTES de extração, auditorias, resumos e plots.
    saveRDS(list(fingerprint=fingerprint,result=result),checkpoint_file,
            compress=TRUE)
    mark_stage(run_dir,"Checkpoint bruto salvo")
  }
  elapsed <- proc.time()[["elapsed"]]-start
  mark_stage(run_dir,"Extraindo e validando proporções")
  p <- extract_proportions(result,inputs)
  mark_stage(run_dir,"Calculando resumos, variância e degenerescência")
  out <- summarize_bisque(p,inputs$meta,dataset,inputs$disease,inputs$role,marker_set)

  # v5 — Auditoria da ancoragem no perfil de referência Mathys.
  # Quando use.overlap=FALSE, a transformação semissupervisionada do Bisque
  # centraliza cada gene na média de sc.ref %*% sc.props (código original).
  # A comparação NÃO comprova erro, mas alerta para possível prior/calibração.
  sc_props_fit <- result$sc.props
  if (is.matrix(sc_props_fit) &&
      setequal(rownames(sc_props_fit),EXPECTED_CELL_TYPES) &&
      all(is.finite(sc_props_fit))) {
    prior_audit <- data.table(
      dataset=dataset, reference_version="main24", marker_set=marker_set,
      cell_type=EXPECTED_CELL_TYPES,
      mean_bulk_estimate=as.numeric(colMeans(p)[EXPECTED_CELL_TYPES]),
      mean_reference_donor_fraction=as.numeric(rowMeans(sc_props_fit[EXPECTED_CELL_TYPES,,drop=FALSE])),
      pooled_reference_nucleus_fraction=as.numeric(
        reference_support$n_nuclei[match(EXPECTED_CELL_TYPES,reference_support$cell_type)] / nrow(ref_meta)))
    prior_audit[,`:=`(
      delta_bulk_minus_donor_pp=100*(mean_bulk_estimate-mean_reference_donor_fraction),
      delta_bulk_minus_pooled_pp=100*(mean_bulk_estimate-pooled_reference_nucleus_fraction))]
    write_tsv(prior_audit,file.path(run_dir,"07_reference_anchor_audit.tsv"))
    if (max(abs(prior_audit$delta_bulk_minus_donor_pp))<1.0) {
      mark_stage(run_dir,"AVISO: médias bulk próximas das médias dos doadores da referência")
    }
  }

  # 'rnorm' é resíduo do ajuste, não medida normalizada de acurácia biológica.
  residuals <- if (!is.null(result$rnorm)) result$rnorm[inputs$meta$sample_id] else
    setNames(rep(NA_real_,nrow(inputs$meta)),inputs$meta$sample_id)
  out$sample_diag[,residual_norm:=as.numeric(residuals[sample_id])]

  genes_final <- if (is.null(result$genes.used)) character() else as.character(result$genes.used)
  summary <- data.table(dataset=dataset,disease=inputs$disease,cohort_role=inputs$role,
                        reference_version="main24",marker_set=marker_set,
                        reference_donors=EXPECTED_MAIN_DONORS,
                        reference_nuclei=ncol(inputs$sc),
                        genes_input=length(inputs$genes),
                        genes_used=length(genes_final),
                        bulk_samples=nrow(p),
                        controls=sum(inputs$meta$group=="Control"),
                        cases=sum(inputs$meta$group==inputs$disease),
                        max_proportion_sum_error=max(abs(rowSums(p)-1)),
                        min_proportion=min(p),max_proportion=max(p),
                        degenerate_uniform=out$qc$degenerate_uniform,
                        degenerate_constant=out$qc$degenerate_constant,
                        low_variance_celltypes=out$qc$low_variance_celltypes,
                        unique_vectors_6dec=out$qc$n_unique_vectors,
                        fraction_samples_uniform=out$qc$fraction_samples_uniform,
                        median_uniform_rmse_pp=out$qc$median_uniform_rmse_pp,
                        median_celltype_sd_pp=out$qc$median_celltype_sd_pp,
                        elapsed_seconds=elapsed)

  status <- rbindlist(list(
    make_status(dataset,"main24",marker_set,"execution","PASS","completed","completed"),
    make_status(dataset,"main24",marker_set,"eight_celltypes","PASS",ncol(p),8L),
    make_status(dataset,"main24",marker_set,"all_samples_returned","PASS",nrow(p),nrow(inputs$meta)),
    make_status(dataset,"main24",marker_set,"sum_to_one","PASS",summary$max_proportion_sum_error,PROP_SUM_TOL),
    make_status(dataset,"main24",marker_set,"no_technical_feature_used","PASS",
                !("__LIBRARY_REMAINDER__" %in% genes_final),TRUE),
    make_status(dataset,"main24",marker_set,"constant_vector_degeneracy",
                if (out$qc$degenerate_constant) "WARN" else "PASS",
                out$qc$degenerate_constant,FALSE,
                "Todas as oito proporções quase constantes entre amostras: revisar Bisque."),
    make_status(dataset,"main24",marker_set,"near_uniform_degeneracy",
                if (out$qc$degenerate_uniform) "WARN" else "PASS",
                out$qc$degenerate_uniform,FALSE,
                "Próximo de 12,5% por tipo e sem variabilidade: NÃO usar como evidência independente.")
  ),use.names=TRUE,fill=TRUE)

  write_tsv(out$wide,file.path(run_dir,"proportions_wide.tsv"))
  write_tsv(out$long,file.path(run_dir,"proportions_long.tsv"))
  write_tsv(out$summary,file.path(run_dir,"proportion_summary_by_group.tsv"))
  write_tsv(out$variability,file.path(run_dir,"celltype_variability.tsv"))
  write_tsv(out$sample_diag,file.path(run_dir,"sample_diagnostics.tsv"))
  write_tsv(data.table(gene_id=genes_final),file.path(run_dir,"genes_used.tsv"))
  write_tsv(summary,file.path(run_dir,"run_summary.tsv"))
  write_tsv(status,file.path(run_dir,"status.tsv"))
  saveRDS(result,file.path(run_dir,"bisque_result.rds"),compress=TRUE)
  mark_stage(run_dir,"Tabelas e objeto Bisque final salvos")
  mark_stage(run_dir,"Gerando figuras individuais")
  plot_status <- tryCatch(
    plot_one_dataset(out, inputs$disease, dataset, plot_dir),
    error=function(e) data.table(file="construcao_das_figuras",
                                 status="WARN", error=conditionMessage(e)))
  plot_status[,`:=`(dataset=dataset,reference_version="main24",marker_set=marker_set)]
  setcolorder(plot_status,c("dataset","reference_version","marker_set","file","status","error"))
  write_tsv(plot_status,file.path(run_dir,"07_plot_status.tsv"))
  plots_ok <- nrow(plot_status)==4L && all(plot_status$status=="PASS")
  status <- rbind(status,make_status(dataset,"main24",marker_set,
     "plot_generation",if (plots_ok) "PASS" else "WARN",
     paste0(sum(plot_status$status=="PASS"),"/4 figuras geradas"),"4/4"),fill=TRUE)
  write_tsv(status,file.path(run_dir,"status.tsv"))
  if (!plots_ok) {
    bad_plot <- plot_status[status!="PASS"]
    writeLines(paste(bad_plot$file,bad_plot$error,sep=" | "),
               file.path(run_dir,"07_plot_error.txt"))
    warning(dataset,": uma ou mais figuras falharam; resultados numéricos preservados. ",
            "Consulte 07_plot_status.tsv")
    mark_stage(run_dir,"WARN: algumas figuras falharam; resultados numéricos preservados")
  } else {
    if (file.exists(file.path(run_dir,"07_plot_error.txt")))
      unlink(file.path(run_dir,"07_plot_error.txt"))
    mark_stage(run_dir,"Quatro figuras individuais geradas")
  }
  if (file.exists(file.path(run_dir,"ERROR.tsv")))
    unlink(file.path(run_dir,"ERROR.tsv"))
  mark_stage(run_dir,"Concluído")

  say("Concluído: %s | %d genes finais | degenerescência uniforme=%s",
      dataset,length(genes_final),out$qc$degenerate_uniform)
  if (out$qc$degenerate_constant) {
    warning(dataset,": estimativas quase constantes detectadas; verifique o 07_status.tsv.")
  }
  # Devolver apenas tabelas; grandes ExpressionSets liberados ao sair da função.
  list(summary=summary,long=out$long,group_summary=out$summary,
       variability=out$variability,sample_diag=out$sample_diag,
       status=status,plot_status=plot_status,run_dir=run_dir)
}

# ==============================================================================
# 8. PLANO / EXECUÇÃO PROTEGIDA / CONSOLIDAÇÃO
# ==============================================================================
plan <- CJ(dataset=RUN_DATASETS,marker_set=MARKER_SETS_TO_RUN,unique=TRUE)
plan[,dataset_order:=match(dataset,RUN_DATASETS)]
plan[,marker_order:=match(marker_set,MARKER_SETS_TO_RUN)]
setorder(plan,dataset_order,marker_order)
plan[,c("dataset_order","marker_order"):=NULL]
plan[,reference_version:="main24"]
plan[,run_order:=seq_len(.N)]
write_tsv(plan,file.path(RESULT_ROOT,"07_run_plan.tsv"))
say("Coortes selecionadas nesta execução (%d): %s",length(RUN_DATASETS),
    paste(RUN_DATASETS,collapse=", "))
if (nzchar(RUN_DATASET_FILTER)) {
  say("ATENÇÃO: filtro BISQUE_DATASET ativo: %s",RUN_DATASET_FILTER)
}

runs <- list()
statuses <- list()
errors <- list()
for (i in seq_len(nrow(plan))) {
  row <- plan[i]
  ds <- row$dataset[[1]]
  ms <- row$marker_set[[1]]
  seed <- as.integer(BASE_SEED+i-1L)
  # Cada execução deve ser isolada: falhas não destroem os resultados de outros.
  run_dir <- file.path(RESULT_ROOT,ds,"main24",ms)
  x <- tryCatch(
    withCallingHandlers(
      run_one(ds,ms,seed),
      error=function(e) {
        # Handler executado ANTES do desempilhamento da pilha R.
        calls <- sys.calls()
        trace_lines <- vapply(calls,function(cl) {
          txt <- paste(deparse(cl,width.cutoff=110L),collapse=" ")
          substr(txt,1L,1000L)
        },character(1))
        trace_file <- file.path(run_dir,"07_error_trace.txt")
        dir.create(run_dir,recursive=TRUE,showWarnings=FALSE)
        stages_file <- file.path(run_dir,"07_execution_stages.tsv")
        last_stage <- if (file.exists(stages_file)) tail(readLines(stages_file,warn=FALSE),1L) else NA_character_
        writeLines(c(paste("dataset:",ds),paste("marker:",ms),
                     paste("last_stage:",last_stage),
                     paste("error:",conditionMessage(e)),
                     "--- STACK (chamadas mais internas por último) ---",
                     paste(seq_along(trace_lines),trace_lines,sep=": ")),
                   trace_file)
        message("[07] Traceback gravado em: ",trace_file)
      }
    ),error=function(e)e
  )
  if (inherits(x,"error")) {
    error_dt <- data.table(dataset=ds,reference_version="main24",marker_set=ms,
                           error=conditionMessage(x))
    write_tsv(error_dt,file.path(RESULT_ROOT,ds,"main24",ms,"ERROR.tsv"))
    errors[[length(errors)+1L]] <- error_dt
    statuses[[length(statuses)+1L]] <- make_status(ds,"main24",ms,"execution","FAIL",
                                            conditionMessage(x),"completed")
    message("FALHA em ",ds,"/",ms,": ",conditionMessage(x))
  } else {
    runs[[length(runs)+1L]] <- x
    statuses[[length(statuses)+1L]] <- x$status
  }
  invisible(gc(verbose=FALSE))
}

status_all <- rbindlist(statuses,use.names=TRUE,fill=TRUE)
write_tsv(status_all,file.path(RESULT_ROOT,"07_status.tsv"))

if (length(runs)) {
  run_summary_all <- rbindlist(lapply(runs,`[[`,"summary"),use.names=TRUE,fill=TRUE)
  prop_long_all <- rbindlist(lapply(runs,`[[`,"long"),use.names=TRUE,fill=TRUE)
  group_all <- rbindlist(lapply(runs,`[[`,"group_summary"),use.names=TRUE,fill=TRUE)
  var_all <- rbindlist(lapply(runs,`[[`,"variability"),use.names=TRUE,fill=TRUE)
  sample_all <- rbindlist(lapply(runs,`[[`,"sample_diag"),use.names=TRUE,fill=TRUE)
  write_tsv(run_summary_all,file.path(RESULT_ROOT,"07_run_summary.tsv"))
  write_tsv(prop_long_all,file.path(RESULT_ROOT,"07_proportions_long_all.tsv"))
  write_tsv(group_all,file.path(RESULT_ROOT,"07_proportion_summary_all.tsv"))
  write_tsv(var_all,file.path(RESULT_ROOT,"07_celltype_variability_all.tsv"))
  write_tsv(sample_all,file.path(RESULT_ROOT,"07_sample_diagnostics_all.tsv"))

  # Comparação DESCRITIVA com a etapa 06 (sem declarar concordância de validade).
  if (file.exists(MUSIC2_LONG_FILE)) {
    mus <- fread(MUSIC2_LONG_FILE)
    required <- c("dataset","reference_version","marker_set","sample_id",
                  "cell_type","proportion")
    if (all(required %in% names(mus))) {
      mus <- mus[reference_version=="main24" & marker_set %in% MARKER_SETS_TO_RUN,
                 .(dataset,marker_set,sample_id,cell_type,music2=proportion)]
      bis <- prop_long_all[,.(dataset,marker_set,sample_id,cell_type,bisque=proportion)]
      matched <- merge(bis,mus,by=c("dataset","marker_set","sample_id","cell_type"))
      if (nrow(matched)>0L) {
        comparison <- matched[,.(
          n_samples=.N,
          bisque_mean=mean(bisque),music2_mean=mean(music2),
          mean_abs_difference_pp=100*mean(abs(bisque-music2)),
          spearman=if (sd(bisque)>1e-12 && sd(music2)>1e-12)
            suppressWarnings(cor(bisque,music2,method="spearman")) else NA_real_),
          by=.(dataset,marker_set,cell_type)]
        comparison[,note:="Descritivo; discordância não determina qual método é correto."]
        write_tsv(comparison,file.path(RESULT_ROOT,"07_vs_MuSiC2.tsv"))
      }
    }
  }

  # Gráficos consolidados: porcentagem de zeros e variação por coorte/tipo.
  global_plot_dir <- file.path(RESULT_ROOT,"plots")
  dir.create(global_plot_dir,recursive=TRUE,showWarnings=FALSE)
  v <- copy(var_all[group=="ALL"])
  v[,cell_type:=factor(cell_type,levels=EXPECTED_CELL_TYPES)]
  p0 <- ggplot(v,aes(x=cell_type,y=100*fraction_zero,fill=cell_type))+
    geom_col(show.legend=FALSE)+facet_wrap(~dataset,ncol=2)+
    scale_fill_manual(values=PALETTE_CELL)+scale_y_continuous(limits=c(0,100))+
    labs(title="BisqueRNA • fração de estimativas zero por coorte",
         x="Tipo celular",y="Amostras com zero (%)")+base_theme()
  global_plot_status <- list()
  global_plot_status[[1L]] <- save_plot_07(
    p0,file.path(global_plot_dir,"07_zero_fraction_all.png"),12,8,180)
  pv <- ggplot(v,aes(x=cell_type,y=sd_pp,fill=cell_type))+
    geom_col(show.legend=FALSE)+facet_wrap(~dataset,ncol=2,scales="free_y")+
    geom_hline(yintercept=LOW_SD_PP,linetype="dashed",color="#C92A2A")+
    scale_fill_manual(values=PALETTE_CELL)+
    labs(title="BisqueRNA • desvio-padrão das proporções por coorte",
         subtitle="Linha tracejada = sinalização de variância baixa",
         x="Tipo celular",y="Desvio-padrão (pontos percentuais)")+base_theme()
  global_plot_status[[2L]] <- save_plot_07(
    pv,file.path(global_plot_dir,"07_variability_all.png"),12,8,180)
  global_plot_dt <- rbindlist(global_plot_status,use.names=TRUE,fill=TRUE)
  global_plot_dt[,scope:="global"]
  all_plot_dt <- rbindlist(
    c(lapply(runs, function(x) {z <- copy(x$plot_status);z[,scope:="dataset"];z}),
      list(global_plot_dt)),use.names=TRUE,fill=TRUE)
  write_tsv(all_plot_dt,file.path(RESULT_ROOT,"07_plot_status_all.tsv"))
  if (any(all_plot_dt$status=="WARN")) {
    message("ATENÇÃO: há figuras não geradas. Verifique 07_plot_status_all.tsv; ",
            "os resultados numéricos foram preservados.")
  }
  run_summary_view <- run_summary_all[,.(dataset,marker_set,bulk_samples,genes_used,
                                         degenerate_uniform,degenerate_constant,
                                         unique_vectors_6dec,median_uniform_rmse_pp)]
  say("\n=== RESUMO CONSOLIDADO BISQUE ===")
  print(run_summary_view)
  if (any(run_summary_all$degenerate_constant)) {
    say("ALERTA: resultados quase constantes. NÃO interpretar como composição confiável.")
  }
}

# Manifesto de TODOS os arquivos, inclusive plots e auditorias.
all_files <- list.files(RESULT_ROOT,recursive=TRUE,full.names=TRUE)
all_files <- all_files[file.info(all_files)$isdir==FALSE]
manifest <- data.table(path=all_files,exists=file.exists(all_files),
                       size_MB=as.numeric(file.info(all_files)$size)/1024^2)
write_tsv(manifest,file.path(RESULT_ROOT,"07_output_manifest.tsv"))

sink(file.path(RESULT_ROOT,"07_sessionInfo.txt"))
print(sessionInfo())
sink()

if (nrow(plan)!=length(runs) || any(status_all$status=="FAIL")) {
  fail("07 finalizado com falhas. Consulte ",file.path(RESULT_ROOT,"07_status.tsv"))
}
say("\n=== 07 CONCLUÍDO ===")
say("Outputs: %s",RESULT_ROOT)
say("Coortes concluídas nesta chamada: %s", paste(run_summary_all$dataset,collapse=", "))
say("Revisar: 07_run_summary.tsv, 07_celltype_variability_all.tsv e 07_status.tsv")
say("WARN de degenerescência deve ser investigado antes de usar Bisque como evidência.")
