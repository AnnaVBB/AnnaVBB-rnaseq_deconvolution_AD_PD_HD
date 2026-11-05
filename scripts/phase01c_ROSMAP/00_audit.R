# ==============================================================================
# 00_audit_ROSMAP.R
# ==============================================================================
#
# Fase 01c — ROSMAP
# Auditoria estrutural, técnica e de metadados ANTES da análise diferencial
#
# Objetivos
# ---------
# 1. Registrar os arquivos utilizados e sua proveniência (Synapse IDs).
# 2. Auditar a estrutura da matriz harmonizada de counts do ROSMAP batch 2.
# 3. Identificar as amostras DLPFC e confirmar a região no biospecimen metadata.
# 4. Validar o encadeamento:
#
#       counts/specimenID
#             ↓
#       assay metadata
#             ↓
#       biospecimen metadata
#             ↓
#       individualID
#             ↓
#       clinical metadata
#
# 5. Auditar flags oficiais de exclusão e seus motivos.
# 6. Quantificar os cenários fenotípicos PRÉ-DEFINIDOS:
#
#       Principal:
#         Control = cogdx 1
#         AD      = cogdx 4
#
#       Sensibilidade fenotípica futura:
#         Control = cogdx 1
#         AD      = cogdx 4 ou 5
#
# 7. Auditar covariáveis técnicas/demográficas relevantes antes do DESeq2.
# 8. Congelar e registrar, ANTES de ver os resultados de expressão diferencial,
#    os critérios estatísticos e os critérios de avaliação do ROSMAP.
#
# IMPORTANTE
# ----------
# Este script NÃO:
#   - executa DESeq2;
#   - filtra genes por expressão;
#   - converte Ensembl para GeneID;
#   - escolhe o ROSMAP como dataset definitivo;
#   - usa número de DEGs como critério de escolha.
#
# A finalidade é produzir uma auditoria reproduzível da entrada da Fase 01c.
#
# Critério principal já definido desde a Fase 01:
#
#       Benjamini-Hochberg FDR < 0.05
#       |log2FoldChange| > 0.58
#
# O valor 0.58 corresponde aproximadamente a fold-change de 1.5:
#
#       2^0.58 ≈ 1.50
#
# Critério estrito de efeito:
#
#       Benjamini-Hochberg FDR < 0.05
#       |log2FoldChange| > 1
#
# Critério nominal comparativo:
#
#       p < 0.05
#       |log2FoldChange| > 1
#
# O critério nominal é mantido apenas para comparabilidade histórica com as
# análises anteriores e NÃO substitui o critério principal baseado em FDR.
#
# ==============================================================================


# ==============================================================================
# 0. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Pacotes necessários não instalados: ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
})


# ==============================================================================
# 1. PARÂMETROS PRÉ-DEFINIDOS
# ==============================================================================

# --------------------------------------------------------------------------
# 1.1 Critérios de DEG
# --------------------------------------------------------------------------

PRIMARY_ALPHA <- 0.05
PRIMARY_LFC <- 0.58

STRICT_ALPHA <- 0.05
STRICT_LFC <- 1.00

NOMINAL_P <- 0.05
NOMINAL_LFC <- 1.00


# --------------------------------------------------------------------------
# 1.2 Definição fenotípica
# --------------------------------------------------------------------------
#
# Principal:
#   cogdx 1 = sem comprometimento cognitivo → Control
#   cogdx 4 = Alzheimer's dementia sem outra causa → AD
#
# Sensibilidade fenotípica futura:
#   Control = cogdx 1
#   AD      = cogdx 4 ou 5
#
# Nesta auditoria, esses cenários são apenas QUANTIFICADOS.
# O script 01_prepare_metadata_ROSMAP.R construirá o metadata analítico final.
# --------------------------------------------------------------------------

COGDX_CONTROL <- 1
COGDX_AD_STRICT <- 4
COGDX_AD_BROAD <- c(
  4,
  5
)


# --------------------------------------------------------------------------
# 1.3 Região principal
# --------------------------------------------------------------------------

TARGET_TISSUE <- "dorsolateral prefrontal cortex"
TARGET_SUFFIX <- "-DLPFC"


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

raw_dir <- file.path(
  "data",
  "raw",
  "ROSMAP"
)

legacy_test_dir <- file.path(
  "data",
  "raw",
  "ROSMAP_test"
)

if (
  !dir.exists(raw_dir) &&
  dir.exists(legacy_test_dir)
) {
  stop(
    paste0(
      "A pasta formal data/raw/ROSMAP ainda não existe, mas ",
      "data/raw/ROSMAP_test foi encontrada.\n\n",
      "Antes de executar a Fase 01c, renomeie:\n\n",
      "  mv data/raw/ROSMAP_test data/raw/ROSMAP\n"
    )
  )
}

if (!dir.exists(raw_dir)) {
  stop(
    "Diretório não encontrado: ",
    raw_dir
  )
}

counts_file <- file.path(
  raw_dir,
  "ROSMAP_batch2_gene_all_counts_matrix_clean.txt"
)

assay_file <- file.path(
  raw_dir,
  "ROSMAP_assay_rnaSeq_metadata.csv"
)

biospecimen_file <- file.path(
  raw_dir,
  "ROSMAP_biospecimen_metadata.csv"
)

clinical_file <- file.path(
  raw_dir,
  "ROSMAP_clinical.csv"
)

candidate_dir <- file.path(
  "results",
  "candidate_ROSMAP"
)

audit_dir <- file.path(
  candidate_dir,
  "audit"
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

audit_report_file <- file.path(
  audit_dir,
  "00_audit_ROSMAP.txt"
)

file_manifest_out <- file.path(
  audit_dir,
  "00_file_manifest.tsv"
)

synapse_manifest_out <- file.path(
  audit_dir,
  "00_synapse_manifest.tsv"
)

predefined_criteria_out <- file.path(
  audit_dir,
  "00_predefined_selection_criteria.tsv"
)

region_summary_out <- file.path(
  audit_dir,
  "00_counts_region_summary.tsv"
)

mapping_summary_out <- file.path(
  audit_dir,
  "00_metadata_mapping_summary.tsv"
)

cogdx_summary_out <- file.path(
  audit_dir,
  "00_cogdx_distribution.tsv"
)

exclusion_summary_out <- file.path(
  audit_dir,
  "00_exclusion_summary.tsv"
)

covariate_coverage_out <- file.path(
  audit_dir,
  "00_covariate_coverage.tsv"
)

strict_group_batch_out <- file.path(
  audit_dir,
  "00_strict_group_by_sequencingBatch.tsv"
)

strict_group_sex_out <- file.path(
  audit_dir,
  "00_strict_group_by_msex.tsv"
)

strict_braak_out <- file.path(
  audit_dir,
  "00_strict_group_by_braak.tsv"
)

strict_cerad_out <- file.path(
  audit_dir,
  "00_strict_group_by_cerad.tsv"
)

feature_audit_out <- file.path(
  audit_dir,
  "00_feature_id_audit.tsv"
)

technical_features_out <- file.path(
  audit_dir,
  "00_technical_features.tsv"
)

unmatched_out <- file.path(
  audit_dir,
  "00_unmatched_DLPFC_specimens.tsv"
)


# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================

assert_file <- function(path) {

  if (!file.exists(path)) {
    stop(
      "Arquivo não encontrado:\n",
      path
    )
  }

  invisible(TRUE)
}


assert_columns <- function(
  dt,
  required,
  object_name
) {

  missing <- setdiff(
    required,
    names(dt)
  )

  if (length(missing) > 0) {
    stop(
      object_name,
      ": colunas obrigatórias ausentes: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}


clean_character_columns <- function(dt) {

  char_cols <- names(dt)[
    vapply(
      dt,
      function(x) {
        is.character(x) ||
          inherits(
            x,
            "factor"
          )
      },
      logical(1)
    )
  ]

  for (col in char_cols) {

    set(
      dt,
      j = col,
      value = trimws(
        as.character(
          dt[[col]]
        )
      )
    )

    dt[
      get(col) == "",
      (col) := NA_character_
    ]
  }

  invisible(dt)
}


as_exclude_flag <- function(x) {

  if (is.logical(x)) {
    return(
      fifelse(
        is.na(x),
        FALSE,
        x
      )
    )
  }

  y <- tolower(
    trimws(
      as.character(x)
    )
  )

  y %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}


safe_numeric <- function(x) {

  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


fmt_int <- function(x) {

  format(
    as.integer(
      round(
        x
      )
    ),
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


fmt_num <- function(
  x,
  digits = 3
) {

  format(
    round(
      as.numeric(x),
      digits
    ),
    nsmall = digits,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


count_unique_nonmissing <- function(x) {

  uniqueN(
    x[
      !is.na(x)
    ]
  )
}


classify_header_region <- function(sample_ids) {

  fifelse(
    grepl(
      "-DLPFC$",
      sample_ids
    ),
    "DLPFC",
    fifelse(
      grepl(
        "-PCC$",
        sample_ids
      ),
      "PCC",
      fifelse(
        grepl(
          "-AC$",
          sample_ids
        ),
        "AC",
        "Other"
      )
    )
  )
}


count_fractional_values <- function(matrix_values) {

  sum(
    abs(
      matrix_values -
        round(
          matrix_values
        )
    ) > 1e-8,
    na.rm = TRUE
  )
}


write_crosstab <- function(
  dt,
  row_var,
  col_var,
  path
) {

  tmp <- dt[
    !is.na(
      get(row_var)
    ) &
      !is.na(
        get(col_var)
      )
  ]

  tab <- dcast(
    tmp,
    formula = as.formula(
      paste(
        row_var,
        "~",
        col_var
      )
    ),
    fun.aggregate = length,
    value.var = col_var
  )

  fwrite(
    tab,
    path,
    sep = "\t"
  )

  tab
}


# ==============================================================================
# 4. VALIDAR ARQUIVOS
# ==============================================================================

input_files <- c(
  counts_file,
  assay_file,
  biospecimen_file,
  clinical_file
)

invisible(
  lapply(
    input_files,
    assert_file
  )
)


# ==============================================================================
# 5. REGISTRAR PROVENIÊNCIA DOS ARQUIVOS
# ==============================================================================

synapse_manifest <- data.table(
  role = c(
    "ROSMAP harmonized RNA-seq counts — batch 2",
    "ROSMAP clinical metadata",
    "ROSMAP biospecimen metadata",
    "ROSMAP RNA-seq assay metadata"
  ),
  synapse_id = c(
    "syn22301601",
    "syn3191087",
    "syn21323366",
    "syn21088596"
  ),
  version_number = c(
    NA_integer_,
    NA_integer_,
    NA_integer_,
    NA_integer_
  ),
  access = rep(
    "controlled access; READ + DOWNLOAD permission confirmed",
    4
  ),
  local_file = c(
    basename(counts_file),
    basename(clinical_file),
    basename(biospecimen_file),
    basename(assay_file)
  )
)

fwrite(
  synapse_manifest,
  synapse_manifest_out,
  sep = "\t",
  na = "NA"
)

file_info <- file.info(
  input_files
)

file_manifest <- data.table(
  role = c(
    "counts",
    "assay_metadata",
    "biospecimen_metadata",
    "clinical_metadata"
  ),
  path = normalizePath(
    input_files,
    winslash = "/",
    mustWork = TRUE
  ),
  size_bytes = file_info$size,
  modified_time = as.character(
    file_info$mtime
  ),
  md5 = unname(
    tools::md5sum(
      input_files
    )
  )
)

fwrite(
  file_manifest,
  file_manifest_out,
  sep = "\t"
)


# ==============================================================================
# 6. CONGELAR CRITÉRIOS ANTES DO DESeq2
# ==============================================================================

predefined_criteria <- data.table(
  category = c(
    "DEG principal",
    "DEG estrito / efeito forte",
    "Comparação nominal histórica",
    "Fenótipo principal ROSMAP",
    "Sensibilidade fenotípica futura",
    "Seleção do candidato — qualidade técnica",
    "Seleção do candidato — estabilidade",
    "Seleção do candidato — replicabilidade externa",
    "Seleção do candidato — coerência funcional",
    "Seleção do candidato — compatibilidade anatômica",
    "Critério explicitamente NÃO utilizado"
  ),
  definition = c(
    paste0(
      "Benjamini-Hochberg FDR < ",
      PRIMARY_ALPHA,
      " e |log2FC| > ",
      PRIMARY_LFC
    ),
    paste0(
      "Benjamini-Hochberg FDR < ",
      STRICT_ALPHA,
      " e |log2FC| > ",
      STRICT_LFC
    ),
    paste0(
      "p < ",
      NOMINAL_P,
      " e |log2FC| > ",
      NOMINAL_LFC,
      "; apenas para comparabilidade histórica"
    ),
    "DLPFC; exclude != TRUE; Control = cogdx 1; AD = cogdx 4",
    "DLPFC; exclude != TRUE; Control = cogdx 1; AD = cogdx 4 ou 5",
    paste0(
      "Mapeamento counts↔metadata íntegro, QC adequado e ausência de ",
      "confundimento completo diagnóstico↔variáveis técnicas"
    ),
    paste0(
      "Concordância global dos log2FC e manutenção predominante de direção ",
      "entre modelo principal e modelo ajustado por covariáveis"
    ),
    paste0(
      "Concordância direcional/log2FC com GSE125583 no universo comparável; ",
      "não exigir replicação gene-a-gene completa por diferença de região"
    ),
    "Convergência de processos/vias entre ROSMAP e GSE125583 por ORA/GSEA",
    paste0(
      "DLPFC do ROSMAP favorece alinhamento com referência Mathys 2019 ",
      "e com a região pré-frontal dos datasets de HD/PD"
    ),
    paste0(
      "Número absoluto de DEGs não será usado isoladamente para escolher ",
      "o dataset principal de AD"
    )
  )
)

fwrite(
  predefined_criteria,
  predefined_criteria_out,
  sep = "\t"
)


# ==============================================================================
# 7. AUDITAR CABEÇALHO DA MATRIZ DE COUNTS
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "00 — AUDITORIA ROSMAP\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Lendo cabeçalho da matriz de counts...\n"
)

counts_header <- fread(
  counts_file,
  nrows = 0
)

header_names <- names(
  counts_header
)

if (length(header_names) < 2) {
  stop(
    "A matriz de counts possui menos de duas colunas."
  )
}

feature_col <- header_names[1]

if (feature_col != "feature") {
  warning(
    "A primeira coluna da matriz se chama '",
    feature_col,
    "', e não 'feature'."
  )
}

count_sample_ids <- header_names[
  -1
]

header_region <- classify_header_region(
  count_sample_ids
)

region_summary <- data.table(
  region_by_header = c(
    "DLPFC",
    "PCC",
    "AC",
    "Other"
  ),
  n_samples = c(
    sum(
      header_region == "DLPFC"
    ),
    sum(
      header_region == "PCC"
    ),
    sum(
      header_region == "AC"
    ),
    sum(
      header_region == "Other"
    )
  )
)

fwrite(
  region_summary,
  region_summary_out,
  sep = "\t"
)

dlpfc_ids <- count_sample_ids[
  header_region == "DLPFC"
]

cat(
  "Colunas de amostra na matriz:",
  fmt_int(
    length(
      count_sample_ids
    )
  ),
  "\n"
)

cat(
  "DLPFC pelo cabeçalho:",
  fmt_int(
    length(
      dlpfc_ids
    )
  ),
  "\n\n"
)


# ==============================================================================
# 8. AUDITAR FEATURES SEM CARREGAR TODAS AS COLUNAS
# ==============================================================================

cat(
  "Auditando coluna de features...\n"
)

feature_dt <- fread(
  counts_file,
  select = 1,
  na.strings = c(
    "NA",
    ""
  )
)

setnames(
  feature_dt,
  1,
  "feature"
)

feature_dt[
  ,
  feature := as.character(feature)
]

n_rows_total <- nrow(
  feature_dt
)

n_feature_missing <- sum(
  is.na(
    feature_dt$feature
  ) |
    trimws(
      feature_dt$feature
    ) == ""
)

n_feature_duplicates <- anyDuplicated(
  feature_dt$feature
)

technical_mask <- grepl(
  "^N_",
  feature_dt$feature
)

technical_features <- feature_dt[
  technical_mask
]

gene_features <- feature_dt[
  !technical_mask
]

n_gene_features <- nrow(
  gene_features
)

n_gene_feature_duplicates <- anyDuplicated(
  gene_features$feature
)

ensembl_versioned_pattern <- "^ENSG[0-9]+\\.[0-9]+$"
ensembl_unversioned_pattern <- "^ENSG[0-9]+$"

n_ensembl_versioned <- sum(
  grepl(
    ensembl_versioned_pattern,
    gene_features$feature
  )
)

n_ensembl_unversioned <- sum(
  grepl(
    ensembl_unversioned_pattern,
    gene_features$feature
  )
)

n_other_gene_features <- (
  n_gene_features -
    n_ensembl_versioned -
    n_ensembl_unversioned
)

feature_audit <- data.table(
  metric = c(
    "rows_after_header",
    "technical_N_rows",
    "gene_feature_rows",
    "missing_feature_ids",
    "duplicated_feature_ids_any",
    "duplicated_gene_feature_ids_any",
    "versioned_Ensembl_gene_ids",
    "unversioned_Ensembl_gene_ids",
    "other_gene_feature_ids"
  ),
  value = c(
    n_rows_total,
    nrow(
      technical_features
    ),
    n_gene_features,
    n_feature_missing,
    as.integer(
      n_feature_duplicates > 0
    ),
    as.integer(
      n_gene_feature_duplicates > 0
    ),
    n_ensembl_versioned,
    n_ensembl_unversioned,
    n_other_gene_features
  )
)

fwrite(
  feature_audit,
  feature_audit_out,
  sep = "\t"
)

fwrite(
  technical_features,
  technical_features_out,
  sep = "\t"
)


# ==============================================================================
# 9. LER APENAS AS COLUNAS DLPFC E AUDITAR COUNTS
# ==============================================================================

cat(
  "Lendo apenas as ",
  length(
    dlpfc_ids
  ),
  " colunas DLPFC para auditoria de integridade...\n",
  sep = ""
)

dlpfc_counts <- fread(
  counts_file,
  select = c(
    feature_col,
    dlpfc_ids
  ),
  na.strings = c(
    "NA",
    ""
  )
)

setnames(
  dlpfc_counts,
  1,
  "feature"
)

if (
  nrow(
    dlpfc_counts
  ) != n_rows_total
) {
  stop(
    "Número de linhas mudou ao ler as colunas DLPFC."
  )
}

dlpfc_technical <- grepl(
  "^N_",
  as.character(
    dlpfc_counts$feature
  )
)

dlpfc_gene_counts <- dlpfc_counts[
  !dlpfc_technical
]

count_matrix <- as.matrix(
  dlpfc_gene_counts[
    ,
    -1,
    with = FALSE
  ]
)

storage.mode(
  count_matrix
) <- "numeric"

n_count_na <- sum(
  is.na(
    count_matrix
  )
)

n_count_nonfinite <- sum(
  !is.finite(
    count_matrix
  ),
  na.rm = TRUE
)

n_count_negative <- sum(
  count_matrix < 0,
  na.rm = TRUE
)

n_count_fractional <- count_fractional_values(
  count_matrix
)

n_all_zero_genes <- sum(
  rowSums(
    count_matrix,
    na.rm = TRUE
  ) == 0
)

n_all_zero_samples <- sum(
  colSums(
    count_matrix,
    na.rm = TRUE
  ) == 0
)

library_sizes <- colSums(
  count_matrix,
  na.rm = TRUE
)

detected_genes <- colSums(
  count_matrix > 0,
  na.rm = TRUE
)


# ==============================================================================
# 10. CARREGAR E VALIDAR METADADOS
# ==============================================================================

cat(
  "Lendo metadados...\n"
)

assay <- fread(
  assay_file,
  na.strings = c(
    "NA",
    ""
  )
)

biospecimen <- fread(
  biospecimen_file,
  na.strings = c(
    "NA",
    ""
  )
)

clinical <- fread(
  clinical_file,
  na.strings = c(
    "NA",
    ""
  )
)

clean_character_columns(
  assay
)

clean_character_columns(
  biospecimen
)

clean_character_columns(
  clinical
)

assert_columns(
  assay,
  c(
    "specimenID",
    "RIN",
    "sequencingBatch"
  ),
  "assay metadata"
)

assert_columns(
  biospecimen,
  c(
    "individualID",
    "specimenID",
    "organ",
    "tissue",
    "BrodmannArea",
    "exclude",
    "excludeReason"
  ),
  "biospecimen metadata"
)

assert_columns(
  clinical,
  c(
    "individualID",
    "msex",
    "age_death",
    "pmi",
    "braaksc",
    "ceradsc",
    "cogdx"
  ),
  "clinical metadata"
)


# ==============================================================================
# 11. AUDITAR UNICIDADE DAS CHAVES
# ==============================================================================
#
# IMPORTANTE:
# O arquivo global de biospecimen do ROSMAP contém múltiplos registros por
# specimenID em outras regiões/assays. Isso é esperado no arquivo completo e
# NÃO constitui, por si só, erro para a análise DLPFC.
#
# A unicidade relevante deve ser verificada APENAS entre os 210 specimenIDs
# DLPFC presentes na matriz de counts.
# ==============================================================================

assay_dup_global <- assay[
  !is.na(specimenID),
  .N,
  by = specimenID
][
  N > 1
]

biospecimen_dup_global <- biospecimen[
  !is.na(specimenID),
  .N,
  by = specimenID
][
  N > 1
]

clinical_dup_global <- clinical[
  !is.na(individualID),
  .N,
  by = individualID
][
  N > 1
]

# O assay metadata precisa ser unívoco para specimenID.
if (nrow(assay_dup_global) > 0) {
  stop(
    "assay metadata: existem specimenID duplicados no arquivo completo."
  )
}

# O clinical metadata precisa ser unívoco por individualID.
if (nrow(clinical_dup_global) > 0) {
  stop(
    "clinical metadata: existem individualID duplicados."
  )
}

# No biospecimen, duplicatas globais são registradas, mas só se tornam
# problemáticas se atingirem as amostras DLPFC efetivamente analisadas.
bio_dlpfc <- biospecimen[
  specimenID %in% dlpfc_ids
]

assay_dlpfc <- assay[
  specimenID %in% dlpfc_ids
]

biospecimen_dup_dlpfc <- bio_dlpfc[
  !is.na(specimenID),
  .N,
  by = specimenID
][
  N > 1
]

assay_dup_dlpfc <- assay_dlpfc[
  !is.na(specimenID),
  .N,
  by = specimenID
][
  N > 1
]

if (nrow(biospecimen_dup_dlpfc) > 0) {
  stop(
    "biospecimen metadata: existem specimenID duplicados entre as amostras DLPFC."
  )
}

if (nrow(assay_dup_dlpfc) > 0) {
  stop(
    "assay metadata: existem specimenID duplicados entre as amostras DLPFC."
  )
}

# Salvar resumo das duplicatas globais e DLPFC para documentação.
duplicate_key_summary_out <- file.path(
  audit_dir,
  "00_duplicate_key_summary.tsv"
)

duplicate_key_summary <- data.table(
  source = c(
    "biospecimen_global",
    "biospecimen_DLPFC",
    "assay_global",
    "assay_DLPFC",
    "clinical_global"
  ),
  n_duplicated_ids = c(
    nrow(biospecimen_dup_global),
    nrow(biospecimen_dup_dlpfc),
    nrow(assay_dup_global),
    nrow(assay_dup_dlpfc),
    nrow(clinical_dup_global)
  )
)

fwrite(
  duplicate_key_summary,
  duplicate_key_summary_out,
  sep = "\t"
)


# ==============================================================================
# 12. MAPEAMENTO COUNTS → ASSAY / BIOSPECIMEN
# ==============================================================================

sample_dt <- data.table(
  specimenID = as.character(
    dlpfc_ids
  )
)

assay_ids <- as.character(
  assay$specimenID
)

biospecimen_ids <- as.character(
  biospecimen$specimenID
)

matched_assay <- sample_dt$specimenID %in% assay_ids
matched_biospecimen <- sample_dt$specimenID %in% biospecimen_ids

mapping_summary <- data.table(
  metric = c(
    "DLPFC columns in counts",
    "matched to assay metadata",
    "unmatched to assay metadata",
    "matched to biospecimen metadata",
    "unmatched to biospecimen metadata"
  ),
  n = c(
    nrow(
      sample_dt
    ),
    sum(
      matched_assay
    ),
    sum(
      !matched_assay
    ),
    sum(
      matched_biospecimen
    ),
    sum(
      !matched_biospecimen
    )
  )
)

fwrite(
  mapping_summary,
  mapping_summary_out,
  sep = "\t"
)

unmatched_dt <- data.table(
  specimenID = sample_dt$specimenID,
  matched_assay = matched_assay,
  matched_biospecimen = matched_biospecimen
)[
  !matched_assay |
    !matched_biospecimen
]

fwrite(
  unmatched_dt,
  unmatched_out,
  sep = "\t"
)


# ==============================================================================
# 13. CONSTRUIR TABELA AUDITADA DLPFC
# ==============================================================================

audit_meta <- merge(
  sample_dt,
  assay,
  by = "specimenID",
  all.x = TRUE,
  sort = FALSE
)

audit_meta <- merge(
  audit_meta,
  bio_dlpfc,
  by = "specimenID",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c(
    "_assay",
    "_biospecimen"
  )
)

n_individual_missing <- sum(
  is.na(
    audit_meta$individualID
  )
)

n_unique_individuals <- uniqueN(
  audit_meta$individualID[
    !is.na(
      audit_meta$individualID
    )
  ]
)

individual_counts <- audit_meta[
  !is.na(
    individualID
  ),
  .N,
  by = individualID
]

n_individuals_multiple_dlpfc <- sum(
  individual_counts$N > 1
)

audit_meta <- merge(
  audit_meta,
  clinical,
  by = "individualID",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c(
    "",
    "_clinical"
  )
)


# ==============================================================================
# 14. CONFIRMAR REGIÃO NO BIOSPECIMEN METADATA
# ==============================================================================

tissue_normalized <- tolower(
  trimws(
    as.character(
      audit_meta$tissue
    )
  )
)

n_target_tissue <- sum(
  tissue_normalized == tolower(
    TARGET_TISSUE
  ),
  na.rm = TRUE
)

n_non_target_tissue <- sum(
  !is.na(
    tissue_normalized
  ) &
    tissue_normalized != tolower(
      TARGET_TISSUE
    )
)

n_tissue_missing <- sum(
  is.na(
    audit_meta$tissue
  )
)

organ_distribution <- audit_meta[
  ,
  .N,
  by = organ
][
  order(
    -N
  )
]

tissue_distribution <- audit_meta[
  ,
  .N,
  by = tissue
][
  order(
    -N
  )
]

brodmann_distribution <- audit_meta[
  ,
  .N,
  by = BrodmannArea
][
  order(
    -N
  )
]


# ==============================================================================
# 15. FLAGS DE EXCLUSÃO
# ==============================================================================

audit_meta[
  ,
  exclude_flag := as_exclude_flag(
    exclude
  )
]

exclusion_summary <- audit_meta[
  ,
  .N,
  by = .(
    exclude = as.character(
      exclude
    ),
    excludeReason
  )
][
  order(
    -N
  )
]

fwrite(
  exclusion_summary,
  exclusion_summary_out,
  sep = "\t",
  na = "NA"
)

n_excluded <- sum(
  audit_meta$exclude_flag
)

n_not_excluded <- sum(
  !audit_meta$exclude_flag
)


# ==============================================================================
# 16. COGDX E CENÁRIOS FENOTÍPICOS
# ==============================================================================

audit_meta[
  ,
  cogdx_numeric := safe_numeric(
    cogdx
  )
]

cogdx_distribution <- audit_meta[
  ,
  .N,
  by = .(
    exclude_flag,
    cogdx = cogdx_numeric
  )
][
  order(
    exclude_flag,
    cogdx
  )
]

fwrite(
  cogdx_distribution,
  cogdx_summary_out,
  sep = "\t",
  na = "NA"
)

valid_dlpfc <- audit_meta[
  exclude_flag == FALSE
]

strict_cohort <- valid_dlpfc[
  cogdx_numeric %in% c(
    COGDX_CONTROL,
    COGDX_AD_STRICT
  )
]

strict_cohort[
  ,
  group := fifelse(
    cogdx_numeric == COGDX_CONTROL,
    "Control",
    "AD"
  )
]

strict_group_counts <- strict_cohort[
  ,
  .N,
  by = group
][
  order(
    group
  )
]

broad_cohort <- valid_dlpfc[
  cogdx_numeric %in% c(
    COGDX_CONTROL,
    COGDX_AD_BROAD
  )
]

broad_cohort[
  ,
  group := fifelse(
    cogdx_numeric == COGDX_CONTROL,
    "Control",
    "AD"
  )
]

broad_group_counts <- broad_cohort[
  ,
  .N,
  by = group
][
  order(
    group
  )
]


# ==============================================================================
# 17. COBERTURA DE COVARIÁVEIS
# ==============================================================================

covariates <- c(
  "msex",
  "age_death",
  "pmi",
  "braaksc",
  "ceradsc",
  "RIN",
  "rnaBatch",
  "libraryBatch",
  "sequencingBatch"
)

covariate_coverage <- rbindlist(
  lapply(
    covariates,
    function(variable) {

      if (!variable %in% names(audit_meta)) {

        return(
          data.table(
            variable = variable,
            n_total_DLPFC = nrow(
              audit_meta
            ),
            n_nonmissing = 0L,
            n_missing = nrow(
              audit_meta
            ),
            n_unique_nonmissing = 0L,
            present_in_metadata = FALSE
          )
        )
      }

      x <- audit_meta[[variable]]

      data.table(
        variable = variable,
        n_total_DLPFC = nrow(
          audit_meta
        ),
        n_nonmissing = sum(
          !is.na(x)
        ),
        n_missing = sum(
          is.na(x)
        ),
        n_unique_nonmissing = count_unique_nonmissing(
          x
        ),
        present_in_metadata = TRUE
      )
    }
  )
)

fwrite(
  covariate_coverage,
  covariate_coverage_out,
  sep = "\t"
)


# ==============================================================================
# 18. AUDITORIA DE age_death
# ==============================================================================

age_death_text <- as.character(
  strict_cohort$age_death
)

n_age_90plus_strict <- sum(
  trimws(
    age_death_text
  ) == "90+",
  na.rm = TRUE
)

age_exact_numeric <- suppressWarnings(
  as.numeric(
    age_death_text[
      trimws(
        age_death_text
      ) != "90+"
    ]
  )
)

age_exact_numeric <- age_exact_numeric[
  !is.na(
    age_exact_numeric
  )
]


# ==============================================================================
# 19. CROSSTABS NO CENÁRIO PRINCIPAL
# ==============================================================================

strict_batch_table <- write_crosstab(
  strict_cohort,
  "sequencingBatch",
  "group",
  strict_group_batch_out
)

strict_sex_table <- write_crosstab(
  strict_cohort,
  "msex",
  "group",
  strict_group_sex_out
)

strict_braak_table <- write_crosstab(
  strict_cohort,
  "braaksc",
  "group",
  strict_braak_out
)

strict_cerad_table <- write_crosstab(
  strict_cohort,
  "ceradsc",
  "group",
  strict_cerad_out
)


# ==============================================================================
# 20. VERIFICAÇÕES ESTRUTURAIS
# ==============================================================================

if (
  any(
    !matched_assay
  )
) {
  warning(
    "Existem amostras DLPFC sem correspondência no assay metadata."
  )
}

if (
  any(
    !matched_biospecimen
  )
) {
  warning(
    "Existem amostras DLPFC sem correspondência no biospecimen metadata."
  )
}

if (
  n_non_target_tissue > 0 ||
    n_tissue_missing > 0
) {
  warning(
    paste0(
      "Nem todas as colunas selecionadas como DLPFC pelo cabeçalho foram ",
      "confirmadas como '",
      TARGET_TISSUE,
      "' no biospecimen metadata."
    )
  )
}

if (n_count_na > 0) {
  stop(
    "Foram encontrados valores NA na matriz DLPFC."
  )
}

if (n_count_nonfinite > 0) {
  stop(
    "Foram encontrados valores não finitos na matriz DLPFC."
  )
}

if (n_count_negative > 0) {
  stop(
    "Foram encontrados counts negativos na matriz DLPFC."
  )
}

if (n_count_fractional > 0) {
  stop(
    "Foram encontrados counts fracionários na matriz DLPFC."
  )
}

if (n_all_zero_samples > 0) {
  stop(
    "Foram encontradas amostras DLPFC totalmente zeradas."
  )
}

if (n_individual_missing > 0) {
  warning(
    "Existem amostras DLPFC sem individualID."
  )
}

if (n_individuals_multiple_dlpfc > 0) {
  warning(
    "Existem indivíduos com mais de uma amostra DLPFC."
  )
}

batch_by_group <- table(
  strict_cohort$sequencingBatch,
  strict_cohort$group,
  useNA = "no"
)

perfect_batch_separation <- FALSE

if (
  nrow(
    batch_by_group
  ) > 0 &&
    ncol(
      batch_by_group
    ) == 2
) {

  rows_with_both_groups <- rowSums(
    batch_by_group > 0
  ) == 2

  if (
    !any(
      rows_with_both_groups
    )
  ) {
    perfect_batch_separation <- TRUE
  }
}


# ==============================================================================
# 21. ESCREVER RELATÓRIO TXT
# ==============================================================================

sink(
  audit_report_file
)

cat(
  "============================================================\n"
)
cat(
  "FASE 01c — 00_audit_ROSMAP.R\n"
)
cat(
  "AUDITORIA ESTRUTURAL PRÉ-DESeq2\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "1. PROVENIÊNCIA / SYNAPSE\n"
)
cat(
  "------------------------------------------------------------\n"
)

for (
  i in seq_len(
    nrow(
      synapse_manifest
    )
  )
) {

  cat(
    synapse_manifest$role[i],
    "\n"
  )

  cat(
    "  Synapse ID: ",
    synapse_manifest$synapse_id[i],
    "\n",
    sep = ""
  )

  cat(
    "  Version number: ",
    ifelse(
      is.na(
        synapse_manifest$version_number[i]
      ),
      "NA — preencher após consultar a entidade no Synapse",
      synapse_manifest$version_number[i]
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Arquivo local: ",
    synapse_manifest$local_file[i],
    "\n",
    sep = ""
  )

  cat(
    "  Acesso: ",
    synapse_manifest$access[i],
    "\n\n",
    sep = ""
  )
}

cat(
  "------------------------------------------------------------\n"
)
cat(
  "2. CRITÉRIOS DEFINIDOS ANTES DO DESeq2 DO ROSMAP\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Critério principal de DEG:\n"
)
cat(
  "  Benjamini-Hochberg FDR < ",
  PRIMARY_ALPHA,
  "\n",
  sep = ""
)
cat(
  "  |log2FC| > ",
  PRIMARY_LFC,
  "\n\n",
  sep = ""
)

cat(
  "Critério estrito / efeito forte:\n"
)
cat(
  "  Benjamini-Hochberg FDR < ",
  STRICT_ALPHA,
  "\n",
  sep = ""
)
cat(
  "  |log2FC| > ",
  STRICT_LFC,
  "\n\n",
  sep = ""
)

cat(
  "Critério nominal comparativo:\n"
)
cat(
  "  p < ",
  NOMINAL_P,
  "\n",
  sep = ""
)
cat(
  "  |log2FC| > ",
  NOMINAL_LFC,
  "\n",
  sep = ""
)
cat(
  "  Uso: apenas comparabilidade histórica.\n\n"
)

cat(
  "Fenótipo principal ROSMAP:\n"
)
cat(
  "  DLPFC\n"
)
cat(
  "  excluir somente biospecimen com exclude=TRUE\n"
)
cat(
  "  Control = cogdx 1\n"
)
cat(
  "  AD      = cogdx 4\n\n"
)

cat(
  "Sensibilidade fenotípica futura:\n"
)
cat(
  "  Control = cogdx 1\n"
)
cat(
  "  AD      = cogdx 4 ou 5\n\n"
)

cat(
  "O número absoluto de DEGs NÃO será usado isoladamente para escolher\n"
)
cat(
  "o dataset principal de Alzheimer.\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "3. ESTRUTURA DA MATRIZ DE COUNTS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Arquivo: ",
  basename(
    counts_file
  ),
  "\n",
  sep = ""
)

cat(
  "Linhas após cabeçalho: ",
  fmt_int(
    n_rows_total
  ),
  "\n",
  sep = ""
)

cat(
  "Colunas de amostra: ",
  fmt_int(
    length(
      count_sample_ids
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Distribuição pelo sufixo das colunas:\n"
)

print(
  region_summary
)

cat(
  "\nFeatures técnicas N_: ",
  fmt_int(
    nrow(
      technical_features
    )
  ),
  "\n",
  sep = ""
)

if (
  nrow(
    technical_features
  ) > 0
) {

  cat(
    "Features técnicas:\n"
  )

  cat(
    paste0(
      "  ",
      technical_features$feature,
      collapse = "\n"
    ),
    "\n"
  )
}

cat(
  "\nFeatures gênicas após remover N_: ",
  fmt_int(
    n_gene_features
  ),
  "\n",
  sep = ""
)

cat(
  "Ensembl versionados: ",
  fmt_int(
    n_ensembl_versioned
  ),
  "\n",
  sep = ""
)

cat(
  "Ensembl sem versão: ",
  fmt_int(
    n_ensembl_unversioned
  ),
  "\n",
  sep = ""
)

cat(
  "Outros IDs entre features gênicas: ",
  fmt_int(
    n_other_gene_features
  ),
  "\n\n",
  sep = ""
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "4. INTEGRIDADE DOS COUNTS — DLPFC\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Amostras DLPFC: ",
  fmt_int(
    ncol(
      count_matrix
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Features gênicas: ",
  fmt_int(
    nrow(
      count_matrix
    )
  ),
  "\n",
  sep = ""
)

cat(
  "NA: ",
  fmt_int(
    n_count_na
  ),
  "\n",
  sep = ""
)

cat(
  "Não finitos: ",
  fmt_int(
    n_count_nonfinite
  ),
  "\n",
  sep = ""
)

cat(
  "Negativos: ",
  fmt_int(
    n_count_negative
  ),
  "\n",
  sep = ""
)

cat(
  "Fracionários: ",
  fmt_int(
    n_count_fractional
  ),
  "\n",
  sep = ""
)

cat(
  "Genes totalmente zerados: ",
  fmt_int(
    n_all_zero_genes
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras totalmente zeradas: ",
  fmt_int(
    n_all_zero_samples
  ),
  "\n\n",
  sep = ""
)

cat(
  "Library size:\n"
)

cat(
  "  min = ",
  fmt_int(
    min(
      library_sizes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  mediana = ",
  fmt_int(
    median(
      library_sizes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  max = ",
  fmt_int(
    max(
      library_sizes
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Genes detectados por amostra (>0):\n"
)

cat(
  "  min = ",
  fmt_int(
    min(
      detected_genes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  mediana = ",
  fmt_int(
    median(
      detected_genes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "  max = ",
  fmt_int(
    max(
      detected_genes
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "5. AUDITORIA DE UNICIDADE DAS CHAVES\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Biospecimen — specimenIDs duplicados no arquivo completo: ",
  fmt_int(
    nrow(
      biospecimen_dup_global
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Biospecimen — specimenIDs duplicados entre DLPFC: ",
  fmt_int(
    nrow(
      biospecimen_dup_dlpfc
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Assay — specimenIDs duplicados no arquivo completo: ",
  fmt_int(
    nrow(
      assay_dup_global
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Assay — specimenIDs duplicados entre DLPFC: ",
  fmt_int(
    nrow(
      assay_dup_dlpfc
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Clinical — individualIDs duplicados no arquivo completo: ",
  fmt_int(
    nrow(
      clinical_dup_global
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Interpretação: duplicatas globais no biospecimen não impedem a análise,\n"
)
cat(
  "pois nenhuma delas atinge os 210 specimenIDs DLPFC presentes nos counts.\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "6. MAPEAMENTO COUNTS ↔ METADATA\n"
)
cat(
  "------------------------------------------------------------\n"
)

print(
  mapping_summary
)

cat(
  "\nDLPFC com individualID ausente: ",
  fmt_int(
    n_individual_missing
  ),
  "\n",
  sep = ""
)

cat(
  "Indivíduos DLPFC únicos: ",
  fmt_int(
    n_unique_individuals
  ),
  "\n",
  sep = ""
)

cat(
  "Indivíduos com >1 amostra DLPFC: ",
  fmt_int(
    n_individuals_multiple_dlpfc
  ),
  "\n\n",
  sep = ""
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "6. REGIÃO ANATÔMICA CONFIRMADA NO BIOSPECIMEN METADATA\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Target tissue: ",
  TARGET_TISSUE,
  "\n\n",
  sep = ""
)

cat(
  "DLPFC confirmados pelo campo tissue: ",
  fmt_int(
    n_target_tissue
  ),
  " / ",
  fmt_int(
    nrow(
      audit_meta
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "organ:\n"
)

print(
  organ_distribution
)

cat(
  "\ntissue:\n"
)

print(
  tissue_distribution
)

cat(
  "\nBrodmannArea:\n"
)

print(
  brodmann_distribution
)

cat(
  "\nOBSERVAÇÃO: se BrodmannArea estiver ausente, o relatório final deve\n"
)
cat(
  "descrever estas amostras como DLPFC, sem afirmar BA9 com base neste\n"
)
cat(
  "metadata específico.\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "7. FLAGS OFICIAIS DE EXCLUSÃO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "DLPFC total: ",
  fmt_int(
    nrow(
      audit_meta
    )
  ),
  "\n",
  sep = ""
)

cat(
  "exclude=TRUE: ",
  fmt_int(
    n_excluded
  ),
  "\n",
  sep = ""
)

cat(
  "DLPFC após exclusões oficiais: ",
  fmt_int(
    n_not_excluded
  ),
  "\n\n",
  sep = ""
)

print(
  exclusion_summary
)

cat(
  "\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "8. DISTRIBUIÇÃO cogdx\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Antes/depois das flags oficiais:\n"
)

print(
  cogdx_distribution
)

cat(
  "\nCenário principal após exclusões:\n"
)

print(
  strict_group_counts
)

cat(
  "\nTotal principal AD + Control: ",
  fmt_int(
    nrow(
      strict_cohort
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Cenário amplo após exclusões (sensibilidade fenotípica futura):\n"
)

print(
  broad_group_counts
)

cat(
  "\nTotal amplo AD + Control: ",
  fmt_int(
    nrow(
      broad_cohort
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "9. COBERTURA DE COVARIÁVEIS — 210 DLPFC\n"
)
cat(
  "------------------------------------------------------------\n"
)

print(
  covariate_coverage
)

cat(
  "\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "10. age_death NO CENÁRIO PRINCIPAL\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Amostras no cenário principal: ",
  fmt_int(
    nrow(
      strict_cohort
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Valores codificados como 90+: ",
  fmt_int(
    n_age_90plus_strict
  ),
  "\n",
  sep = ""
)

if (
  length(
    age_exact_numeric
  ) > 0
) {

  cat(
    "Entre idades exatas (<90+), min/mediana/max: ",
    fmt_num(
      min(
        age_exact_numeric
      ),
      2
    ),
    " / ",
    fmt_num(
      median(
        age_exact_numeric
      ),
      2
    ),
    " / ",
    fmt_num(
      max(
        age_exact_numeric
      ),
      2
    ),
    "\n",
    sep = ""
  )
}

cat(
  "\nA codificação 90+ NÃO foi convertida para idade exata nesta auditoria.\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "11. CENÁRIO PRINCIPAL — sequencingBatch × group\n"
)
cat(
  "------------------------------------------------------------\n"
)

print(
  strict_batch_table
)

cat(
  "\nConfundimento completo group↔sequencingBatch: ",
  perfect_batch_separation,
  "\n\n",
  sep = ""
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "12. CENÁRIO PRINCIPAL — msex × group\n"
)
cat(
  "------------------------------------------------------------\n"
)

print(
  strict_sex_table
)

cat(
  "\nOBSERVAÇÃO: a codificação 0/1 de msex não é renomeada neste script;\n"
)
cat(
  "o significado dos níveis deve ser confirmado no codebook antes do modelo.\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "13. CENÁRIO PRINCIPAL — Braak × group\n"
)
cat(
  "------------------------------------------------------------\n"
)

print(
  strict_braak_table
)

cat(
  "\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "14. CENÁRIO PRINCIPAL — CERAD × group\n"
)
cat(
  "------------------------------------------------------------\n"
)

print(
  strict_cerad_table
)

cat(
  "\nBraak e CERAD são auditados, mas NÃO são definidos aqui como\n"
)
cat(
  "covariáveis do modelo principal de AD, pois são características\n"
)
cat(
  "neuropatológicas diretamente relacionadas ao fenótipo de doença.\n\n"
)

cat(
  "------------------------------------------------------------\n"
)
cat(
  "15. CONCLUSÃO DA AUDITORIA\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "- A auditoria foi executada ANTES do DESeq2 do ROSMAP.\n"
)

cat(
  "- O critério principal permanece BH-FDR < 0,05 + |log2FC| > 0,58.\n"
)

cat(
  "- A região-alvo formal é DLPFC.\n"
)

cat(
  "- O biospecimen metadata possui duplicatas globais, mas nenhuma entre\n"
)
cat(
  "  os 210 specimenIDs DLPFC presentes na matriz de counts.\n"
)

cat(
  "- O cenário fenotípico principal é cogdx 1 vs 4, após exclusões oficiais.\n"
)

cat(
  "- O cenário cogdx 1 vs 4+5 fica reservado para sensibilidade fenotípica.\n"
)

cat(
  "- O número absoluto de DEGs não será usado isoladamente para escolher\n"
)
cat(
  "  o dataset principal de Alzheimer.\n"
)

cat(
  "- A escolha posterior deverá considerar conjuntamente qualidade técnica,\n"
)
cat(
  "  estabilidade com covariáveis, concordância de efeitos com GSE125583,\n"
)
cat(
  "  coerência funcional e compatibilidade anatômica.\n"
)

if (
  perfect_batch_separation
) {

  cat(
    "- ATENÇÃO: foi detectada separação completa entre group e sequencingBatch.\n"
  )

} else {

  cat(
    "- Não foi detectado confundimento COMPLETO entre group e sequencingBatch.\n"
  )
}

cat(
  "- Synapse IDs foram registrados; versionNumber ainda deve ser preenchido\n"
)
cat(
  "  a partir das entidades correspondentes no Synapse antes do relatório final.\n\n"
)

cat(
  "============================================================\n"
)
cat(
  "00_audit_ROSMAP.R CONCLUÍDO\n"
)
cat(
  "============================================================\n"
)

sink()


# ==============================================================================
# 22. RESUMO NO TERMINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "AUDITORIA ROSMAP CONCLUÍDA\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Counts:\n"
)
cat(
  "  amostras totais = ",
  fmt_int(
    length(
      count_sample_ids
    )
  ),
  "\n",
  sep = ""
)
cat(
  "  DLPFC = ",
  fmt_int(
    length(
      dlpfc_ids
    )
  ),
  "\n",
  sep = ""
)
cat(
  "  features gênicas = ",
  fmt_int(
    n_gene_features
  ),
  "\n",
  sep = ""
)
cat(
  "  genes all-zero em DLPFC = ",
  fmt_int(
    n_all_zero_genes
  ),
  "\n\n",
  sep = ""
)

cat(
  "Mapeamento:\n"
)
cat(
  "  DLPFC → assay = ",
  fmt_int(
    sum(
      matched_assay
    )
  ),
  "/",
  fmt_int(
    length(
      dlpfc_ids
    )
  ),
  "\n",
  sep = ""
)
cat(
  "  DLPFC → biospecimen = ",
  fmt_int(
    sum(
      matched_biospecimen
    )
  ),
  "/",
  fmt_int(
    length(
      dlpfc_ids
    )
  ),
  "\n",
  sep = ""
)
cat(
  "  indivíduos únicos = ",
  fmt_int(
    n_unique_individuals
  ),
  "\n\n",
  sep = ""
)

cat(
  "Exclusões oficiais:\n"
)
cat(
  "  exclude=TRUE = ",
  fmt_int(
    n_excluded
  ),
  "\n",
  sep = ""
)
cat(
  "  DLPFC válidos = ",
  fmt_int(
    n_not_excluded
  ),
  "\n\n",
  sep = ""
)

cat(
  "Cenário principal cogdx 1 vs 4:\n"
)

print(
  strict_group_counts
)

cat(
  "  total = ",
  fmt_int(
    nrow(
      strict_cohort
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Cenário amplo cogdx 1 vs 4+5:\n"
)

print(
  broad_group_counts
)

cat(
  "  total = ",
  fmt_int(
    nrow(
      broad_cohort
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Confundimento completo group↔sequencingBatch: ",
  perfect_batch_separation,
  "\n\n",
  sep = ""
)

cat(
  "Relatório:\n  ",
  audit_report_file,
  "\n\n",
  sep = ""
)

cat(
  "Critérios pré-definidos:\n  ",
  predefined_criteria_out,
  "\n\n",
  sep = ""
)

cat(
  "Manifesto Synapse:\n  ",
  synapse_manifest_out,
  "\n\n",
  sep = ""
)

cat(
  "Próxima etapa, após revisar esta auditoria:\n"
)
cat(
  "  01_prepare_metadata_ROSMAP.R\n"
)
