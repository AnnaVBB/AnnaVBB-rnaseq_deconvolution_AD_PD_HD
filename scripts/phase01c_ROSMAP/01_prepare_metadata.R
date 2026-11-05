# ==============================================================================
# 01_prepare_metadata_ROSMAP.R
# ==============================================================================
#
# Fase 01c — ROSMAP
# Preparação do metadata analítico para a análise de expressão diferencial
#
# Este script deve ser executado APÓS:
#
#   scripts/phase01c_ROSMAP/00_audit_ROSMAP.R
#
# Objetivos
# ---------
# 1. Reconstruir, de forma reproduzível, o encadeamento:
#
#       counts
#         ↓
#       assay metadata
#         ↓
#       biospecimen metadata
#         ↓
#       clinical metadata
#
# 2. Preservar EXATAMENTE a ordem das amostras DLPFC na matriz de counts.
#
# 3. Aplicar apenas as exclusões e definições fenotípicas que foram
#    pré-especificadas ANTES do DESeq2:
#
#       a) remover amostras com exclude = TRUE no biospecimen metadata;
#
#       b) análise principal:
#            Control = cogdx 1
#            AD      = cogdx 4
#
#       c) sensibilidade fenotípica futura:
#            Control = cogdx 1
#            AD      = cogdx 4 ou 5
#
# 4. Preparar as covariáveis que poderão ser usadas posteriormente no
#    modelo ajustado de robustez:
#
#       sequencingBatch
#       msex
#       RIN
#       pmi
#
# 5. Preservar age_death SEM transformar "90+" em idade exata.
#
# 6. Produzir:
#
#       metadata_ROSMAP.tsv
#           metadata FINAL da análise principal (cogdx 1 vs 4)
#
#       metadata_ROSMAP_broad_sensitivity.tsv
#           metadata para sensibilidade fenotípica (cogdx 1 vs 4+5)
#
#       metadata_ROSMAP_all_DLPFC.tsv
#           todas as 210 amostras DLPFC com decisão documentada
#
#       metadata_ROSMAP_valid_DLPFC.tsv
#           DLPFC após exclusões oficiais
#
#       metadata_ROSMAP_exclusion_log.tsv
#           log explícito de inclusão/exclusão
#
#       metadata_ROSMAP_summary.tsv
#           resumo quantitativo da preparação
#
#       01_prepare_metadata_ROSMAP.txt
#           relatório textual da execução
#
# IMPORTANTE
# ----------
# Este script NÃO:
#   - executa DESeq2;
#   - filtra genes;
#   - transforma counts;
#   - converte Ensembl para GeneID;
#   - interpreta msex 0/1 como masculino/feminino sem codebook;
#   - converte age_death = "90+" para 90;
#   - usa Braak ou CERAD para redefinir diagnóstico;
#   - usa Braak/CERAD como covariáveis do modelo principal.
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

# ------------------------------------------------------------------------------
# Região alvo
# ------------------------------------------------------------------------------

TARGET_TISSUE <- "dorsolateral prefrontal cortex"
TARGET_SUFFIX <- "-DLPFC"


# ------------------------------------------------------------------------------
# Fenótipo principal
# ------------------------------------------------------------------------------

COGDX_CONTROL <- 1
COGDX_AD_MAIN <- 4


# ------------------------------------------------------------------------------
# Sensibilidade fenotípica futura
# ------------------------------------------------------------------------------

COGDX_AD_BROAD <- c(
  4,
  5
)


# ------------------------------------------------------------------------------
# Covariáveis pré-definidas para o modelo ajustado de robustez
#
# Modelo futuro planejado:
#
#   ~ sequencingBatch + msex + RIN + pmi + group
#
# age_death fica fora do modelo ajustado principal porque grande parte dos
# indivíduos está codificada como "90+", e portanto a idade exata é censurada.
# ------------------------------------------------------------------------------

ADJUSTED_MODEL_COVARIATES <- c(
  "sequencingBatch",
  "msex",
  "RIN",
  "pmi"
)


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

raw_dir <- file.path(
  "data",
  "raw",
  "ROSMAP"
)

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


# ------------------------------------------------------------------------------
# Resultados
# ------------------------------------------------------------------------------

candidate_dir <- file.path(
  "results",
  "candidate_ROSMAP"
)

metadata_dir <- file.path(
  candidate_dir,
  "metadata"
)

dir.create(
  metadata_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Arquivos de saída
# ------------------------------------------------------------------------------

metadata_main_out <- file.path(
  metadata_dir,
  "metadata_ROSMAP.tsv"
)

metadata_broad_out <- file.path(
  metadata_dir,
  "metadata_ROSMAP_broad_sensitivity.tsv"
)

metadata_all_dlpfc_out <- file.path(
  metadata_dir,
  "metadata_ROSMAP_all_DLPFC.tsv"
)

metadata_valid_dlpfc_out <- file.path(
  metadata_dir,
  "metadata_ROSMAP_valid_DLPFC.tsv"
)

exclusion_log_out <- file.path(
  metadata_dir,
  "metadata_ROSMAP_exclusion_log.tsv"
)

summary_out <- file.path(
  metadata_dir,
  "metadata_ROSMAP_summary.tsv"
)

batch_table_out <- file.path(
  metadata_dir,
  "metadata_ROSMAP_main_group_by_sequencingBatch.tsv"
)

sex_table_out <- file.path(
  metadata_dir,
  "metadata_ROSMAP_main_group_by_msex.tsv"
)

report_out <- file.path(
  metadata_dir,
  "01_prepare_metadata_ROSMAP.txt"
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
      is.character,
      logical(1)
    )
  ]

  for (col in char_cols) {

    set(
      dt,
      j = col,
      value = trimws(
        dt[[col]]
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
        as.numeric(x)
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
  digits = 2
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


# ------------------------------------------------------------------------------
# Construir crosstab sem alterar os dados.
# ------------------------------------------------------------------------------

make_crosstab <- function(
  dt,
  row_var,
  col_var
) {

  tmp <- dt[
    !is.na(
      get(row_var)
    ) &
      !is.na(
        get(col_var)
      )
  ]

  dcast(
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
}


# ------------------------------------------------------------------------------
# Reordenar metadata de acordo com a matriz de counts.
#
# A ordem da matriz é a referência canônica.
# ------------------------------------------------------------------------------

order_by_counts <- function(
  metadata,
  counts_sample_order,
  object_name
) {

  idx <- match(
    metadata$sample_id,
    counts_sample_order
  )

  if (anyNA(idx)) {

    stop(
      object_name,
      ": existem sample_id ausentes na matriz de counts."
    )
  }

  metadata[
    order(idx)
  ]
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
# 5. LER CABEÇALHO DOS COUNTS E DEFINIR A ORDEM CANÔNICA DLPFC
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "01 — PREPARAÇÃO DO METADATA ROSMAP\n"
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

count_names <- names(
  counts_header
)

if (length(count_names) < 2) {

  stop(
    "A matriz de counts possui menos de duas colunas."
  )
}

counts_sample_ids <- count_names[
  -1
]

dlpfc_sample_ids <- grep(
  paste0(
    TARGET_SUFFIX,
    "$"
  ),
  counts_sample_ids,
  value = TRUE
)

if (
  length(
    dlpfc_sample_ids
  ) == 0
) {

  stop(
    "Nenhuma coluna DLPFC foi identificada na matriz de counts."
  )
}

cat(
  "Amostras totais na matriz: ",
  fmt_int(
    length(
      counts_sample_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras DLPFC: ",
  fmt_int(
    length(
      dlpfc_sample_ids
    )
  ),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 6. LER METADADOS
# ==============================================================================

cat(
  "Lendo assay metadata...\n"
)

assay <- fread(
  assay_file,
  na.strings = c(
    "NA",
    ""
  )
)

cat(
  "Lendo biospecimen metadata...\n"
)

biospecimen <- fread(
  biospecimen_file,
  na.strings = c(
    "NA",
    ""
  )
)

cat(
  "Lendo clinical metadata...\n"
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


# ==============================================================================
# 7. VALIDAR COLUNAS NECESSÁRIAS
# ==============================================================================

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
# 8. RESTRINGIR ASSAY E BIOSPECIMEN ÀS 210 DLPFC
# ==============================================================================

assay_dlpfc <- assay[
  specimenID %in% dlpfc_sample_ids
]

biospecimen_dlpfc <- biospecimen[
  specimenID %in% dlpfc_sample_ids
]


# ------------------------------------------------------------------------------
# A auditoria 00 mostrou duplicatas globais no biospecimen, mas nenhuma entre
# os 210 specimenIDs DLPFC. Revalidamos isso aqui para o script ser independente.
# ------------------------------------------------------------------------------

if (
  anyDuplicated(
    assay_dlpfc$specimenID
  )
) {

  stop(
    "Existem specimenID duplicados no assay metadata entre as DLPFC."
  )
}

if (
  anyDuplicated(
    biospecimen_dlpfc$specimenID
  )
) {

  stop(
    "Existem specimenID duplicados no biospecimen metadata entre as DLPFC."
  )
}

if (
  anyDuplicated(
    clinical$individualID[
      !is.na(
        clinical$individualID
      )
    ]
  )
) {

  stop(
    "Existem individualID duplicados no clinical metadata."
  )
}


# ------------------------------------------------------------------------------
# Conferir cobertura completa dos 210 IDs.
# ------------------------------------------------------------------------------

missing_assay <- setdiff(
  dlpfc_sample_ids,
  assay_dlpfc$specimenID
)

missing_biospecimen <- setdiff(
  dlpfc_sample_ids,
  biospecimen_dlpfc$specimenID
)

if (length(missing_assay) > 0) {

  stop(
    "Há amostras DLPFC sem assay metadata: ",
    paste(
      missing_assay,
      collapse = ", "
    )
  )
}

if (length(missing_biospecimen) > 0) {

  stop(
    "Há amostras DLPFC sem biospecimen metadata: ",
    paste(
      missing_biospecimen,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 9. CONSTRUIR METADATA BASE PRESERVANDO A ORDEM DOS COUNTS
# ==============================================================================

# ------------------------------------------------------------------------------
# Em vez de usar merge() e depois tentar recuperar a ordem, usamos match().
#
# Isso garante:
#
#   metadata_base$sample_id
#
# na MESMA ordem das 210 colunas DLPFC da matriz de counts.
# ------------------------------------------------------------------------------

assay_idx <- match(
  dlpfc_sample_ids,
  assay_dlpfc$specimenID
)

biospecimen_idx <- match(
  dlpfc_sample_ids,
  biospecimen_dlpfc$specimenID
)

assay_aligned <- assay_dlpfc[
  assay_idx
]

biospecimen_aligned <- biospecimen_dlpfc[
  biospecimen_idx
]


# ------------------------------------------------------------------------------
# Confirmar que match() realmente preservou a correspondência 1:1.
# ------------------------------------------------------------------------------

stopifnot(
  identical(
    as.character(
      assay_aligned$specimenID
    ),
    as.character(
      dlpfc_sample_ids
    )
  )
)

stopifnot(
  identical(
    as.character(
      biospecimen_aligned$specimenID
    ),
    as.character(
      dlpfc_sample_ids
    )
  )
)


# ==============================================================================
# 10. MAPEAR INDIVÍDUO → CLINICAL METADATA
# ==============================================================================

clinical_idx <- match(
  biospecimen_aligned$individualID,
  clinical$individualID
)

if (anyNA(clinical_idx)) {

  missing_individuals <- unique(
    biospecimen_aligned$individualID[
      is.na(
        clinical_idx
      )
    ]
  )

  stop(
    "Existem indivíduos DLPFC sem clinical metadata: ",
    paste(
      missing_individuals,
      collapse = ", "
    )
  )
}

clinical_aligned <- clinical[
  clinical_idx
]


# ------------------------------------------------------------------------------
# Validar correspondência individualID.
# ------------------------------------------------------------------------------

stopifnot(
  identical(
    as.character(
      biospecimen_aligned$individualID
    ),
    as.character(
      clinical_aligned$individualID
    )
  )
)


# ==============================================================================
# 11. CONSTRUIR METADATA DLPFC COMPLETO
# ==============================================================================

metadata_all <- data.table(

  # Identificadores principais
  sample_id = as.character(
    dlpfc_sample_ids
  ),

  individual_id = as.character(
    biospecimen_aligned$individualID
  ),

  # Ordem original da coluna na matriz DLPFC.
  counts_order = seq_along(
    dlpfc_sample_ids
  ),

  # Biospecimen
  organ = as.character(
    biospecimen_aligned$organ
  ),

  tissue = as.character(
    biospecimen_aligned$tissue
  ),

  BrodmannArea = as.character(
    biospecimen_aligned$BrodmannArea
  ),

  exclude_raw = as.character(
    biospecimen_aligned$exclude
  ),

  excludeReason = as.character(
    biospecimen_aligned$excludeReason
  ),

  # Fenótipo clínico
  cogdx = safe_numeric(
    clinical_aligned$cogdx
  ),

  msex = safe_numeric(
    clinical_aligned$msex
  ),

  # age_death é preservado como texto.
  age_death_raw = as.character(
    clinical_aligned$age_death
  ),

  pmi = safe_numeric(
    clinical_aligned$pmi
  ),

  braaksc = safe_numeric(
    clinical_aligned$braaksc
  ),

  ceradsc = safe_numeric(
    clinical_aligned$ceradsc
  ),

  # Qualidade / batch do RNA-seq
  RIN = safe_numeric(
    assay_aligned$RIN
  ),

  sequencingBatch = as.character(
    assay_aligned$sequencingBatch
  )
)


# ==============================================================================
# 12. ADICIONAR CAMPOS OPCIONAIS ÚTEIS, SE EXISTIREM
# ==============================================================================
#
# Esses campos NÃO são necessários para o modelo principal, mas são preservados
# caso estejam disponíveis para auditorias e análises futuras.
# ==============================================================================

optional_clinical_fields <- c(
  "Study",
  "apoe_genotype",
  "age_at_visit_max",
  "age_first_ad_dx",
  "cts_mmse30_first_ad_dx",
  "cts_mmse30_lv",
  "dcfdx_lv"
)

for (field in optional_clinical_fields) {

  if (field %in% names(clinical_aligned)) {

    metadata_all[
      ,
      (field) := clinical_aligned[[field]]
    ]
  }
}


# ==============================================================================
# 13. NORMALIZAR FLAGS E CAMPOS DERIVADOS
# ==============================================================================

metadata_all[
  ,
  exclude_flag := as_exclude_flag(
    exclude_raw
  )
]


# ------------------------------------------------------------------------------
# Confirmar tissue.
# ------------------------------------------------------------------------------

metadata_all[
  ,
  tissue_confirmed_DLPFC :=
    tolower(
      trimws(
        tissue
      )
    ) ==
    tolower(
      TARGET_TISSUE
    )
]


# ------------------------------------------------------------------------------
# age_death
#
# IMPORTANTE:
#   "90+" NÃO é convertido para 90.
#
# Criamos:
#
#   age_death_90plus
#       TRUE/FALSE
#
#   age_death_exact
#       valor numérico apenas quando uma idade exata está disponível;
#       NA para "90+".
# ------------------------------------------------------------------------------

metadata_all[
  ,
  age_death_90plus :=
    trimws(
      age_death_raw
    ) == "90+"
]

metadata_all[
  ,
  age_death_exact :=
    fifelse(
      age_death_90plus,
      NA_real_,
      safe_numeric(
        age_death_raw
      )
    )
]


# ==============================================================================
# 14. DEFINIR DECISÃO ANALÍTICA PARA CADA UMA DAS 210 DLPFC
# ==============================================================================

metadata_all[
  ,
  analysis_decision := NA_character_
]


# ------------------------------------------------------------------------------
# 14.1 Exclusões oficiais primeiro
# ------------------------------------------------------------------------------

metadata_all[
  exclude_flag == TRUE,
  analysis_decision := "exclude_official"
]


# ------------------------------------------------------------------------------
# 14.2 Amostras válidas com cogdx principal
# ------------------------------------------------------------------------------

metadata_all[
  exclude_flag == FALSE &
    cogdx == COGDX_CONTROL,
  analysis_decision := "include_main_Control"
]

metadata_all[
  exclude_flag == FALSE &
    cogdx == COGDX_AD_MAIN,
  analysis_decision := "include_main_AD"
]


# ------------------------------------------------------------------------------
# 14.3 cogdx 5:
#
# Não entra na análise principal.
# É explicitamente reservado para sensibilidade fenotípica futura.
# ------------------------------------------------------------------------------

metadata_all[
  exclude_flag == FALSE &
    cogdx == 5,
  analysis_decision := "reserve_broad_sensitivity_AD"
]


# ------------------------------------------------------------------------------
# 14.4 Demais categorias cogdx
# ------------------------------------------------------------------------------

metadata_all[
  exclude_flag == FALSE &
    cogdx == 2,
  analysis_decision := "exclude_phenotype_cogdx_2"
]

metadata_all[
  exclude_flag == FALSE &
    cogdx == 3,
  analysis_decision := "exclude_phenotype_cogdx_3"
]

metadata_all[
  exclude_flag == FALSE &
    cogdx == 6,
  analysis_decision := "exclude_phenotype_cogdx_6"
]

metadata_all[
  exclude_flag == FALSE &
    is.na(
      cogdx
    ),
  analysis_decision := "exclude_missing_cogdx"
]


# ------------------------------------------------------------------------------
# Garantir que TODAS as 210 amostras receberam uma decisão.
# ------------------------------------------------------------------------------

if (
  anyNA(
    metadata_all$analysis_decision
  )
) {

  unknown <- metadata_all[
    is.na(
      analysis_decision
    )
  ]

  print(
    unknown
  )

  stop(
    "Há amostras DLPFC sem decisão analítica definida."
  )
}


# ==============================================================================
# 15. CONSTRUIR METADATA DA ANÁLISE PRINCIPAL
# ==============================================================================

metadata_main <- metadata_all[
  exclude_flag == FALSE &
    cogdx %in% c(
      COGDX_CONTROL,
      COGDX_AD_MAIN
    )
]


# ------------------------------------------------------------------------------
# Definir grupo explicitamente.
#
# Control será a referência futura no DESeq2.
# No TSV o campo é salvo como texto; no script 05 ele será convertido em factor
# com níveis c("Control", "AD").
# ------------------------------------------------------------------------------

metadata_main[
  ,
  group := fifelse(
    cogdx == COGDX_CONTROL,
    "Control",
    "AD"
  )
]


# ------------------------------------------------------------------------------
# Colocar group logo após individual_id para facilitar inspeção.
# ------------------------------------------------------------------------------

setcolorder(
  metadata_main,
  c(
    "sample_id",
    "individual_id",
    "group",
    "cogdx",
    setdiff(
      names(
        metadata_main
      ),
      c(
        "sample_id",
        "individual_id",
        "group",
        "cogdx"
      )
    )
  )
)


# ==============================================================================
# 16. CONSTRUIR METADATA DE SENSIBILIDADE FENOTÍPICA
# ==============================================================================

metadata_broad <- metadata_all[
  exclude_flag == FALSE &
    cogdx %in% c(
      COGDX_CONTROL,
      COGDX_AD_BROAD
    )
]

metadata_broad[
  ,
  group := fifelse(
    cogdx == COGDX_CONTROL,
    "Control",
    "AD"
  )
]

metadata_broad[
  ,
  phenotype_definition :=
    fifelse(
      cogdx == COGDX_CONTROL,
      "Control_cogdx1",
      fifelse(
        cogdx == 4,
        "AD_cogdx4",
        "AD_cogdx5"
      )
    )
]

setcolorder(
  metadata_broad,
  c(
    "sample_id",
    "individual_id",
    "group",
    "phenotype_definition",
    "cogdx",
    setdiff(
      names(
        metadata_broad
      ),
      c(
        "sample_id",
        "individual_id",
        "group",
        "phenotype_definition",
        "cogdx"
      )
    )
  )
)


# ==============================================================================
# 17. METADATA DLPFC VÁLIDO APÓS EXCLUSÕES OFICIAIS
# ==============================================================================

metadata_valid <- metadata_all[
  exclude_flag == FALSE
]


# ==============================================================================
# 18. VALIDAR REGIÃO
# ==============================================================================

if (
  any(
    !metadata_all$tissue_confirmed_DLPFC,
    na.rm = TRUE
  ) ||
    any(
      is.na(
        metadata_all$tissue_confirmed_DLPFC
      )
    )
) {

  stop(
    "Nem todas as 210 amostras foram confirmadas como DLPFC pelo campo tissue."
  )
}


# ==============================================================================
# 19. VALIDAR IDENTIFICADORES E TAMANHOS
# ==============================================================================

# ------------------------------------------------------------------------------
# 19.1 DLPFC completo
# ------------------------------------------------------------------------------

if (
  nrow(
    metadata_all
  ) != length(
    dlpfc_sample_ids
  )
) {

  stop(
    "metadata_all não contém exatamente todas as colunas DLPFC dos counts."
  )
}

if (
  anyDuplicated(
    metadata_all$sample_id
  )
) {

  stop(
    "metadata_all contém sample_id duplicado."
  )
}

if (
  anyDuplicated(
    metadata_all$individual_id
  )
) {

  stop(
    "metadata_all contém individual_id duplicado."
  )
}


# ------------------------------------------------------------------------------
# 19.2 Principal
# ------------------------------------------------------------------------------

if (
  anyDuplicated(
    metadata_main$sample_id
  )
) {

  stop(
    "metadata_main contém sample_id duplicado."
  )
}

if (
  anyDuplicated(
    metadata_main$individual_id
  )
) {

  stop(
    "metadata_main contém individual_id duplicado."
  )
}


# ------------------------------------------------------------------------------
# 19.3 Sensibilidade fenotípica
# ------------------------------------------------------------------------------

if (
  anyDuplicated(
    metadata_broad$sample_id
  )
) {

  stop(
    "metadata_broad contém sample_id duplicado."
  )
}


# ==============================================================================
# 20. GARANTIR ORDEM IDÊNTICA À MATRIZ DE COUNTS
# ==============================================================================

metadata_all <- order_by_counts(
  metadata_all,
  dlpfc_sample_ids,
  "metadata_all"
)

metadata_valid <- order_by_counts(
  metadata_valid,
  dlpfc_sample_ids,
  "metadata_valid"
)

metadata_main <- order_by_counts(
  metadata_main,
  dlpfc_sample_ids,
  "metadata_main"
)

metadata_broad <- order_by_counts(
  metadata_broad,
  dlpfc_sample_ids,
  "metadata_broad"
)


# ------------------------------------------------------------------------------
# Validar explicitamente que a ordem do metadata principal corresponde à ordem
# relativa dessas mesmas amostras no cabeçalho dos counts.
# ------------------------------------------------------------------------------

expected_main_order <- dlpfc_sample_ids[
  dlpfc_sample_ids %in% metadata_main$sample_id
]

if (
  !identical(
    as.character(
      metadata_main$sample_id
    ),
    as.character(
      expected_main_order
    )
  )
) {

  stop(
    "A ordem final do metadata principal não corresponde à ordem dos counts."
  )
}


# ==============================================================================
# 21. COMPLETUDE DAS COVARIÁVEIS DO MODELO AJUSTADO
# ==============================================================================

missing_adjusted_covariates <- setdiff(
  ADJUSTED_MODEL_COVARIATES,
  names(
    metadata_main
  )
)

if (
  length(
    missing_adjusted_covariates
  ) > 0
) {

  stop(
    "Covariáveis planejadas ausentes do metadata principal: ",
    paste(
      missing_adjusted_covariates,
      collapse = ", "
    )
  )
}

metadata_main[
  ,
  adjusted_model_complete :=
    complete.cases(
      .SD
    ),
  .SDcols = c(
    "group",
    ADJUSTED_MODEL_COVARIATES
  )
]

metadata_broad[
  ,
  adjusted_model_complete :=
    complete.cases(
      .SD
    ),
  .SDcols = c(
    "group",
    ADJUSTED_MODEL_COVARIATES
  )
]


# ==============================================================================
# 22. TABELAS DE RESUMO
# ==============================================================================

main_group_counts <- metadata_main[
  ,
  .N,
  by = group
][
  order(
    group
  )
]

broad_group_counts <- metadata_broad[
  ,
  .N,
  by = group
][
  order(
    group
  )
]

decision_counts <- metadata_all[
  ,
  .N,
  by = analysis_decision
][
  order(
    -N
  )
]


# ------------------------------------------------------------------------------
# Batch × group
# ------------------------------------------------------------------------------

main_batch_table <- make_crosstab(
  metadata_main,
  "sequencingBatch",
  "group"
)

fwrite(
  main_batch_table,
  batch_table_out,
  sep = "\t"
)


# ------------------------------------------------------------------------------
# msex × group
# ------------------------------------------------------------------------------

main_sex_table <- make_crosstab(
  metadata_main,
  "msex",
  "group"
)

fwrite(
  main_sex_table,
  sex_table_out,
  sep = "\t"
)


# ==============================================================================
# 23. LOG DE EXCLUSÃO / INCLUSÃO
# ==============================================================================

exclusion_log <- metadata_all[
  ,
  .(
    sample_id,
    individual_id,
    counts_order,
    cogdx,
    exclude_flag,
    excludeReason,
    analysis_decision
  )
]

fwrite(
  exclusion_log,
  exclusion_log_out,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 24. RESUMO QUANTITATIVO
# ==============================================================================

n_official_excluded <- metadata_all[
  exclude_flag == TRUE,
  .N
]

n_valid <- metadata_all[
  exclude_flag == FALSE,
  .N
]

n_main_control <- metadata_main[
  group == "Control",
  .N
]

n_main_ad <- metadata_main[
  group == "AD",
  .N
]

n_broad_control <- metadata_broad[
  group == "Control",
  .N
]

n_broad_ad <- metadata_broad[
  group == "AD",
  .N
]

n_main_adjusted_complete <- metadata_main[
  adjusted_model_complete == TRUE,
  .N
]

n_age_90plus_main <- metadata_main[
  age_death_90plus == TRUE,
  .N
]


metadata_summary <- data.table(
  metric = c(
    "DLPFC_in_counts",
    "DLPFC_confirmed_by_tissue",
    "official_excluded",
    "valid_DLPFC_after_official_exclusion",
    "main_Control_cogdx1",
    "main_AD_cogdx4",
    "main_total",
    "broad_Control_cogdx1",
    "broad_AD_cogdx4_or_5",
    "broad_total",
    "main_adjusted_model_complete_cases",
    "main_age_death_90plus"
  ),
  n = c(
    length(
      dlpfc_sample_ids
    ),
    sum(
      metadata_all$tissue_confirmed_DLPFC
    ),
    n_official_excluded,
    n_valid,
    n_main_control,
    n_main_ad,
    nrow(
      metadata_main
    ),
    n_broad_control,
    n_broad_ad,
    nrow(
      metadata_broad
    ),
    n_main_adjusted_complete,
    n_age_90plus_main
  )
)

fwrite(
  metadata_summary,
  summary_out,
  sep = "\t"
)


# ==============================================================================
# 25. SALVAR METADADOS
# ==============================================================================

fwrite(
  metadata_all,
  metadata_all_dlpfc_out,
  sep = "\t",
  na = "NA"
)

fwrite(
  metadata_valid,
  metadata_valid_dlpfc_out,
  sep = "\t",
  na = "NA"
)

fwrite(
  metadata_main,
  metadata_main_out,
  sep = "\t",
  na = "NA"
)

fwrite(
  metadata_broad,
  metadata_broad_out,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 26. RELATÓRIO TEXTUAL
# ==============================================================================

sink(
  report_out
)

cat(
  "============================================================\n"
)
cat(
  "FASE 01c — 01_prepare_metadata_ROSMAP.R\n"
)
cat(
  "PREPARAÇÃO DO METADATA ANALÍTICO\n"
)
cat(
  "============================================================\n\n"
)


# ------------------------------------------------------------------------------
# 1. Entrada
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "1. ENTRADA E REGIÃO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Amostras totais na matriz de counts: ",
  fmt_int(
    length(
      counts_sample_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras DLPFC pelo cabeçalho: ",
  fmt_int(
    length(
      dlpfc_sample_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "DLPFC confirmadas pelo campo tissue: ",
  fmt_int(
    sum(
      metadata_all$tissue_confirmed_DLPFC
    )
  ),
  " / ",
  fmt_int(
    nrow(
      metadata_all
    )
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 2. Exclusões oficiais
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "2. EXCLUSÕES OFICIAIS DO BIOSPECIMEN\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "DLPFC total: ",
  fmt_int(
    nrow(
      metadata_all
    )
  ),
  "\n",
  sep = ""
)

cat(
  "exclude=TRUE: ",
  fmt_int(
    n_official_excluded
  ),
  "\n",
  sep = ""
)

cat(
  "DLPFC após exclusões oficiais: ",
  fmt_int(
    n_valid
  ),
  "\n\n",
  sep = ""
)

official_reasons <- metadata_all[
  exclude_flag == TRUE,
  .N,
  by = excludeReason
][
  order(
    -N
  )
]

print(
  official_reasons
)

cat(
  "\n"
)


# ------------------------------------------------------------------------------
# 3. Decisões
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "3. DECISÃO ANALÍTICA DAS 210 DLPFC\n"
)
cat(
  "------------------------------------------------------------\n"
)

print(
  decision_counts
)

cat(
  "\n"
)


# ------------------------------------------------------------------------------
# 4. Principal
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "4. METADATA PRINCIPAL — cogdx 1 vs 4\n"
)
cat(
  "------------------------------------------------------------\n"
)

print(
  main_group_counts
)

cat(
  "\nTotal principal: ",
  fmt_int(
    nrow(
      metadata_main
    )
  ),
  "\n",
  sep = ""
)

cat(
  "sample_id únicos: ",
  fmt_int(
    uniqueN(
      metadata_main$sample_id
    )
  ),
  "\n",
  sep = ""
)

cat(
  "individual_id únicos: ",
  fmt_int(
    uniqueN(
      metadata_main$individual_id
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Ordem compatível com counts: TRUE\n\n"
)


# ------------------------------------------------------------------------------
# 5. Sensibilidade fenotípica
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "5. SENSIBILIDADE FENOTÍPICA — cogdx 1 vs 4+5\n"
)
cat(
  "------------------------------------------------------------\n"
)

print(
  broad_group_counts
)

cat(
  "\nTotal: ",
  fmt_int(
    nrow(
      metadata_broad
    )
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 6. Covariáveis
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "6. COVARIÁVEIS DO MODELO AJUSTADO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Modelo ajustado planejado:\n"
)
cat(
  "  ~ sequencingBatch + msex + RIN + pmi + group\n\n"
)

for (variable in ADJUSTED_MODEL_COVARIATES) {

  cat(
    variable,
    " — ausentes no metadata principal: ",
    fmt_int(
      sum(
        is.na(
          metadata_main[[variable]]
        )
      )
    ),
    "\n",
    sep = ""
  )
}

cat(
  "\nCasos completos para o modelo ajustado: ",
  fmt_int(
    n_main_adjusted_complete
  ),
  " / ",
  fmt_int(
    nrow(
      metadata_main
    )
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 7. Batch
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "7. sequencingBatch × group — ANÁLISE PRINCIPAL\n"
)
cat(
  "------------------------------------------------------------\n"
)

print(
  main_batch_table
)

cat(
  "\n"
)


# ------------------------------------------------------------------------------
# 8. msex
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "8. msex × group — ANÁLISE PRINCIPAL\n"
)
cat(
  "------------------------------------------------------------\n"
)

print(
  main_sex_table
)

cat(
  "\nA codificação 0/1 de msex foi PRESERVADA sem interpretação nominal.\n"
)
cat(
  "O significado dos níveis deve ser confirmado no codebook antes da\n"
)
cat(
  "descrição textual como masculino/feminino.\n\n"
)


# ------------------------------------------------------------------------------
# 9. age_death
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "9. age_death\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Amostras principais com age_death='90+': ",
  fmt_int(
    n_age_90plus_main
  ),
  " / ",
  fmt_int(
    nrow(
      metadata_main
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Esses valores NÃO foram convertidos para idade exata.\n"
)
cat(
  "age_death_exact = NA quando age_death_raw = '90+'.\n\n"
)


# ------------------------------------------------------------------------------
# 10. Braak / CERAD
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "10. BRAAK / CERAD\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "braaksc e ceradsc foram preservados no metadata para descrição e auditoria.\n"
)
cat(
  "Eles NÃO redefinem group e NÃO fazem parte do modelo principal de AD,\n"
)
cat(
  "pois representam características neuropatológicas associadas à doença.\n\n"
)


# ------------------------------------------------------------------------------
# 11. Arquivos
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "11. ARQUIVOS GERADOS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Principal:\n  ",
  metadata_main_out,
  "\n\n",
  sep = ""
)

cat(
  "Sensibilidade fenotípica:\n  ",
  metadata_broad_out,
  "\n\n",
  sep = ""
)

cat(
  "Todas DLPFC:\n  ",
  metadata_all_dlpfc_out,
  "\n\n",
  sep = ""
)

cat(
  "DLPFC válidas:\n  ",
  metadata_valid_dlpfc_out,
  "\n\n",
  sep = ""
)

cat(
  "Log de exclusão:\n  ",
  exclusion_log_out,
  "\n\n",
  sep = ""
)

cat(
  "Resumo:\n  ",
  summary_out,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 12. Conclusão
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "12. CONCLUSÃO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "- O metadata principal foi construído exclusivamente com DLPFC válidas.\n"
)

cat(
  "- As exclusões oficiais foram aplicadas antes da definição fenotípica.\n"
)

cat(
  "- O grupo principal permanece pré-definido como cogdx 1 vs 4.\n"
)

cat(
  "- cogdx 5 foi preservado somente para sensibilidade fenotípica futura.\n"
)

cat(
  "- A ordem das amostras foi preservada em relação à matriz de counts.\n"
)

cat(
  "- sequencingBatch, msex, RIN e pmi foram preservados para o modelo\n"
)
cat(
  "  ajustado de robustez.\n"
)

cat(
  "- age_death='90+' não foi transformado artificialmente em idade exata.\n"
)

cat(
  "- Braak e CERAD foram preservados apenas como características clínicas/\n"
)
cat(
  "  neuropatológicas e não redefinem o grupo analítico.\n\n"
)

cat(
  "============================================================\n"
)
cat(
  "01_prepare_metadata_ROSMAP.R CONCLUÍDO\n"
)
cat(
  "============================================================\n"
)

sink()


# ==============================================================================
# 27. RESUMO NO TERMINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "METADATA ROSMAP PREPARADO\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "DLPFC inicial: ",
  fmt_int(
    nrow(
      metadata_all
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Exclusões oficiais: ",
  fmt_int(
    n_official_excluded
  ),
  "\n",
  sep = ""
)

cat(
  "DLPFC válidas: ",
  fmt_int(
    n_valid
  ),
  "\n\n",
  sep = ""
)

cat(
  "Análise principal:\n"
)

print(
  main_group_counts
)

cat(
  "Total principal: ",
  fmt_int(
    nrow(
      metadata_main
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Sensibilidade fenotípica:\n"
)

print(
  broad_group_counts
)

cat(
  "Total sensibilidade: ",
  fmt_int(
    nrow(
      metadata_broad
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Casos completos para modelo ajustado: ",
  fmt_int(
    n_main_adjusted_complete
  ),
  "/",
  fmt_int(
    nrow(
      metadata_main
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Metadata principal:\n  ",
  metadata_main_out,
  "\n\n",
  sep = ""
)

cat(
  "Relatório:\n  ",
  report_out,
  "\n\n",
  sep = ""
)

cat(
  "Próxima etapa, após revisar este metadata:\n"
)
cat(
  "  02_prepare_counts_ROSMAP.R\n"
)
