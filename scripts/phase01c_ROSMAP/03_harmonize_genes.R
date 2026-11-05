# ==============================================================================
# 03_harmonize_genes_ROSMAP.R
# ==============================================================================
#
# Fase 01c — ROSMAP
# Harmonização dos identificadores gênicos:
#
#   Ensembl versionado/original
#          ↓
#   Ensembl stable ID
#          ↓
#   NCBI GeneID / Entrez Gene ID
#          ↓
#   matriz com GeneID único
#
# Este script deve ser executado APÓS:
#
#   00_audit_ROSMAP.R
#   01_prepare_metadata_ROSMAP.R
#   02_prepare_counts_ROSMAP.R
#
# Entrada principal:
#
#   results/candidate_ROSMAP/counts/counts_ROSMAP_raw_main.tsv.gz
#
# Metadata:
#
#   results/candidate_ROSMAP/metadata/metadata_ROSMAP.tsv
#
# Objetivos
# ---------
# 1. Preservar os identificadores originais do ROSMAP em uma tabela de auditoria.
#
# 2. Extrair o Ensembl stable ID dos identificadores originais.
#
#    Exemplos:
#
#       ENSG00000186092.6
#             ↓
#       ENSG00000186092
#
#    O procedimento também consegue reconhecer IDs Ensembl com sufixos
#    adicionais, desde que comecem por "ENSG" seguido de dígitos.
#
# 3. Mapear Ensembl stable ID → NCBI GeneID (Entrez Gene ID) usando:
#
#       AnnotationDbi
#       org.Hs.eg.db
#
# 4. Auditar explicitamente:
#
#       - IDs Ensembl reconhecidos;
#       - IDs não Ensembl;
#       - Ensembl sem mapeamento para GeneID;
#       - Ensembl com mapeamento ambíguo para >1 GeneID;
#       - vários Ensembl/original features convergindo para o mesmo GeneID.
#
# 5. Política de harmonização:
#
#       a) Ensembl → exatamente 1 GeneID:
#            manter;
#
#       b) Ensembl → nenhum GeneID:
#            não entra na matriz harmonizada;
#            registrar como unmapped;
#
#       c) Ensembl → >1 GeneID:
#            não duplicar counts;
#            excluir da matriz harmonizada e registrar como ambiguous;
#
#       d) vários features originais → mesmo GeneID:
#            somar counts por amostra;
#            registrar o colapso em tabela específica.
#
# 6. Produzir uma matriz final com:
#
#       GeneID + 144 amostras
#
#    contendo um único registro por NCBI GeneID.
#
# IMPORTANTE
# ----------
# Este script NÃO:
#   - filtra genes por baixa expressão;
#   - normaliza counts;
#   - executa DESeq2;
#   - arredonda counts;
#   - remove genes all-zero apenas por serem all-zero;
#   - escolhe genes por significância.
#
# A filtragem de baixa expressão será feita apenas no:
#
#   04_qc_counts_ROSMAP.R
#
# Critério planejado:
#
#   count >= 10 em pelo menos 67 amostras
#
# ==============================================================================


# ==============================================================================
# 0. PACOTES
# ==============================================================================

required_packages <- c(
  "data.table",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (
  length(
    missing_packages
  ) > 0
) {

  stop(
    paste0(
      "Pacotes necessários não instalados: ",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n\n",
      "Este script não instala pacotes automaticamente. ",
      "Instale/registre as dependências no ambiente do projeto antes de continuar."
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})


# ==============================================================================
# 1. CAMINHOS
# ==============================================================================

candidate_dir <- file.path(
  "results",
  "candidate_ROSMAP"
)

counts_dir <- file.path(
  candidate_dir,
  "counts"
)

metadata_dir <- file.path(
  candidate_dir,
  "metadata"
)

harmonization_dir <- file.path(
  candidate_dir,
  "harmonization"
)

dir.create(
  harmonization_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Entradas
# ------------------------------------------------------------------------------

counts_input <- file.path(
  counts_dir,
  "counts_ROSMAP_raw_main.tsv.gz"
)

metadata_input <- file.path(
  metadata_dir,
  "metadata_ROSMAP.tsv"
)


# ------------------------------------------------------------------------------
# Saídas
# ------------------------------------------------------------------------------

counts_output <- file.path(
  harmonization_dir,
  "counts_ROSMAP_GeneID_main.tsv.gz"
)

mapping_detailed_output <- file.path(
  harmonization_dir,
  "gene_mapping_ROSMAP_detailed.tsv"
)

mapping_summary_output <- file.path(
  harmonization_dir,
  "gene_mapping_ROSMAP_summary.tsv"
)

unmapped_output <- file.path(
  harmonization_dir,
  "gene_mapping_ROSMAP_unmapped.tsv"
)

ambiguous_output <- file.path(
  harmonization_dir,
  "gene_mapping_ROSMAP_ambiguous.tsv"
)

collapsed_output <- file.path(
  harmonization_dir,
  "gene_mapping_ROSMAP_collapsed_to_GeneID.tsv"
)

gene_annotation_output <- file.path(
  harmonization_dir,
  "gene_annotation_ROSMAP_GeneID.tsv"
)

annotation_metadata_output <- file.path(
  harmonization_dir,
  "annotation_database_metadata.tsv"
)

dimensions_output <- file.path(
  harmonization_dir,
  "harmonization_dimensions.tsv"
)

count_balance_output <- file.path(
  harmonization_dir,
  "harmonization_count_balance.tsv"
)

report_output <- file.path(
  harmonization_dir,
  "03_harmonize_genes_ROSMAP.txt"
)


# ==============================================================================
# 2. EXPECTATIVAS DA FASE 01c
# ==============================================================================

EXPECTED_MAIN_N <- 144L
EXPECTED_CONTROL_N <- 67L
EXPECTED_AD_N <- 77L


# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================

assert_file <- function(path) {

  if (
    !file.exists(
      path
    )
  ) {

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
    names(
      dt
    )
  )

  if (
    length(
      missing
    ) > 0
  ) {

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


fmt_int <- function(x) {

  format(
    round(
      as.numeric(x),
      digits = 0
    ),
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

fmt_num <- function(
  x,
  digits = 4
) {

  format(
    round(
      as.numeric(
        x
      ),
      digits
    ),
    nsmall = digits,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


collapse_unique <- function(
  x,
  sep = ";"
) {

  x <- unique(
    as.character(
      x[
        !is.na(
          x
        ) &
          as.character(
            x
          ) != ""
      ]
    )
  )

  if (
    length(
      x
    ) == 0
  ) {

    return(
      NA_character_
    )
  }

  paste(
    sort(
      x
    ),
    collapse = sep
  )
}


# ------------------------------------------------------------------------------
# Extrair stable Ensembl Gene ID.
#
# Exemplos reconhecidos:
#
#   ENSG00000186092.6
#       → ENSG00000186092
#
#   ENSG00000124333.15_PAR_Y
#       → ENSG00000124333
#
# Qualquer feature que NÃO comece por ENSG + dígitos recebe NA.
# ------------------------------------------------------------------------------

extract_ensembl_stable <- function(x) {

  x <- as.character(
    x
  )

  is_ensembl_like <- grepl(
    "^ENSG[0-9]+",
    x
  )

  out <- rep(
    NA_character_,
    length(
      x
    )
  )

  out[
    is_ensembl_like
  ] <- sub(
    "^(ENSG[0-9]+).*$",
    "\\1",
    x[
      is_ensembl_like
    ]
  )

  out
}


# ------------------------------------------------------------------------------
# Classificar formato original apenas para auditoria.
# ------------------------------------------------------------------------------

classify_feature_format <- function(x) {

  fifelse(
    grepl(
      "^ENSG[0-9]+\\.[0-9]+$",
      x
    ),
    "Ensembl_versioned_standard",
    fifelse(
      grepl(
        "^ENSG[0-9]+$",
        x
      ),
      "Ensembl_unversioned_standard",
      fifelse(
        grepl(
          "^ENSG[0-9]+\\.[0-9]+_PAR_Y$",
          x
        ),
        "Ensembl_versioned_PAR_Y",
        fifelse(
          grepl(
            "^ENSG[0-9]+",
            x
          ),
          "Ensembl_other_format",
          "non_Ensembl_format"
        )
      )
    )
  )
}


# ------------------------------------------------------------------------------
# Contar valores fracionários.
# ------------------------------------------------------------------------------

count_fractional_values <- function(x) {

  sum(
    abs(
      x -
        round(
          x
        )
    ) > 1e-8,
    na.rm = TRUE
  )
}


# ------------------------------------------------------------------------------
# Salvar TSV comprimido sem depender de zlib compilado no data.table.
#
# Estratégia:
#   1. fwrite() sem compressão
#   2. gzip do sistema
# ------------------------------------------------------------------------------

write_tsv_gzip <- function(
  dt,
  path_gz
) {

  if (
    !grepl(
      "\\.gz$",
      path_gz
    )
  ) {

    stop(
      "write_tsv_gzip exige caminho terminado em .gz."
    )
  }

  path_plain <- sub(
    "\\.gz$",
    "",
    path_gz
  )

  fwrite(
    dt,
    path_plain,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  gzip_bin <- Sys.which(
    "gzip"
  )

  if (
    identical(
      gzip_bin,
      ""
    )
  ) {

    stop(
      paste0(
        "O executável gzip não foi encontrado. ",
        "Arquivo sem compressão salvo em:\n",
        path_plain
      )
    )
  }

  status <- system2(
    command = gzip_bin,
    args = c(
      "-f",
      path_plain
    )
  )

  if (
    status != 0
  ) {

    stop(
      "Falha ao comprimir arquivo com gzip. Código: ",
      status
    )
  }

  if (
    !file.exists(
      path_gz
    )
  ) {

    stop(
      "Arquivo .gz esperado não foi criado: ",
      path_gz
    )
  }

  invisible(
    path_gz
  )
}


# ==============================================================================
# 4. VALIDAR ENTRADAS
# ==============================================================================

assert_file(
  counts_input
)

assert_file(
  metadata_input
)


# ==============================================================================
# 5. LER METADATA PRINCIPAL
# ==============================================================================

cat(
  "\n============================================================\n"
)
cat(
  "03 — HARMONIZAÇÃO GÊNICA ROSMAP\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Lendo metadata principal...\n"
)

metadata <- fread(
  metadata_input,
  na.strings = c(
    "NA",
    ""
  )
)

assert_columns(
  metadata,
  c(
    "sample_id",
    "individual_id",
    "group"
  ),
  "metadata_ROSMAP.tsv"
)

metadata[
  ,
  sample_id := trimws(
    as.character(
      sample_id
    )
  )
]

metadata[
  ,
  group := trimws(
    as.character(
      group
    )
  )
]


# ------------------------------------------------------------------------------
# Validar tamanho e grupos.
# ------------------------------------------------------------------------------

if (
  nrow(
    metadata
  ) != EXPECTED_MAIN_N
) {

  stop(
    "Metadata principal deveria conter ",
    EXPECTED_MAIN_N,
    " amostras; observado: ",
    nrow(
      metadata
    ),
    "."
  )
}

if (
  anyDuplicated(
    metadata$sample_id
  )
) {

  stop(
    "metadata principal contém sample_id duplicado."
  )
}

group_counts <- metadata[
  ,
  .N,
  by = group
]

n_ad <- group_counts[
  group == "AD",
  N
]

n_control <- group_counts[
  group == "Control",
  N
]

if (
  length(
    n_ad
  ) != 1 ||
    n_ad != EXPECTED_AD_N
) {

  stop(
    "Número inesperado de casos AD no metadata."
  )
}

if (
  length(
    n_control
  ) != 1 ||
    n_control != EXPECTED_CONTROL_N
) {

  stop(
    "Número inesperado de controles no metadata."
  )
}


# ==============================================================================
# 6. LER MATRIZ PREPARADA DO 02
# ==============================================================================

cat(
  "Lendo matriz preparada do 02_prepare_counts...\n"
)

counts <- fread(
  counts_input,
  na.strings = c(
    "NA",
    ""
  )
)

assert_columns(
  counts,
  "feature",
  "counts_ROSMAP_raw_main.tsv.gz"
)

counts[
  ,
  feature := trimws(
    as.character(
      feature
    )
  )
]


# ------------------------------------------------------------------------------
# Validar sample IDs e ordem.
# ------------------------------------------------------------------------------

count_sample_ids <- names(
  counts
)[
  -1
]

if (
  length(
    count_sample_ids
  ) != EXPECTED_MAIN_N
) {

  stop(
    "A matriz preparada não possui 144 colunas de amostras."
  )
}

if (
  !identical(
    count_sample_ids,
    metadata$sample_id
  )
) {

  stop(
    paste0(
      "A ordem das amostras na matriz de counts não corresponde exatamente ",
      "à ordem do metadata principal."
    )
  )
}

if (
  anyDuplicated(
    counts$feature
  )
) {

  stop(
    "Existem feature IDs duplicados antes da harmonização."
  )
}

if (
  anyNA(
    counts$feature
  ) ||
    any(
      counts$feature == ""
    )
) {

  stop(
    "Existem feature IDs ausentes/vazios antes da harmonização."
  )
}

n_features_input <- nrow(
  counts
)

cat(
  "Matriz de entrada: ",
  fmt_int(
    n_features_input
  ),
  " features × ",
  fmt_int(
    length(
      count_sample_ids
    )
  ),
  " amostras.\n\n",
  sep = ""
)


# ==============================================================================
# 7. VALIDAR INTEGRIDADE NUMÉRICA ANTES DA HARMONIZAÇÃO
# ==============================================================================

count_matrix <- as.matrix(
  counts[
    ,
    -1,
    with = FALSE
  ]
)

storage.mode(
  count_matrix
) <- "numeric"

n_na_input <- sum(
  is.na(
    count_matrix
  )
)

n_nonfinite_input <- sum(
  !is.finite(
    count_matrix
  ),
  na.rm = TRUE
)

n_negative_input <- sum(
  count_matrix < 0,
  na.rm = TRUE
)

n_fractional_input <- count_fractional_values(
  count_matrix
)

if (
  n_na_input > 0
) {

  stop(
    "Foram encontrados NA nos counts de entrada."
  )
}

if (
  n_nonfinite_input > 0
) {

  stop(
    "Foram encontrados valores não finitos nos counts de entrada."
  )
}

if (
  n_negative_input > 0
) {

  stop(
    "Foram encontrados counts negativos."
  )
}

if (
  n_fractional_input > 0
) {

  stop(
    "Foram encontrados counts fracionários. A harmonização foi interrompida."
  )
}


# ==============================================================================
# 8. CONSTRUIR TABELA DE FEATURES ORIGINAIS
# ==============================================================================

feature_map <- data.table(
  row_index = seq_len(
    nrow(
      counts
    )
  ),
  original_feature = counts$feature
)

feature_map[
  ,
  feature_format := classify_feature_format(
    original_feature
  )
]

feature_map[
  ,
  ensembl_id := extract_ensembl_stable(
    original_feature
  )
]


# ------------------------------------------------------------------------------
# Auditoria de normalização.
# ------------------------------------------------------------------------------

n_ensembl_extractable <- sum(
  !is.na(
    feature_map$ensembl_id
  )
)

n_non_ensembl <- sum(
  is.na(
    feature_map$ensembl_id
  )
)

format_summary <- feature_map[
  ,
  .N,
  by = feature_format
][
  order(
    -N
  )
]

cat(
  "Ensembl stable ID extraído de ",
  fmt_int(
    n_ensembl_extractable
  ),
  " / ",
  fmt_int(
    n_features_input
  ),
  " features.\n",
  sep = ""
)

cat(
  "Features sem Ensembl reconhecível: ",
  fmt_int(
    n_non_ensembl
  ),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 9. AUDITAR DUPLICAÇÃO APÓS REMOVER VERSÃO/SUFIXO
# ==============================================================================

stable_id_multiplicity <- feature_map[
  !is.na(
    ensembl_id
  ),
  .(
    n_original_features = .N,
    original_features = collapse_unique(
      original_feature
    )
  ),
  by = ensembl_id
]

n_stable_ids_multiple_features <- stable_id_multiplicity[
  n_original_features > 1,
  .N
]

n_features_in_stable_duplicates <- stable_id_multiplicity[
  n_original_features > 1,
  sum(
    n_original_features
  )
]

cat(
  "Stable Ensembl IDs associados a >1 feature original: ",
  fmt_int(
    n_stable_ids_multiple_features
  ),
  "\n",
  sep = ""
)


# ==============================================================================
# 10. CONSULTAR org.Hs.eg.db
# ==============================================================================

cat(
  "Consultando org.Hs.eg.db para Ensembl → GeneID...\n"
)

ensembl_query <- sort(
  unique(
    feature_map$ensembl_id[
      !is.na(
        feature_map$ensembl_id
      )
    ]
  )
)

db_ensembl_keys <- AnnotationDbi::keys(
  org.Hs.eg.db,
  keytype = "ENSEMBL"
)

ensembl_in_db <- intersect(
  ensembl_query,
  db_ensembl_keys
)

ensembl_not_in_db <- setdiff(
  ensembl_query,
  db_ensembl_keys
)

cat(
  "Stable Ensembl únicos consultados: ",
  fmt_int(
    length(
      ensembl_query
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Presentes como chave ENSEMBL em org.Hs.eg.db: ",
  fmt_int(
    length(
      ensembl_in_db
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Ausentes como chave ENSEMBL em org.Hs.eg.db: ",
  fmt_int(
    length(
      ensembl_not_in_db
    )
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# Consultar apenas chaves válidas.
# ------------------------------------------------------------------------------

if (
  length(
    ensembl_in_db
  ) == 0
) {

  stop(
    "Nenhum Ensembl ID da matriz foi encontrado em org.Hs.eg.db."
  )
}

mapping_raw <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ensembl_in_db,
    keytype = "ENSEMBL",
    columns = c(
      "ENTREZID"
    )
  )
)

mapping_raw <- as.data.table(
  mapping_raw
)

setnames(
  mapping_raw,
  old = c(
    "ENSEMBL",
    "ENTREZID"
  ),
  new = c(
    "ensembl_id",
    "GeneID"
  )
)

mapping_raw[
  ,
  ensembl_id := as.character(
    ensembl_id
  )
]

mapping_raw[
  ,
  GeneID := as.character(
    GeneID
  )
]


# ==============================================================================
# 11. RESUMIR Ensembl → GeneID
# ==============================================================================

mapping_by_ensembl <- mapping_raw[
  !is.na(
    GeneID
  ) &
    GeneID != "",
  .(
    n_GeneID = uniqueN(
      GeneID
    ),
    GeneID = collapse_unique(
      GeneID
    )
  ),
  by = ensembl_id
]


# ------------------------------------------------------------------------------
# Adicionar todos os Ensembl consultados, inclusive os sem GeneID retornado.
# ------------------------------------------------------------------------------

ensembl_status <- data.table(
  ensembl_id = ensembl_query
)

ensembl_status <- merge(
  ensembl_status,
  mapping_by_ensembl,
  by = "ensembl_id",
  all.x = TRUE,
  sort = FALSE
)

ensembl_status[
  is.na(
    n_GeneID
  ),
  n_GeneID := 0L
]

ensembl_status[
  ,
  mapping_status := fifelse(
    n_GeneID == 0,
    "unmapped",
    fifelse(
      n_GeneID == 1,
      "mapped_unique",
      "ambiguous_multiple_GeneID"
    )
  )
]


# ------------------------------------------------------------------------------
# Para mapeamento único, GeneID deve ser um único valor.
# Para ambiguidades, GeneID permanece como lista textual separada por ";".
# ------------------------------------------------------------------------------

n_ensembl_unique_mapped <- ensembl_status[
  mapping_status == "mapped_unique",
  .N
]

n_ensembl_unmapped <- ensembl_status[
  mapping_status == "unmapped",
  .N
]

n_ensembl_ambiguous <- ensembl_status[
  mapping_status == "ambiguous_multiple_GeneID",
  .N
]


# ==============================================================================
# 12. ANEXAR MAPEAMENTO A CADA FEATURE ORIGINAL
# ==============================================================================

feature_map <- merge(
  feature_map,
  ensembl_status,
  by = "ensembl_id",
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------------------------
# O merge pode alterar ordem; restauramos pela row_index.
# ------------------------------------------------------------------------------

setorder(
  feature_map,
  row_index
)


# ------------------------------------------------------------------------------
# Features não Ensembl recebem status próprio.
# ------------------------------------------------------------------------------

feature_map[
  is.na(
    ensembl_id
  ),
  `:=`(
    n_GeneID = 0L,
    GeneID = NA_character_,
    mapping_status = "non_Ensembl_format"
  )
]


# ------------------------------------------------------------------------------
# Garantir que mapeamento único não contenha ";"
# ------------------------------------------------------------------------------

if (
  any(
    feature_map$mapping_status == "mapped_unique" &
      grepl(
        ";",
        feature_map$GeneID,
        fixed = TRUE
      ),
    na.rm = TRUE
  )
) {

  stop(
    "Inconsistência interna: mapping_status='mapped_unique' com múltiplos GeneID."
  )
}


# ==============================================================================
# 13. TABELAS DE UNMAPPED E AMBIGUOUS
# ==============================================================================

unmapped_features <- feature_map[
  mapping_status %in% c(
    "unmapped",
    "non_Ensembl_format"
  )
]

ambiguous_features <- feature_map[
  mapping_status == "ambiguous_multiple_GeneID"
]

fwrite(
  unmapped_features,
  unmapped_output,
  sep = "\t",
  na = "NA"
)

fwrite(
  ambiguous_features,
  ambiguous_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 14. PREPARAR FEATURES COM MAPEAMENTO UNÍVOCO
# ==============================================================================

keep_unique_mapping <- feature_map$mapping_status == "mapped_unique"

if (
  !any(
    keep_unique_mapping
  )
) {

  stop(
    "Nenhuma feature possui mapeamento unívoco Ensembl → GeneID."
  )
}

mapped_features <- feature_map[
  keep_unique_mapping
]


# ------------------------------------------------------------------------------
# Validar GeneID.
# ------------------------------------------------------------------------------

if (
  anyNA(
    mapped_features$GeneID
  ) ||
    any(
      mapped_features$GeneID == ""
    )
) {

  stop(
    "Há GeneID ausente entre os mapeamentos classificados como únicos."
  )
}


# ==============================================================================
# 15. AUDITAR MUITOS Ensembl/FEATURES → UM GeneID
# ==============================================================================

geneid_multiplicity <- mapped_features[
  ,
  .(
    n_original_features = .N,
    n_stable_ensembl = uniqueN(
      ensembl_id
    ),
    original_features = collapse_unique(
      original_feature
    ),
    ensembl_ids = collapse_unique(
      ensembl_id
    )
  ),
  by = GeneID
]

collapsed_geneids <- geneid_multiplicity[
  n_original_features > 1
]

fwrite(
  collapsed_geneids,
  collapsed_output,
  sep = "\t",
  na = "NA"
)

n_geneids_collapsed <- nrow(
  collapsed_geneids
)

n_input_features_in_collapses <- collapsed_geneids[
  ,
  sum(
    n_original_features
  )
]


# ==============================================================================
# 16. COLAPSAR COUNTS PARA GeneID ÚNICO
# ==============================================================================
#
# Cada feature com mapeamento unívoco possui exatamente um GeneID.
#
# Se mais de uma feature converge para o mesmo GeneID:
#
#   counts_GeneID = soma das features correspondentes, por amostra
#
# rowsum() evita converter a matriz para formato longo.
# ==============================================================================

cat(
  "Colapsando counts por GeneID...\n"
)

mapped_row_indices <- mapped_features$row_index

mapped_count_matrix <- count_matrix[
  mapped_row_indices,
  ,
  drop = FALSE
]

mapped_geneids <- mapped_features$GeneID


# ------------------------------------------------------------------------------
# Total de counts que deverá ser preservado após o colapso.
# ------------------------------------------------------------------------------

total_counts_input <- sum(
  count_matrix
)

total_counts_retained_before_collapse <- sum(
  mapped_count_matrix
)

counts_harmonized_matrix <- rowsum(
  mapped_count_matrix,
  group = mapped_geneids,
  reorder = FALSE
)


# ------------------------------------------------------------------------------
# rowsum usa os grupos como rownames.
# ------------------------------------------------------------------------------

geneids_harmonized <- rownames(
  counts_harmonized_matrix
)

if (
  anyDuplicated(
    geneids_harmonized
  )
) {

  stop(
    "GeneID duplicado permaneceu após rowsum()."
  )
}


# ------------------------------------------------------------------------------
# Ordenar por GeneID numérico para saída determinística.
# ------------------------------------------------------------------------------

geneid_numeric_order <- suppressWarnings(
  as.numeric(
    geneids_harmonized
  )
)

if (
  anyNA(
    geneid_numeric_order
  )
) {

  stop(
    "Foram encontrados GeneID não numéricos após harmonização."
  )
}

ord <- order(
  geneid_numeric_order
)

counts_harmonized_matrix <- counts_harmonized_matrix[
  ord,
  ,
  drop = FALSE
]

geneids_harmonized <- rownames(
  counts_harmonized_matrix
)


# ==============================================================================
# 17. VALIDAR BALANÇO DE COUNTS
# ==============================================================================

total_counts_output <- sum(
  counts_harmonized_matrix
)

counts_excluded_from_harmonized <- (
  total_counts_input -
    total_counts_retained_before_collapse
)

if (
  !isTRUE(
    all.equal(
      total_counts_retained_before_collapse,
      total_counts_output,
      tolerance = 0
    )
  )
) {

  stop(
    paste0(
      "A soma dos counts mudou durante o colapso GeneID. ",
      "Isso não deveria acontecer."
    )
  )
}

percent_counts_retained <- (
  100 *
    total_counts_output /
    total_counts_input
)

percent_counts_excluded <- (
  100 -
    percent_counts_retained
)


# ==============================================================================
# 18. VALIDAR INTEGRIDADE DA MATRIZ HARMONIZADA
# ==============================================================================

n_na_output <- sum(
  is.na(
    counts_harmonized_matrix
  )
)

n_nonfinite_output <- sum(
  !is.finite(
    counts_harmonized_matrix
  ),
  na.rm = TRUE
)

n_negative_output <- sum(
  counts_harmonized_matrix < 0,
  na.rm = TRUE
)

n_fractional_output <- count_fractional_values(
  counts_harmonized_matrix
)

n_all_zero_output <- sum(
  rowSums(
    counts_harmonized_matrix
  ) == 0
)

n_all_zero_samples_output <- sum(
  colSums(
    counts_harmonized_matrix
  ) == 0
)

if (
  n_na_output > 0 ||
    n_nonfinite_output > 0 ||
    n_negative_output > 0 ||
    n_fractional_output > 0
) {

  stop(
    "A matriz harmonizada apresentou problema de integridade numérica."
  )
}

if (
  n_all_zero_samples_output > 0
) {

  stop(
    "A matriz harmonizada contém amostra all-zero."
  )
}

if (
  !identical(
    colnames(
      counts_harmonized_matrix
    ),
    metadata$sample_id
  )
) {

  stop(
    "A ordem das amostras mudou durante a harmonização."
  )
}


# ==============================================================================
# 19. CRIAR ANOTAÇÃO GeneID → SYMBOL / GENENAME
# ==============================================================================
#
# A anotação é descritiva.
# A chave analítica da matriz continua sendo GeneID.
# ==============================================================================

cat(
  "Anotando GeneID com SYMBOL/GENENAME...\n"
)

geneid_keys <- as.character(
  geneids_harmonized
)

gene_annotation_raw <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = geneid_keys,
    keytype = "ENTREZID",
    columns = c(
      "SYMBOL",
      "GENENAME",
      "ENSEMBL"
    )
  )
)

gene_annotation_raw <- as.data.table(
  gene_annotation_raw
)

setnames(
  gene_annotation_raw,
  "ENTREZID",
  "GeneID"
)

gene_annotation <- gene_annotation_raw[
  ,
  .(
    symbol = collapse_unique(
      SYMBOL
    ),
    gene_name = collapse_unique(
      GENENAME
    ),
    ensembl_annotation = collapse_unique(
      ENSEMBL
    )
  ),
  by = GeneID
]


# ------------------------------------------------------------------------------
# Garantir que todos os GeneID da matriz apareçam na tabela de anotação.
# ------------------------------------------------------------------------------

gene_annotation <- merge(
  data.table(
    GeneID = geneid_keys
  ),
  gene_annotation,
  by = "GeneID",
  all.x = TRUE,
  sort = FALSE
)

gene_annotation[
  ,
  GeneID_numeric := suppressWarnings(
    as.numeric(
      GeneID
    )
  )
]

setorder(
  gene_annotation,
  GeneID_numeric
)

gene_annotation[
  ,
  GeneID_numeric := NULL
]

fwrite(
  gene_annotation,
  gene_annotation_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 20. SALVAR METADADOS DA BASE DE ANOTAÇÃO
# ==============================================================================

annotation_db_metadata <- as.data.table(
  AnnotationDbi::metadata(
    org.Hs.eg.db
  )
)

annotation_package_metadata <- data.table(
  name = c(
    "R_version",
    "AnnotationDbi_version",
    "org.Hs.eg.db_version"
  ),
  value = c(
    as.character(
      getRversion()
    ),
    as.character(
      packageVersion(
        "AnnotationDbi"
      )
    ),
    as.character(
      packageVersion(
        "org.Hs.eg.db"
      )
    )
  )
)

annotation_metadata_combined <- rbind(
  annotation_package_metadata,
  annotation_db_metadata[
    ,
    .(
      name,
      value
    )
  ],
  fill = TRUE
)

fwrite(
  annotation_metadata_combined,
  annotation_metadata_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 21. COMPLETAR TABELA DETALHADA DE MAPEAMENTO
# ==============================================================================

# ------------------------------------------------------------------------------
# Informar quantas features originais convergem para o GeneID final.
# ------------------------------------------------------------------------------

geneid_n_original <- mapped_features[
  ,
  .(
    n_original_features_to_GeneID = .N,
    n_stable_ensembl_to_GeneID = uniqueN(
      ensembl_id
    )
  ),
  by = GeneID
]

feature_map <- merge(
  feature_map,
  geneid_n_original,
  by = "GeneID",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  feature_map,
  row_index
)

feature_map[
  ,
  retained_in_harmonized_matrix :=
    mapping_status == "mapped_unique"
]

fwrite(
  feature_map,
  mapping_detailed_output,
  sep = "\t",
  na = "NA"
)


# ==============================================================================
# 22. RESUMO DE MAPEAMENTO
# ==============================================================================

mapping_summary <- data.table(
  metric = c(
    "input_gene_features",
    "Ensembl_extractable_features",
    "non_Ensembl_features",
    "unique_stable_Ensembl_ids",
    "stable_Ensembl_ids_with_multiple_original_features",
    "stable_Ensembl_unique_mapped_to_one_GeneID",
    "stable_Ensembl_unmapped",
    "stable_Ensembl_ambiguous_multiple_GeneID",
    "original_features_retained_unique_mapping",
    "original_features_excluded_unmapped_or_non_Ensembl",
    "original_features_excluded_ambiguous",
    "unique_GeneID_before_collapse",
    "GeneID_with_multiple_original_features",
    "final_unique_GeneID",
    "all_zero_GeneID_after_harmonization"
  ),
  value = c(
    n_features_input,
    n_ensembl_extractable,
    n_non_ensembl,
    length(
      ensembl_query
    ),
    n_stable_ids_multiple_features,
    n_ensembl_unique_mapped,
    n_ensembl_unmapped,
    n_ensembl_ambiguous,
    sum(
      feature_map$mapping_status == "mapped_unique",
      na.rm = TRUE
    ),
    sum(
      feature_map$mapping_status %in% c(
        "unmapped",
        "non_Ensembl_format"
      ),
      na.rm = TRUE
    ),
    sum(
      feature_map$mapping_status == "ambiguous_multiple_GeneID",
      na.rm = TRUE
    ),
    uniqueN(
      mapped_features$GeneID
    ),
    n_geneids_collapsed,
    nrow(
      counts_harmonized_matrix
    ),
    n_all_zero_output
  )
)

fwrite(
  mapping_summary,
  mapping_summary_output,
  sep = "\t"
)


# ==============================================================================
# 23. BALANÇO DE COUNTS
# ==============================================================================

count_balance <- data.table(
  metric = c(
    "total_counts_input_before_mapping",
    "total_counts_retained_unique_mapping_before_collapse",
    "total_counts_output_after_GeneID_collapse",
    "total_counts_excluded_unmapped_or_ambiguous",
    "percent_counts_retained",
    "percent_counts_excluded"
  ),
  value = c(
    total_counts_input,
    total_counts_retained_before_collapse,
    total_counts_output,
    counts_excluded_from_harmonized,
    percent_counts_retained,
    percent_counts_excluded
  )
)

fwrite(
  count_balance,
  count_balance_output,
  sep = "\t"
)


# ==============================================================================
# 24. DIMENSÕES
# ==============================================================================

dimensions <- data.table(
  object = c(
    "input_features",
    "input_samples",
    "features_with_unique_mapping",
    "final_unique_GeneID",
    "final_samples"
  ),
  n = c(
    n_features_input,
    ncol(
      count_matrix
    ),
    nrow(
      mapped_features
    ),
    nrow(
      counts_harmonized_matrix
    ),
    ncol(
      counts_harmonized_matrix
    )
  )
)

fwrite(
  dimensions,
  dimensions_output,
  sep = "\t"
)


# ==============================================================================
# 25. CONSTRUIR E SALVAR MATRIZ GeneID
# ==============================================================================

counts_harmonized <- as.data.table(
  counts_harmonized_matrix,
  keep.rownames = "GeneID"
)

counts_harmonized[
  ,
  GeneID := as.character(
    GeneID
  )
]


# ------------------------------------------------------------------------------
# Validar unicidade final.
# ------------------------------------------------------------------------------

if (
  anyDuplicated(
    counts_harmonized$GeneID
  )
) {

  stop(
    "A matriz final contém GeneID duplicado."
  )
}

if (
  ncol(
    counts_harmonized
  ) != (
    EXPECTED_MAIN_N + 1L
  )
) {

  stop(
    "Número inesperado de colunas na matriz harmonizada."
  )
}

if (
  !identical(
    names(
      counts_harmonized
    )[
      -1
    ],
    metadata$sample_id
  )
) {

  stop(
    "Validação final falhou: samples da matriz GeneID != metadata."
  )
}


# ------------------------------------------------------------------------------
# Salvar comprimido.
# ------------------------------------------------------------------------------

cat(
  "Salvando matriz harmonizada por GeneID...\n"
)

write_tsv_gzip(
  counts_harmonized,
  counts_output
)


# ==============================================================================
# 26. RELATÓRIO TEXTUAL
# ==============================================================================

sink(
  report_output
)

cat(
  "============================================================\n"
)
cat(
  "FASE 01c — 03_harmonize_genes_ROSMAP.R\n"
)
cat(
  "HARMONIZAÇÃO Ensembl → NCBI GeneID\n"
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
  "1. ENTRADA\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Matriz de entrada:\n  ",
  counts_input,
  "\n\n",
  sep = ""
)

cat(
  "Features de entrada: ",
  fmt_int(
    n_features_input
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras: ",
  fmt_int(
    ncol(
      count_matrix
    )
  ),
  "\n",
  sep = ""
)

cat(
  "AD: ",
  fmt_int(
    n_ad
  ),
  "\n",
  sep = ""
)

cat(
  "Control: ",
  fmt_int(
    n_control
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 2. Fonte da anotação
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "2. FONTE DE ANOTAÇÃO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Mapeamento utilizado:\n"
)
cat(
  "  Ensembl stable gene ID → NCBI GeneID/Entrez Gene ID\n\n"
)

cat(
  "AnnotationDbi version: ",
  as.character(
    packageVersion(
      "AnnotationDbi"
    )
  ),
  "\n",
  sep = ""
)

cat(
  "org.Hs.eg.db version: ",
  as.character(
    packageVersion(
      "org.Hs.eg.db"
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Os metadados completos da base foram salvos em:\n  ",
  annotation_metadata_output,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 3. Normalização dos Ensembl
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "3. NORMALIZAÇÃO DOS IDENTIFICADORES ENSEMBL\n"
)
cat(
  "------------------------------------------------------------\n"
)

print(
  format_summary
)

cat(
  "\nFeatures com stable Ensembl extraível: ",
  fmt_int(
    n_ensembl_extractable
  ),
  " / ",
  fmt_int(
    n_features_input
  ),
  "\n",
  sep = ""
)

cat(
  "Features sem Ensembl reconhecível: ",
  fmt_int(
    n_non_ensembl
  ),
  "\n",
  sep = ""
)

cat(
  "Stable Ensembl IDs únicos: ",
  fmt_int(
    length(
      ensembl_query
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Stable Ensembl com >1 feature original: ",
  fmt_int(
    n_stable_ids_multiple_features
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 4. Ensembl → GeneID
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "4. MAPEAMENTO Ensembl → GeneID\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Stable Ensembl mapeados de forma unívoca: ",
  fmt_int(
    n_ensembl_unique_mapped
  ),
  "\n",
  sep = ""
)

cat(
  "Stable Ensembl sem GeneID: ",
  fmt_int(
    n_ensembl_unmapped
  ),
  "\n",
  sep = ""
)

cat(
  "Stable Ensembl ambíguos (>1 GeneID): ",
  fmt_int(
    n_ensembl_ambiguous
  ),
  "\n\n",
  sep = ""
)

cat(
  "Política:\n"
)
cat(
  "  1 GeneID  → manter\n"
)
cat(
  "  0 GeneID  → registrar e excluir da matriz harmonizada\n"
)
cat(
  "  >1 GeneID → registrar como ambíguo e excluir, sem duplicar counts\n\n"
)


# ------------------------------------------------------------------------------
# 5. Colapso
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "5. COLAPSO PARA GeneID ÚNICO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Features originais com mapeamento unívoco: ",
  fmt_int(
    nrow(
      mapped_features
    )
  ),
  "\n",
  sep = ""
)

cat(
  "GeneID únicos antes/depois do colapso: ",
  fmt_int(
    uniqueN(
      mapped_features$GeneID
    )
  ),
  "\n",
  sep = ""
)

cat(
  "GeneID recebendo >1 feature original: ",
  fmt_int(
    n_geneids_collapsed
  ),
  "\n",
  sep = ""
)

cat(
  "Features envolvidas nesses colapsos: ",
  fmt_int(
    ifelse(
      is.na(
        n_input_features_in_collapses
      ),
      0,
      n_input_features_in_collapses
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Quando várias features convergiram para o mesmo GeneID, os counts foram\n"
)
cat(
  "somados por amostra. Nenhuma média ou arredondamento foi aplicado.\n\n"
)


# ------------------------------------------------------------------------------
# 6. Balanço
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "6. BALANÇO DE COUNTS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Counts totais na matriz de entrada: ",
  fmt_int(
    total_counts_input
  ),
  "\n",
  sep = ""
)

cat(
  "Counts retidos após seleção de mapeamentos unívocos: ",
  fmt_int(
    total_counts_retained_before_collapse
  ),
  "\n",
  sep = ""
)

cat(
  "Counts após colapso GeneID: ",
  fmt_int(
    total_counts_output
  ),
  "\n",
  sep = ""
)

cat(
  "Counts excluídos por ausência/ambiguidade de mapeamento: ",
  fmt_int(
    counts_excluded_from_harmonized
  ),
  "\n",
  sep = ""
)

cat(
  "Percentual de counts retidos: ",
  fmt_num(
    percent_counts_retained,
    4
  ),
  "%\n",
  sep = ""
)

cat(
  "Percentual de counts excluídos: ",
  fmt_num(
    percent_counts_excluded,
    4
  ),
  "%\n\n",
  sep = ""
)

cat(
  "A soma dos counts das features retidas foi preservada exatamente durante\n"
)
cat(
  "o colapso por GeneID.\n\n"
)


# ------------------------------------------------------------------------------
# 7. Matriz final
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "7. MATRIZ HARMONIZADA FINAL\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "GeneID únicos: ",
  fmt_int(
    nrow(
      counts_harmonized_matrix
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras: ",
  fmt_int(
    ncol(
      counts_harmonized_matrix
    )
  ),
  "\n",
  sep = ""
)

cat(
  "GeneID duplicados: 0\n"
)

cat(
  "NA: ",
  fmt_int(
    n_na_output
  ),
  "\n",
  sep = ""
)

cat(
  "Não finitos: ",
  fmt_int(
    n_nonfinite_output
  ),
  "\n",
  sep = ""
)

cat(
  "Negativos: ",
  fmt_int(
    n_negative_output
  ),
  "\n",
  sep = ""
)

cat(
  "Fracionários: ",
  fmt_int(
    n_fractional_output
  ),
  "\n",
  sep = ""
)

cat(
  "Genes all-zero: ",
  fmt_int(
    n_all_zero_output
  ),
  "\n",
  sep = ""
)

cat(
  "Amostras all-zero: ",
  fmt_int(
    n_all_zero_samples_output
  ),
  "\n",
  sep = ""
)

cat(
  "Ordem counts ↔ metadata preservada: TRUE\n\n"
)


# ------------------------------------------------------------------------------
# 8. Baixa expressão
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "8. FILTRAGEM DE BAIXA EXPRESSÃO — NÃO APLICADA\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Nenhum gene foi removido por baixa expressão nesta etapa.\n"
)

cat(
  "O critério pré-definido para o 04_qc_counts_ROSMAP.R permanece:\n"
)

cat(
  "  count >= 10 em pelo menos ",
  fmt_int(
    EXPECTED_CONTROL_N
  ),
  " amostras.\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 9. Arquivos
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "9. ARQUIVOS GERADOS\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "Matriz GeneID:\n  ",
  counts_output,
  "\n\n",
  sep = ""
)

cat(
  "Mapeamento detalhado:\n  ",
  mapping_detailed_output,
  "\n\n",
  sep = ""
)

cat(
  "Resumo do mapeamento:\n  ",
  mapping_summary_output,
  "\n\n",
  sep = ""
)

cat(
  "Unmapped:\n  ",
  unmapped_output,
  "\n\n",
  sep = ""
)

cat(
  "Ambíguos:\n  ",
  ambiguous_output,
  "\n\n",
  sep = ""
)

cat(
  "Colapsos múltiplas features → GeneID:\n  ",
  collapsed_output,
  "\n\n",
  sep = ""
)

cat(
  "Anotação GeneID:\n  ",
  gene_annotation_output,
  "\n\n",
  sep = ""
)

cat(
  "Balanço de counts:\n  ",
  count_balance_output,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 10. Conclusão
# ------------------------------------------------------------------------------

cat(
  "------------------------------------------------------------\n"
)
cat(
  "10. CONCLUSÃO\n"
)
cat(
  "------------------------------------------------------------\n"
)

cat(
  "- Os identificadores originais foram preservados em tabela de auditoria.\n"
)

cat(
  "- A versão/sufixos dos Ensembl foram removidos apenas para obter o stable ID.\n"
)

cat(
  "- O stable Ensembl foi mapeado para NCBI GeneID usando org.Hs.eg.db.\n"
)

cat(
  "- Mapeamentos ausentes e ambíguos foram documentados e não duplicaram counts.\n"
)

cat(
  "- Várias features convergindo para o mesmo GeneID foram somadas por amostra.\n"
)

cat(
  "- A matriz final possui GeneID único e mantém a ordem das 144 amostras.\n"
)

cat(
  "- Nenhum filtro de baixa expressão foi aplicado nesta etapa.\n"
)

cat(
  "- A matriz está pronta para 04_qc_counts_ROSMAP.R.\n\n"
)

cat(
  "============================================================\n"
)
cat(
  "03_harmonize_genes_ROSMAP.R CONCLUÍDO\n"
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
  "HARMONIZAÇÃO ROSMAP CONCLUÍDA\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Features de entrada: ",
  fmt_int(
    n_features_input
  ),
  "\n",
  sep = ""
)

cat(
  "Stable Ensembl extraíveis: ",
  fmt_int(
    n_ensembl_extractable
  ),
  "\n",
  sep = ""
)

cat(
  "Mapeamentos unívocos Ensembl → GeneID: ",
  fmt_int(
    n_ensembl_unique_mapped
  ),
  "\n",
  sep = ""
)

cat(
  "Ensembl sem GeneID: ",
  fmt_int(
    n_ensembl_unmapped
  ),
  "\n",
  sep = ""
)

cat(
  "Ensembl ambíguos: ",
  fmt_int(
    n_ensembl_ambiguous
  ),
  "\n",
  sep = ""
)

cat(
  "GeneID finais únicos: ",
  fmt_int(
    nrow(
      counts_harmonized_matrix
    )
  ),
  "\n",
  sep = ""
)

cat(
  "GeneID com colapso de >1 feature: ",
  fmt_int(
    n_geneids_collapsed
  ),
  "\n\n",
  sep = ""
)

cat(
  "Counts retidos: ",
  fmt_num(
    percent_counts_retained,
    4
  ),
  "%\n",
  sep = ""
)

cat(
  "Genes all-zero após harmonização: ",
  fmt_int(
    n_all_zero_output
  ),
  "\n\n",
  sep = ""
)

cat(
  "Matriz harmonizada:\n  ",
  counts_output,
  "\n\n",
  sep = ""
)

cat(
  "Relatório:\n  ",
  report_output,
  "\n\n",
  sep = ""
)

cat(
  "Próxima etapa, após revisar este relatório:\n"
)
cat(
  "  04_qc_counts_ROSMAP.R\n"
)
