#!/usr/bin/env Rscript

# ==============================================================================
# 01_prepare_metadata_GSE125583.R
#
# FASE 01b — PREPARAÇÃO DO METADATA DO DATASET CANDIDATO DE AD (GSE125583)
#
# OBJETIVO GERAL
# ------------------------------------------------------------------------------
# Construir uma tabela de metadados padronizada, uma linha por amostra, a partir
# do arquivo SOFT do GEO e validá-la contra as colunas da matriz de raw counts
# padronizada pelo NCBI.
#
# A etapa 00 já mostrou que:
#   - a matriz possui 39.376 genes e 289 amostras;
#   - counts e SOFT possuem as mesmas amostras;
#   - existem 219 amostras de Alzheimer e 70 controles;
#   - há 1 sexo ausente, 1 idade ausente e 2 Braak scores ausentes;
#   - todos os tecidos são "fusiform gyrus";
#   - não há case IDs ausentes ou duplicados.
#
# Esta etapa transforma essas informações em um metadata canônico que será
# utilizado nas etapas seguintes do pipeline.
#
#
# O QUE ESTE SCRIPT FAZ
# ------------------------------------------------------------------------------
# 1. Lê o arquivo GSE125583_family.soft.gz.
# 2. Extrai, para cada amostra:
#      - sample_id
#      - title
#      - case_id
#      - diagnosis
#      - tissue
#      - age
#      - sex
#      - braak_score
# 3. Padroniza o grupo experimental:
#      "Alzheimer's disease" -> "AD"
#      "control"             -> "Control"
# 4. Converte idade para valor numérico.
# 5. Mantém a informação original do diagnóstico em diagnosis_original.
# 6. Valida unicidade de sample_id e case_id.
# 7. Valida correspondência 1:1 com as amostras da matriz de counts.
# 8. Reordena o metadata para a MESMA ordem das colunas da matriz.
# 9. Audita:
#      - tamanho dos grupos;
#      - group x sex;
#      - group x Braak;
#      - idade por grupo;
#      - amostras com covariáveis ausentes.
# 10. Salva o metadata final e um relatório textual de auditoria.
#
#
# O QUE ESTE SCRIPT NÃO FAZ
# ------------------------------------------------------------------------------
# - Não exclui nenhuma amostra.
# - Não decide quais covariáveis entrarão no DESeq2.
# - Não filtra genes.
# - Não modifica a matriz de counts.
# - Não executa PCA.
# - Não executa DESeq2.
# - Não usa Braak para redefinir diagnóstico.
#
# O diagnóstico clínico informado pelo estudo é preservado como definição de
# grupo nesta etapa. Braak, idade e sexo são mantidos como covariáveis para
# avaliação posterior.
#
#
# ENTRADAS
# ------------------------------------------------------------------------------
# data/raw/GSE125583/GSE125583_family.soft.gz
# data/raw/GSE125583/GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz
#
#
# SAÍDAS
# ------------------------------------------------------------------------------
# results/candidate_GSE125583/metadata/metadata_GSE125583.tsv
# results/candidate_GSE125583/metadata/01_prepare_metadata_GSE125583.txt
# ==============================================================================


# ==============================================================================
# 1. PACOTES
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
# 2. CAMINHOS
# ==============================================================================

soft_file <- file.path(
  "data",
  "raw",
  "GSE125583",
  "GSE125583_family.soft.gz"
)

counts_file <- file.path(
  "data",
  "raw",
  "GSE125583",
  "GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)

output_dir <- file.path(
  "results",
  "candidate_GSE125583",
  "metadata"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

metadata_file <- file.path(
  output_dir,
  "metadata_GSE125583.tsv"
)

audit_file <- file.path(
  output_dir,
  "01_prepare_metadata_GSE125583.txt"
)


# ==============================================================================
# 3. CONSTANTES ESPERADAS
# ==============================================================================
#
# Estes valores não são usados para "forçar" o resultado.
# Eles funcionam como verificações de consistência baseadas na auditoria da
# etapa 00. Se alguma dessas condições mudar inesperadamente, o script para
# para evitar que o pipeline siga silenciosamente com uma coorte diferente.

EXPECTED_N_SAMPLES <- 289L
EXPECTED_N_AD <- 219L
EXPECTED_N_CONTROL <- 70L

EXPECTED_TISSUE <- "fusiform gyrus"

VALID_DIAGNOSES <- c(
  "Alzheimer's disease",
  "control"
)

VALID_GROUPS <- c(
  "Control",
  "AD"
)

VALID_SEX <- c(
  "F",
  "M"
)

VALID_BRAAK <- c(
  "I",
  "II",
  "III",
  "IV",
  "V",
  "VI"
)


# ==============================================================================
# 4. FUNÇÕES AUXILIARES
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 Validar arquivo
# ------------------------------------------------------------------------------
# Interrompe o script se o arquivo não existir ou estiver vazio.

validate_file <- function(path) {

  if (!file.exists(path)) {
    stop(
      "Arquivo não encontrado: ",
      path
    )
  }

  info <- file.info(path)

  if (is.na(info$size) || info$size == 0) {
    stop(
      "Arquivo vazio ou tamanho indisponível: ",
      path
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 4.2 Normalizar valores ausentes
# ------------------------------------------------------------------------------
# O SOFT pode representar ausência como "NA", "N/A", string vazia etc.
# Todos esses casos são convertidos para NA real do R.

normalize_missing_character <- function(x) {

  x <- trimws(
    as.character(x)
  )

  x[
    x %in% c(
      "",
      "NA",
      "N/A",
      "na",
      "n/a"
    )
  ] <- NA_character_

  x
}


# ------------------------------------------------------------------------------
# 4.3 Formatar inteiros para relatório
# ------------------------------------------------------------------------------
# Exemplo:
#   39376 -> "39.376"

format_int <- function(x) {

  format(
    x,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


# ------------------------------------------------------------------------------
# 4.4 Formatar números para relatório
# ------------------------------------------------------------------------------

format_num <- function(
  x,
  digits = 2
) {

  format(
    round(
      x,
      digits
    ),
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


# ------------------------------------------------------------------------------
# 4.5 Extrair um campo simples de um bloco SAMPLE
# ------------------------------------------------------------------------------
# Exemplo:
#   !Sample_title = AD 04-45

extract_single_field <- function(
  block,
  prefix
) {

  hit <- block[
    startsWith(
      block,
      prefix
    )
  ]

  if (length(hit) == 0) {
    return(
      NA_character_
    )
  }

  trimws(
    sub(
      prefix,
      "",
      hit[1],
      fixed = TRUE
    )
  )
}


# ------------------------------------------------------------------------------
# 4.6 Extrair uma característica do SOFT
# ------------------------------------------------------------------------------
# Exemplo:
#   !Sample_characteristics_ch1 = diagnosis: Alzheimer's disease
#
# A função recebe apenas a chave, como "diagnosis", e devolve o valor.

extract_characteristic <- function(
  block,
  key
) {

  prefix <- paste0(
    "!Sample_characteristics_ch1 = ",
    key,
    ":"
  )

  hit <- block[
    startsWith(
      block,
      prefix
    )
  ]

  if (length(hit) == 0) {
    return(
      NA_character_
    )
  }

  trimws(
    sub(
      prefix,
      "",
      hit[1],
      fixed = TRUE
    )
  )
}


# ------------------------------------------------------------------------------
# 4.7 Ler metadata diretamente do arquivo SOFT
# ------------------------------------------------------------------------------
# O arquivo SOFT é dividido em blocos iniciados por:
#   ^SAMPLE = GSM...
#
# Para cada bloco são extraídas apenas as variáveis relevantes para esta fase.

parse_soft_metadata <- function(path) {

  con <- gzfile(
    path,
    open = "rt"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  lines <- readLines(
    con,
    warn = FALSE
  )

  sample_start <- grep(
    "^\\^SAMPLE = ",
    lines
  )

  if (length(sample_start) == 0) {
    stop(
      "Nenhuma seção ^SAMPLE encontrada no arquivo SOFT."
    )
  }

  sample_end <- c(
    sample_start[-1] - 1L,
    length(lines)
  )

  out <- vector(
    "list",
    length(sample_start)
  )

  for (i in seq_along(sample_start)) {

    block <- lines[
      sample_start[i]:sample_end[i]
    ]

    sample_id <- sub(
      "^\\^SAMPLE = ",
      "",
      block[1]
    )

    out[[i]] <- data.table(
      sample_id = trimws(
        sample_id
      ),

      title = extract_single_field(
        block,
        "!Sample_title = "
      ),

      case_id = extract_characteristic(
        block,
        "case id"
      ),

      diagnosis_original = extract_characteristic(
        block,
        "diagnosis"
      ),

      tissue = extract_characteristic(
        block,
        "tissue"
      ),

      age_raw = extract_characteristic(
        block,
        "age"
      ),

      sex = extract_characteristic(
        block,
        "Sex"
      ),

      braak_score = extract_characteristic(
        block,
        "braak.score"
      )
    )
  }

  metadata <- rbindlist(
    out,
    fill = TRUE
  )

  # Normalizar representações textuais de ausência.
  character_columns <- c(
    "sample_id",
    "title",
    "case_id",
    "diagnosis_original",
    "tissue",
    "age_raw",
    "sex",
    "braak_score"
  )

  for (col_name in character_columns) {

    set(
      metadata,
      j = col_name,
      value = normalize_missing_character(
        metadata[[col_name]]
      )
    )
  }

  metadata
}


# ------------------------------------------------------------------------------
# 4.8 Ler apenas os IDs das amostras da matriz de counts
# ------------------------------------------------------------------------------
# Não precisamos carregar toda a matriz em memória nesta etapa.
#
# fread(nrows = 0) lê somente o cabeçalho.

read_count_sample_ids <- function(path) {

  header <- fread(
    path,
    nrows = 0,
    showProgress = FALSE
  )

  if (!"GeneID" %in% names(header)) {
    stop(
      "A matriz de counts não contém a coluna GeneID."
    )
  }

  if (names(header)[1] != "GeneID") {
    stop(
      "A primeira coluna da matriz de counts deveria ser GeneID."
    )
  }

  sample_ids <- setdiff(
    names(header),
    "GeneID"
  )

  if (length(sample_ids) == 0) {
    stop(
      "Nenhuma amostra encontrada na matriz de counts."
    )
  }

  if (anyDuplicated(sample_ids)) {
    stop(
      "Existem sample IDs duplicados no cabeçalho da matriz de counts."
    )
  }

  sample_ids
}


# ==============================================================================
# 5. VALIDAR ENTRADAS
# ==============================================================================

validate_file(
  soft_file
)

validate_file(
  counts_file
)


# ==============================================================================
# 6. LER O SOFT E CONSTRUIR METADATA INICIAL
# ==============================================================================

metadata <- parse_soft_metadata(
  soft_file
)

if (nrow(metadata) != EXPECTED_N_SAMPLES) {
  stop(
    "Número inesperado de amostras no SOFT. Esperado: ",
    EXPECTED_N_SAMPLES,
    "; observado: ",
    nrow(metadata),
    "."
  )
}


# ==============================================================================
# 7. VALIDAR IDENTIFICADORES
# ==============================================================================

# ------------------------------------------------------------------------------
# 7.1 sample_id
# ------------------------------------------------------------------------------

if (any(is.na(metadata$sample_id))) {
  stop(
    "Existem sample IDs ausentes no SOFT."
  )
}

if (anyDuplicated(metadata$sample_id)) {
  stop(
    "Existem sample IDs duplicados no SOFT."
  )
}


# ------------------------------------------------------------------------------
# 7.2 case_id
# ------------------------------------------------------------------------------
# A etapa 00 encontrou zero ausentes e zero duplicados.
# Validamos novamente porque o case_id representa o indivíduo.

if (any(is.na(metadata$case_id))) {
  stop(
    "Existem case IDs ausentes."
  )
}

if (anyDuplicated(metadata$case_id)) {
  stop(
    "Existem case IDs duplicados. ",
    "Revisar possível existência de múltiplas amostras por indivíduo."
  )
}


# ==============================================================================
# 8. PADRONIZAR O GRUPO EXPERIMENTAL
# ==============================================================================

unexpected_diagnoses <- setdiff(
  unique(
    na.omit(
      metadata$diagnosis_original
    )
  ),
  VALID_DIAGNOSES
)

if (length(unexpected_diagnoses) > 0) {
  stop(
    "Diagnósticos inesperados encontrados: ",
    paste(
      unexpected_diagnoses,
      collapse = ", "
    )
  )
}

if (any(is.na(metadata$diagnosis_original))) {
  stop(
    "Existem diagnósticos ausentes."
  )
}

metadata[
  diagnosis_original == "Alzheimer's disease",
  group := "AD"
]

metadata[
  diagnosis_original == "control",
  group := "Control"
]

if (any(is.na(metadata$group))) {
  stop(
    "Falha ao padronizar o grupo experimental."
  )
}

if (!all(metadata$group %in% VALID_GROUPS)) {
  stop(
    "Foram gerados grupos inesperados."
  )
}


# ==============================================================================
# 9. CONVERTER E VALIDAR IDADE
# ==============================================================================

# Guardamos age_raw apenas durante a preparação.
# O metadata final conterá "age" numérica.

metadata[
  ,
  age := suppressWarnings(
    as.numeric(age_raw)
  )
]

# Se o valor original não era NA mas falhou na conversão numérica, isso indica
# uma estrutura inesperada do SOFT e deve ser investigado.

invalid_age_conversion <- (
  !is.na(metadata$age_raw) &
    is.na(metadata$age)
)

if (any(invalid_age_conversion)) {

  stop(
    "Existem valores de idade não ausentes que não puderam ser convertidos ",
    "para número: ",
    paste(
      unique(
        metadata$age_raw[
          invalid_age_conversion
        ]
      ),
      collapse = ", "
    )
  )
}

# Checagem ampla de plausibilidade, sem definir nenhum critério de exclusão.
if (
  any(
    metadata$age < 0 |
      metadata$age > 120,
    na.rm = TRUE
  )
) {
  stop(
    "Foram encontradas idades fora do intervalo plausível de 0 a 120 anos."
  )
}


# ==============================================================================
# 10. VALIDAR SEXO, BRAAK E TECIDO
# ==============================================================================

# ------------------------------------------------------------------------------
# 10.1 Sexo
# ------------------------------------------------------------------------------
# Um valor ausente já era esperado.
# Valores válidos observados no estudo: F e M.

unexpected_sex <- setdiff(
  unique(
    na.omit(
      metadata$sex
    )
  ),
  VALID_SEX
)

if (length(unexpected_sex) > 0) {
  stop(
    "Valores inesperados de sexo encontrados: ",
    paste(
      unexpected_sex,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# 10.2 Braak
# ------------------------------------------------------------------------------
# O score é preservado como categoria ordinal textual.
# Nesta etapa NÃO usamos Braak para redefinir AD/Control.

unexpected_braak <- setdiff(
  unique(
    na.omit(
      metadata$braak_score
    )
  ),
  VALID_BRAAK
)

if (length(unexpected_braak) > 0) {
  stop(
    "Valores inesperados de Braak score encontrados: ",
    paste(
      unexpected_braak,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# 10.3 Tecido
# ------------------------------------------------------------------------------

if (any(is.na(metadata$tissue))) {
  stop(
    "Existem valores de tecido ausentes."
  )
}

unexpected_tissue <- setdiff(
  unique(
    metadata$tissue
  ),
  EXPECTED_TISSUE
)

if (length(unexpected_tissue) > 0) {
  stop(
    "Foram encontrados tecidos inesperados: ",
    paste(
      unexpected_tissue,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 11. VALIDAR TAMANHO DOS GRUPOS
# ==============================================================================

n_ad <- metadata[
  group == "AD",
  .N
]

n_control <- metadata[
  group == "Control",
  .N
]

if (n_ad != EXPECTED_N_AD) {
  stop(
    "Número inesperado de amostras AD. Esperado: ",
    EXPECTED_N_AD,
    "; observado: ",
    n_ad,
    "."
  )
}

if (n_control != EXPECTED_N_CONTROL) {
  stop(
    "Número inesperado de controles. Esperado: ",
    EXPECTED_N_CONTROL,
    "; observado: ",
    n_control,
    "."
  )
}


# ==============================================================================
# 12. VALIDAR COUNTS <-> METADATA
# ==============================================================================

count_sample_ids <- read_count_sample_ids(
  counts_file
)

metadata_sample_ids <- metadata$sample_id

samples_only_counts <- setdiff(
  count_sample_ids,
  metadata_sample_ids
)

samples_only_metadata <- setdiff(
  metadata_sample_ids,
  count_sample_ids
)

if (length(samples_only_counts) > 0) {
  stop(
    "Existem amostras presentes nos counts e ausentes no metadata: ",
    paste(
      samples_only_counts,
      collapse = ", "
    )
  )
}

if (length(samples_only_metadata) > 0) {
  stop(
    "Existem amostras presentes no metadata e ausentes nos counts: ",
    paste(
      samples_only_metadata,
      collapse = ", "
    )
  )
}

if (length(count_sample_ids) != nrow(metadata)) {
  stop(
    "Número de amostras diferente entre counts e metadata."
  )
}

# ------------------------------------------------------------------------------
# Reordenar metadata segundo a matriz de counts
# ------------------------------------------------------------------------------
# Mesmo que a etapa 00 já tenha mostrado que a ordem coincide, fazemos o
# reordenamento explicitamente para tornar o script robusto.
#
# Isso garante que:
#   metadata$sample_id
# seja exatamente igual à ordem das colunas da matriz.

metadata <- metadata[
  match(
    count_sample_ids,
    sample_id
  )
]

if (!identical(
  metadata$sample_id,
  count_sample_ids
)) {
  stop(
    "Falha ao alinhar a ordem do metadata com a matriz de counts."
  )
}


# ==============================================================================
# 13. ORGANIZAR O METADATA FINAL
# ==============================================================================

# age_raw é removida porque foi utilizada apenas para validar a conversão.
#
# A ordem das colunas foi escolhida para deixar:
#   1. identificadores;
#   2. definição do grupo;
#   3. informação original do GEO;
#   4. características biológicas/demográficas.

metadata_final <- metadata[
  ,
  .(
    sample_id,
    title,
    case_id,
    group,
    diagnosis_original,
    tissue,
    age,
    sex,
    braak_score
  )
]


# ==============================================================================
# 14. AUDITORIAS DESCRITIVAS
# ==============================================================================

# ------------------------------------------------------------------------------
# 14.1 Grupo
# ------------------------------------------------------------------------------

group_counts <- metadata_final[
  ,
  .N,
  by = group
][
  match(
    VALID_GROUPS,
    group
  )
]


# ------------------------------------------------------------------------------
# 14.2 Grupo x sexo
# ------------------------------------------------------------------------------

group_sex <- metadata_final[
  ,
  .N,
  by = .(
    group,
    sex
  )
][
  order(
    factor(
      group,
      levels = VALID_GROUPS
    ),
    sex,
    na.last = TRUE
  )
]


# ------------------------------------------------------------------------------
# 14.3 Grupo x Braak
# ------------------------------------------------------------------------------

group_braak <- metadata_final[
  ,
  .N,
  by = .(
    group,
    braak_score
  )
][
  order(
    factor(
      group,
      levels = VALID_GROUPS
    ),
    factor(
      braak_score,
      levels = VALID_BRAAK
    ),
    na.last = TRUE
  )
]


# ------------------------------------------------------------------------------
# 14.4 Idade por grupo
# ------------------------------------------------------------------------------

age_by_group <- metadata_final[
  ,
  .(
    n = .N,
    n_age_missing = sum(
      is.na(age)
    ),
    age_min = if (
      all(
        is.na(age)
      )
    ) {
      NA_real_
    } else {
      min(
        age,
        na.rm = TRUE
      )
    },
    age_median = if (
      all(
        is.na(age)
      )
    ) {
      NA_real_
    } else {
      median(
        age,
        na.rm = TRUE
      )
    },
    age_mean = if (
      all(
        is.na(age)
      )
    ) {
      NA_real_
    } else {
      mean(
        age,
        na.rm = TRUE
      )
    },
    age_sd = if (
      sum(
        !is.na(age)
      ) <= 1
    ) {
      NA_real_
    } else {
      sd(
        age,
        na.rm = TRUE
      )
    },
    age_max = if (
      all(
        is.na(age)
      )
    ) {
      NA_real_
    } else {
      max(
        age,
        na.rm = TRUE
      )
    }
  ),
  by = group
][
  match(
    VALID_GROUPS,
    group
  )
]


# ------------------------------------------------------------------------------
# 14.5 Amostras com informações ausentes
# ------------------------------------------------------------------------------

missing_covariates <- metadata_final[
  is.na(age) |
    is.na(sex) |
    is.na(braak_score)
]


# ==============================================================================
# 15. SALVAR METADATA FINAL
# ==============================================================================

fwrite(
  metadata_final,
  metadata_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 16. GERAR RELATÓRIO TEXTUAL
# ==============================================================================

sink(
  audit_file
)

cat("============================================================\n")
cat("01 — PREPARAÇÃO DO METADATA GSE125583\n")
cat("============================================================\n\n")


cat("------------------------------------------------------------\n")
cat("1. METADATA FINAL\n")
cat("------------------------------------------------------------\n")

cat(
  "Amostras: ",
  format_int(
    nrow(
      metadata_final
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Colunas: ",
  format_int(
    ncol(
      metadata_final
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Sample IDs únicos: ",
  !anyDuplicated(
    metadata_final$sample_id
  ),
  "\n",
  sep = ""
)

cat(
  "Case IDs únicos: ",
  !anyDuplicated(
    metadata_final$case_id
  ),
  "\n",
  sep = ""
)

cat(
  "Counts e metadata alinhados na mesma ordem: ",
  identical(
    metadata_final$sample_id,
    count_sample_ids
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("2. DISTRIBUIÇÃO DOS GRUPOS\n")
cat("------------------------------------------------------------\n")

print(
  group_counts
)

cat("\n")


cat("------------------------------------------------------------\n")
cat("3. GRUPO x SEXO\n")
cat("------------------------------------------------------------\n")

print(
  group_sex
)

cat(
  "\nSexo ausente: ",
  format_int(
    sum(
      is.na(
        metadata_final$sex
      )
    )
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("4. GRUPO x BRAAK SCORE\n")
cat("------------------------------------------------------------\n")

print(
  group_braak
)

cat(
  "\nBraak score ausente: ",
  format_int(
    sum(
      is.na(
        metadata_final$braak_score
      )
    )
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("5. IDADE POR GRUPO\n")
cat("------------------------------------------------------------\n")

print(
  age_by_group
)

cat(
  "\nIdade ausente: ",
  format_int(
    sum(
      is.na(
        metadata_final$age
      )
    )
  ),
  "\n\n",
  sep = ""
)


cat("------------------------------------------------------------\n")
cat("6. TECIDO\n")
cat("------------------------------------------------------------\n")

print(
  metadata_final[
    ,
    .N,
    by = tissue
  ]
)

cat("\n")


cat("------------------------------------------------------------\n")
cat("7. AMOSTRAS COM COVARIÁVEIS AUSENTES\n")
cat("------------------------------------------------------------\n")

if (nrow(missing_covariates) == 0) {

  cat(
    "Nenhuma amostra possui age, sex ou braak_score ausentes.\n"
  )

} else {

  print(
    missing_covariates[
      ,
      .(
        sample_id,
        title,
        case_id,
        group,
        age,
        sex,
        braak_score
      )
    ]
  )
}

cat("\n")


cat("------------------------------------------------------------\n")
cat("8. INTERPRETAÇÃO DA ETAPA\n")
cat("------------------------------------------------------------\n")

cat(
  "- O metadata foi construído diretamente do SOFT do GSE125583.\n"
)

cat(
  "- O diagnóstico original do GEO foi preservado em diagnosis_original.\n"
)

cat(
  "- Os grupos foram padronizados para AD e Control.\n"
)

cat(
  "- Nenhuma amostra foi removida nesta etapa.\n"
)

cat(
  "- Valores ausentes de idade, sexo ou Braak foram preservados como NA.\n"
)

cat(
  "- Braak score não foi utilizado para redefinir diagnóstico.\n"
)

cat(
  "- O metadata foi reordenado para corresponder exatamente à ordem das colunas da matriz de counts.\n"
)

cat(
  "- A decisão sobre covariáveis do modelo DESeq2 será tomada posteriormente, após QC e avaliação das distribuições.\n"
)

cat("\n")
cat("============================================================\n")
cat("01_prepare_metadata_GSE125583.R CONCLUÍDO\n")
cat("============================================================\n")

sink()


# ==============================================================================
# 17. RESUMO NO TERMINAL
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("01 — METADATA GSE125583 CONCLUÍDO\n")
cat("============================================================\n")

cat(
  "Amostras totais: ",
  format_int(
    nrow(
      metadata_final
    )
  ),
  "\n",
  sep = ""
)

cat(
  "AD: ",
  format_int(
    n_ad
  ),
  "\n",
  sep = ""
)

cat(
  "Control: ",
  format_int(
    n_control
  ),
  "\n",
  sep = ""
)

cat(
  "Idade ausente: ",
  format_int(
    sum(
      is.na(
        metadata_final$age
      )
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Sexo ausente: ",
  format_int(
    sum(
      is.na(
        metadata_final$sex
      )
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Braak ausente: ",
  format_int(
    sum(
      is.na(
        metadata_final$braak_score
      )
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Counts <-> metadata alinhados: ",
  identical(
    metadata_final$sample_id,
    count_sample_ids
  ),
  "\n",
  sep = ""
)

cat(
  "Metadata: ",
  metadata_file,
  "\n",
  sep = ""
)

cat(
  "Relatório: ",
  audit_file,
  "\n",
  sep = ""
)

cat("============================================================\n")
