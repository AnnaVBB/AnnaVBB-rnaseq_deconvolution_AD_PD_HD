#!/usr/bin/env Rscript

# ==============================================================================
# 00_audit_raw_data.R
# Auditoria inicial e reprodutivel dos dados de bulk RNA-seq.
#
# Datasets:
#   GSE53697 -> Alzheimer (AD)
#   GSE64810 -> Huntington (HD)
#   GSE68719 -> Parkinson (PD)
#
# Objetivos:
# 1. Verificar integridade e estrutura dos arquivos de entrada.
# 2. Auditar raw counts sem modificar os dados.
# 3. Auditar GeneIDs e a anotacao GRCh38.p13.
# 4. Quantificar a sobreposicao genica entre AD, HD e PD.
# 5. Registrar genes all-zero e valores fracionarios.
# 6. Produzir uma tabela por GeneID para uso posterior na harmonizacao.
#
# Outputs:
#   results/audit/00_audit_raw_data.txt
#   results/audit/00_gene_universe_audit.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

file_ad <- "data/raw/GSE53697/GSE53697_RNAseq_AD.txt.gz"

file_hd <- paste0(
  "data/raw/GSE64810/",
  "GSE64810_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)

file_pd <- paste0(
  "data/raw/GSE68719/",
  "GSE68719_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)

file_annotation <- "data/raw/Human.GRCh38.p13.annot.tsv.gz"

output_dir <- "results/audit"

output_report <- file.path(
  output_dir,
  "00_audit_raw_data.txt"
)

output_gene_universe <- file.path(
  output_dir,
  "00_gene_universe_audit.tsv"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

#Organiza a saída no terminal, deixando o log mais organizado ao criar divisórias usando o sinal = repetidas vezes entre o título (passando-o como argumento)
section <- function(title) {
  cat("\n")
  cat(paste(rep("=", 78), collapse = ""), "\n")
  cat(title, "\n")
  cat(paste(rep("=", 78), collapse = ""), "\n")
}

#Verifica se um arquivo é válido antes de usá-lo
#path é o caminho físico do arquivo, label é o nome amigável para aparecer no log
check_file <- function(path, label) {

  if (!file.exists(path)) { #se o arquivo não existir executamos, interrompemos o script e uma mensagem de erro é gerada
    stop(
      "Arquivo ausente [", label, "]: ",
      path,
      call. = FALSE #deixa a mensagem de erro mais limpa
    )
  }

  info <- file.info(path) #recupera os metadados do arquivo

  if (is.na(info$size) || info$size == 0) { #se o tamanho não puder ser determinado OU o arquivo tiver zero bytes → interrompa 
    stop(
      "Arquivo vazio [", label, "]: ",
      path,
      call. = FALSE
    )
  }

  #Se tudo estiver correto são criadas strings formatadas
  cat(
    sprintf(
      "%-12s OK | %.2f MB | %s\n", #As informações serão inseridas em %. s=tring que reserva 12 caracteres - alinfados à esquerda
      label,
      info$size / 1024^2,
      path
    )
  )
}

#Limpa identificadores gênicos 
normalize_gene_id <- function(x) {
  x <- trimws(as.character(x)) #os IDs são tratados como texto e o trimws remove os espaços no começo e no fim 
  #Nos elementos de x em que x== "" é TRUE coloque NA
  x[
    is.na(x) |
    x == "" |
    toupper(x) == "NA" #transforma o texto em maiúsculas
  ] <- NA_character_
  #NA real ou string vazia ou texto "NA"
  x
}

clean_annotation_field <- function(x) {
  x <- trimws(as.character(x))
  x[
    is.na(x) |
    x == "" |
    toupper(x) == "NA" |
    x == "-" #semelhante à função anterior, mas transforma os campos com - em NA
  ] <- NA_character_
  x
}

#Calcular mínimo e máximo sem quebrar, essa função responde qual é o menor e o maior valor válido
safe_range <- function(x) {
  x <- x[is.finite(x)] #retorna true se o valor for finito
  if (length(x) == 0L) {
    return(c(NA_real_, NA_real_))
  }
  range(x)
}

#Função de auditoria 
#Recebe a matriz de contagens e o nome do dataset, o objetivo é produzir uma linha-resumo sobre a matriz inteira
audit_counts <- function(count_matrix, dataset) {
  values <- as.vector(count_matrix)
  rg <- safe_range(values)
  #criamos uma tabela resumo
  data.table(
    dataset = dataset,
    n_genes = nrow(count_matrix), #numero de linhas= número de genes
    n_samples = ncol(count_matrix), #número de colunas= número de amostras
    n_values = length(values), #numero de gexes * numero de amostras 
    n_na = sum(is.na(values)), #verificando a quantidade de valores NA
    n_negative = sum(values < 0, na.rm = TRUE), #procurando counts negativos
    n_non_integer = sum( 
      abs(values - round(values)) > 1e-8,
      na.rm = TRUE
    ),  #detecta valores não inteiros
    n_zero = sum(values == 0, na.rm = TRUE), #conta a quantidade de zeros
    pct_zero = 100 *
      sum(values == 0, na.rm = TRUE) /
      sum(!is.na(values)),
    min_count = rg[1], #mínimo 
    max_count = rg[2], #máximo 
    total_counts = sum(values, na.rm = TRUE) #soma toda a matriz 
  )
}

#Função para investigar os valores fracionários. 
audit_fractional <- function(count_matrix, dataset) {
  values <- as.vector(count_matrix)
  fractional <- values[
    !is.na(values) &
    abs(values - round(values)) > 1e-8
  ]

  if (length(fractional) == 0L) {
    return(
      data.table(
        dataset = dataset,
        n_fractional = 0L,
        min_fractional = NA_real_,
        median_fractional = NA_real_,
        max_fractional = NA_real_,
        max_distance_to_integer = 0
      )
    )
  }

  data.table(
    dataset = dataset,
    n_fractional = length(fractional),
    min_fractional = min(fractional),
    median_fractional = median(fractional),
    max_fractional = max(fractional),
    max_distance_to_integer =
      max(abs(fractional - round(fractional)))
  )
}

# Verificamos os arquivos e carregamo-os
section("1. VERIFICACAO DOS ARQUIVOS DE ENTRADA")

check_file(file_ad, "AD")
check_file(file_hd, "HD")
check_file(file_pd, "PD")
check_file(file_annotation, "ANNOTATION")

section("2. CARREGANDO RAW DATA")

cat("Lendo AD...\n")
ad <- fread(file_ad)

cat("Lendo HD...\n")
hd <- fread(file_hd)

cat("Lendo PD...\n")
pd <- fread(file_pd)

cat("Lendo anotacao GRCh38.p13...\n")
annotation <- fread(file_annotation)

cat("Arquivos carregados com sucesso.\n")

# Validação da estrutura
section("3. VALIDACAO DA ESTRUTURA")

#verificação das colunas obrigatórias de AD
#names(ad) retorna os nomes das colunas, verificamos se "GeneID" e "GeneSymbol" estão ali, caso não estejam a auditoria é interrompida e imprime uma mensagem de erro
if (!all(c("GeneID", "GeneSymbol") %in% names(ad))) { 
  stop(
    "AD nao possui GeneID e/ou GeneSymbol.",
    call. = FALSE
  )
}
#Para HD e PD verificamos se há o GeneID
if (!"GeneID" %in% names(hd)) {
  stop(
    "HD nao possui GeneID.",
    call. = FALSE
  )
}

if (!"GeneID" %in% names(pd)) {
  stop(
    "PD nao possui GeneID.",
    call. = FALSE
  )
}

#Como AD possui raw_counts e rpkm precisamos separar as duas coisas. Criamos dois vetores e usamos o grep para procurar um padrão:
#a) encontra textos que terminam em _raw (_raw$)
#b) encontra textos que terminam em _rpkm (_rpkm$)
#O grep retorna posições, mas como queremos os nomes das colunas, value=TRUE
ad_raw_cols <- grep(
  "_raw$",
  names(ad),
  value = TRUE
)

ad_rpkm_cols <- grep(
  "_rpkm$",
  names(ad),
  value = TRUE
)

#HD e PD possuem outra estrutura, todas as colunas são amostras, menos a GeneID
#setdiff pega todos os nomes das colunas do dataset menos a "GeneID"
hd_sample_cols <- setdiff(
  names(hd),
  "GeneID"
)

pd_sample_cols <- setdiff(
  names(pd),
  "GeneID"
)

#Depois de descobrir quais colunas são amostras, verificamos se encontramos a quantidade esperada
#O dataset de AD possui um número de colunas diferente de 17? Caso sim, uma mensagem de erro aparece. O mesmo é feito para HD e PD, verificando se possuem 69 e 73 amostras, respectivamente.
if (length(ad_raw_cols) != 17L) {
  stop(
    "AD: numero inesperado de colunas raw: ",
    length(ad_raw_cols),
    ". Esperado: 17.",
    call. = FALSE
  )
}

if (length(ad_rpkm_cols) != 17L) {
  stop(
    "AD: numero inesperado de colunas RPKM: ",
    length(ad_rpkm_cols),
    ". Esperado: 17.",
    call. = FALSE
  )
}

if (length(hd_sample_cols) != 69L) {
  stop(
    "HD: numero inesperado de amostras: ",
    length(hd_sample_cols),
    ". Esperado: 69.",
    call. = FALSE
  )
}

if (length(pd_sample_cols) != 73L) {
  stop(
    "PD: numero inesperado de amostras: ",
    length(pd_sample_cols),
    ". Esperado: 73.",
    call. = FALSE
  )
}

#Criar uma tabela resumo com três linhas (AD, PD e HD) e combina as informações pela posição
file_dimensions <- data.table(
  dataset = c("AD", "HD", "PD"),
  #calcula o numero de linhas de cada tabela
  n_genes = c(
    nrow(ad),
    nrow(hd),
    nrow(pd)
  ),
  #calcula o numero total de colunas do arquivo
  n_columns_total = c(
    ncol(ad),
    ncol(hd),
    ncol(pd)
  ),
  #Contamos as colunas de raw counts
  n_count_samples = c(
    length(ad_raw_cols),
    length(hd_sample_cols),
    length(pd_sample_cols)
  )
)

print(file_dimensions)

cat("\nAD RPKM columns:", length(ad_rpkm_cols), "\n")


# ------------------------------------------------------------------------------
# Construir matrizes temporarias de raw counts
# ------------------------------------------------------------------------------
#ad[, ..ad_raw_cols] seleciona as 17 colunas de raw counts anteriormente identificadas e guardadas em ad_raw_cols
ad_counts <- as.matrix(
  ad[, ..ad_raw_cols]
)

hd_counts <- as.matrix(
  hd[, ..hd_sample_cols]
)

pd_counts <- as.matrix(
  pd[, ..pd_sample_cols]
)
#As matrizes representam números, portanto devem ser tratadas como numeric
storage.mode(ad_counts) <- "numeric"
storage.mode(hd_counts) <- "numeric"
storage.mode(pd_counts) <- "numeric"


# ------------------------------------------------------------------------------
# Auditoria dos counts
# ------------------------------------------------------------------------------
section("4. AUDITORIA DOS RAW COUNTS")

count_summary <- rbindlist(
  list(
    audit_counts(ad_counts, "AD"),
    audit_counts(hd_counts, "HD"),
    audit_counts(pd_counts, "PD")
  )
)

print(count_summary)

#Interromper a execução caso hajam valores NA
if (any(count_summary$n_na > 0)) {
  stop(
    "Foram detectados valores NA nos raw counts.",
    call. = FALSE
  )
}
#Interromper a execução caso hajam valores negativos
if (any(count_summary$n_negative > 0)) {
  stop(
    "Foram detectados valores negativos nos raw counts.",
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# Valores fracionarios
# ------------------------------------------------------------------------------
section("5. VALORES FRACIONARIOS")
#Como o conjunto de dados de AD possuem valores fracionários, vamos verificá-los
fractional_summary <- rbindlist(
  list(
    audit_fractional(ad_counts, "AD"),
    audit_fractional(hd_counts, "HD"),
    audit_fractional(pd_counts, "PD")
  )
)

print(fractional_summary)
#Verifica se algum dataset possui valores não inteiros
if (any(count_summary$n_non_integer > 0)) {

  cat(
    "\nNOTE: valores fracionarios foram detectados. ",
    "Nenhum arredondamento foi realizado nesta etapa.\n",
    sep = ""
  )
}


# ------------------------------------------------------------------------------
# GeneIDs
# ------------------------------------------------------------------------------
section("6. AUDITORIA DOS GeneIDs")
#os IDs são "limpos", removendo espaços extras, transformando em texto e transformando valores ausentes em NA
ad_gene_id <- normalize_gene_id(ad$GeneID)
hd_gene_id <- normalize_gene_id(hd$GeneID)
pd_gene_id <- normalize_gene_id(pd$GeneID)

#Criamos um resumo com a quantidade de GeneIDs, quantos GeneIDs estão ausentes, verifica os valores ausentes e conta quantos valores diferentes existem e se há valores duplicados
gene_id_summary <- data.table(
  dataset = c("AD", "HD", "PD"),

  n_rows = c(
    length(ad_gene_id),
    length(hd_gene_id),
    length(pd_gene_id)
  ),

  n_missing = c(
    sum(is.na(ad_gene_id)),
    sum(is.na(hd_gene_id)),
    sum(is.na(pd_gene_id))
  ),

  n_unique = c(
    uniqueN(ad_gene_id[!is.na(ad_gene_id)]),
    uniqueN(hd_gene_id[!is.na(hd_gene_id)]),
    uniqueN(pd_gene_id[!is.na(pd_gene_id)])
  ),

  n_duplicated = c(
    sum(duplicated(ad_gene_id) & !is.na(ad_gene_id)),
    sum(duplicated(hd_gene_id) & !is.na(hd_gene_id)),
    sum(duplicated(pd_gene_id) & !is.na(pd_gene_id))
  )
)

print(gene_id_summary)

#Se qualquer dataset possuir GeneID ausente ou duplicado, a análise é interrompida
if (
  any(gene_id_summary$n_missing > 0) ||
  any(gene_id_summary$n_duplicated > 0)
) {

  stop(
    "GeneIDs ausentes ou duplicados foram encontrados.",
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# Simbolos do AD
# ------------------------------------------------------------------------------
section("7. SIMBOLOS DO AD")
#Auditoria dos GeneSymbol
ad_symbol <- clean_annotation_field(
  ad$GeneSymbol
)

#calculamos o total de linhas, quantos possuem símbolo e quantos não possuem, quantos símbolos diferentes existem e quantos aparecem repetidamente
ad_symbol_summary <- data.table(
  n_rows = length(ad_symbol),
  n_present = sum(!is.na(ad_symbol)),
  n_missing = sum(is.na(ad_symbol)),
  n_unique = uniqueN(ad_symbol[!is.na(ad_symbol)]),
  n_duplicated = sum(
    duplicated(ad_symbol) &
    !is.na(ad_symbol)
  )
)
print(ad_symbol_summary)


# ------------------------------------------------------------------------------
# Anotacao GRCh38.p13
# ------------------------------------------------------------------------------
section("8. ANOTACAO GRCh38.p13")
#Definindo as colunas obrigatórias
required_annotation <- c(
  "GeneID",
  "Symbol",
  "GeneType",
  "EnsemblGeneID"
)

#Verificar se alguma das colunas obrigatórias está faltando
missing_annotation <- setdiff(
  required_annotation,
  names(annotation)
)
#Se faltar alguma colunas, parar 
if (length(missing_annotation) > 0L) {

  stop(
    "Colunas ausentes na anotacao: ",
    paste(missing_annotation, collapse = ", "),
    call. = FALSE
  )
}

#Limpando as quatro colunas da anotação
annotation_gene_id <- normalize_gene_id(
  annotation$GeneID
)

annotation_symbol <- clean_annotation_field(
  annotation$Symbol
)

annotation_ensembl <- clean_annotation_field(
  annotation$EnsemblGeneID
)

annotation_gene_type <- clean_annotation_field(
  annotation$GeneType
)

annotation_summary <- data.table(
  #número total de registros da anotação
  n_rows = nrow(annotation),

  #quantos GeneIDs únicos existem
  n_unique_gene_id =
    uniqueN(
      annotation_gene_id[
        !is.na(annotation_gene_id)
      ]
    ),

  #quantos GeneIDs estão duplicados 
  n_duplicated_gene_id =
    sum(
      duplicated(annotation_gene_id) &
      !is.na(annotation_gene_id)
    ),

  #Quantas linhas possuem símbolo
  n_symbol_present =
    sum(!is.na(annotation_symbol)),
  #Quantas linhas possuem ou não Ensembl Gene ID
  n_ensembl_present =
    sum(!is.na(annotation_ensembl)),

  n_ensembl_missing =
    sum(is.na(annotation_ensembl)),
  #Quantas possuem ou não informação do tipo gênico
  n_gene_type_present =
    sum(!is.na(annotation_gene_type)),

  n_gene_type_missing =
    sum(is.na(annotation_gene_type))
)

print(annotation_summary)


# ------------------------------------------------------------------------------
# Universo gênico
# ------------------------------------------------------------------------------
section("9. UNIVERSO GENICO")
#Pegamos os GeneIDs válidos e únicos de cada dataset (remove as repetições)
genes_ad <- unique(
  ad_gene_id[!is.na(ad_gene_id)]
)

genes_hd <- unique(
  hd_gene_id[!is.na(hd_gene_id)]
)

genes_pd <- unique(
  pd_gene_id[!is.na(pd_gene_id)]
)
#genes existentes na anotação
annotation_ids <- unique(
  annotation_gene_id[!is.na(annotation_gene_id)]
)

#juntando os três vetores, removendo os genes repetidos e colocando os IDs em ordem
all_gene_ids <- sort(
  unique(
    c(
      genes_ad,
      genes_hd,
      genes_pd
    )
  )
)


gene_universe <- data.table(
  gene_id = all_gene_ids
)
#Modificamos as colunas. 
#Crie a coluna in_AD e coloque TRUE se aquele gene_id existir em genes_ad
gene_universe[
  ,
  in_AD := gene_id %in% genes_ad
]

gene_universe[
  ,
  in_HD := gene_id %in% genes_hd
]

gene_universe[
  ,
  in_PD := gene_id %in% genes_pd
]

gene_universe[
  ,
  in_annotation := gene_id %in% annotation_ids
]

#Colaremos um gene na sua classe a depender em qual dataset estiver presente. Gerando uma tabela com TRUE e FALSE
gene_universe[
  ,
  #Se o gene estiver em AD E HD E PD, classifique como "AD_HD_PD"
  presence_class := fifelse(
    in_AD & in_HD & in_PD,
    "AD_HD_PD",

    fifelse(
      in_AD & in_HD & !in_PD,
      "AD_HD",

      fifelse(
        in_AD & !in_HD & in_PD,
        "AD_PD",

        fifelse(
          !in_AD & in_HD & in_PD,
          "HD_PD",

          fifelse(
            in_AD & !in_HD & !in_PD,
            "AD_only",

            fifelse(
              !in_AD & in_HD & !in_PD,
              "HD_only",

              fifelse(
                !in_AD & !in_HD & in_PD,
                "PD_only",
                "other"
              )
            )
          )
        )
      )
    )
  )
]

#Contando quantos genes existem em cada classe
overlap_summary <- gene_universe[
  ,
  .N,
  by = presence_class #agrupe pela classe de presença
][
  order(presence_class)
]

#renomeamos N para n_genes
setnames(
  overlap_summary,
  "N",
  "n_genes"
)

print(overlap_summary)


cat(
  "\nGenes comuns AD/HD/PD:",
    #Conta os TRUE e obtemos o número de genes compartilhados por AD, HD e PD
  sum(
    gene_universe$in_AD &
    gene_universe$in_HD &
    gene_universe$in_PD
  ),
  "\n"
)

#Conte quantos genes estão em AD e não estão no arquivo de anotação
cat(
  "GeneIDs AD ausentes da anotacao:",
  sum(
    gene_universe$in_AD &
    !gene_universe$in_annotation
  ),
  "\n"
)


# ------------------------------------------------------------------------------
# Genes all-zero
# ------------------------------------------------------------------------------
section("10. GENES ALL-ZERO")
#Encontramos os genes que só possuem valores 0 ao somar o valor de cada linha
all_zero_summary <- data.table(
  dataset = c("AD", "HD", "PD"),

  n_genes = c(
    nrow(ad_counts),
    nrow(hd_counts),
    nrow(pd_counts)
  ),

  n_all_zero = c(
    sum(rowSums(ad_counts) == 0),
    sum(rowSums(hd_counts) == 0),
    sum(rowSums(pd_counts) == 0)
  )
)
#Calculando a porcentagem: %all-0= (n.all-0/n.genes)*100
all_zero_summary[
  ,
  pct_all_zero :=
    100 * n_all_zero / n_genes
]

print(all_zero_summary)


# Salvar a tabela do universo genico
fwrite(
  gene_universe,
  output_gene_universe,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Captura o que apareceria no console e grava em um arquivo
capture.output(
  {
    cat("RAW RNA-SEQ DATA AUDIT\n")
    cat("======================\n\n")

    cat(
      "Audit performed before count preparation, ",
      "gene harmonization, filtering, normalization ",
      "or differential expression analysis.\n\n",
      sep = ""
    )

    cat("1. FILE DIMENSIONS\n")
    cat("------------------\n")
    print(file_dimensions)
    cat("\nAD RPKM columns:", length(ad_rpkm_cols), "\n\n")


    cat("2. RAW COUNT SUMMARY\n")
    cat("--------------------\n")
    print(count_summary)
    cat("\n")


    cat("3. FRACTIONAL VALUES\n")
    cat("--------------------\n")
    print(fractional_summary)

    if (any(count_summary$n_non_integer > 0)) {

      cat(
        "\nNOTE: fractional values were detected. ",
        "No rounding was performed during audit.\n"
      )
    }

    cat("\n")


    cat("4. GENE ID SUMMARY\n")
    cat("------------------\n")
    print(gene_id_summary)
    cat("\n")


    cat("5. AD SYMBOL SUMMARY\n")
    cat("--------------------\n")
    print(ad_symbol_summary)
    cat("\n")


    cat("6. GRCh38.p13 ANNOTATION\n")
    cat("------------------------\n")
    print(annotation_summary)
    cat("\n")


    cat("7. GENE UNIVERSE\n")
    cat("----------------\n")
    print(overlap_summary)

    cat(
      "\nCommon AD/HD/PD GeneIDs:",
      sum(
        gene_universe$in_AD &
        gene_universe$in_HD &
        gene_universe$in_PD
      ),
      "\n"
    )

    cat(
      "AD GeneIDs absent from GRCh38.p13 annotation:",
      sum(
        gene_universe$in_AD &
        !gene_universe$in_annotation
      ),
      "\n\n"
    )


    cat("8. ALL-ZERO GENES\n")
    cat("-----------------\n")
    print(all_zero_summary)
    cat("\n")


    cat("9. AUDIT STATUS\n")
    cat("---------------\n")

    cat(
      "PASS: no missing, negative, or duplicated GeneIDs ",
      "were detected in the input matrices.\n"
    )

    cat(
      "PASS: no missing or negative raw-count values ",
      "were detected.\n"
    )

    if (any(count_summary$n_non_integer > 0)) {

      cat(
        "NOTE: fractional raw values were detected and ",
        "must be handled explicitly during count preparation ",
        "before DESeq2.\n"
      )
    }

    cat(
      "\nNo genes were removed, harmonized, filtered ",
      "or normalized during this audit.\n"
    )

  },
  file = output_report
)



# Mensagem final
section("11. AUDITORIA CONCLUIDA")

cat("Outputs permanentes:\n")
cat("  ", output_report, "\n")
cat("  ", output_gene_universe, "\n\n")

cat(
  "Nenhum gene foi removido, harmonizado, ",
  "filtrado ou normalizado.\n",
  sep = ""
)



