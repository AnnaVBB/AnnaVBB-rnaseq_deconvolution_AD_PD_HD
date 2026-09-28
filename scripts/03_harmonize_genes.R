#!/usr/bin/env Rscript

# ==============================================================================
# 03_harmonize_genes.R
#
# Objetivo:
#   Harmonizar a identidade gênica entre AD, HD e PD utilizando NCBI GeneID
#   como chave primária.
#
# Entradas:
#   results/counts/counts_AD.tsv
#   results/counts/counts_HD.tsv
#   results/counts/counts_PD.tsv
#   data/raw/Human.GRCh38.p13.annot.tsv.gz
#
# Este script:
#   1. lê os GeneIDs das três matrizes preparadas;
#   2. lê a anotação GRCh38.p13;
#   3. constrói a união de todos os GeneIDs;
#   4. classifica presença/ausência em AD, HD e PD;
#   5. incorpora Symbol, EnsemblGeneID e GeneType;
#   6. preserva, quando necessário, o GeneSymbol original de AD;
#   7. identifica genes comuns e específicos;
#   8. audita genes sem anotação;
#   9. cria uma tabela mestre de identidade gênica.
#
# IMPORTANTE:
#   - Nenhum gene é removido das matrizes de contagens.
#   - Nenhuma matriz é equalizada nesta etapa.
#   - Nenhum filtro de expressão é aplicado.
#   - Nenhuma normalização é realizada.
#   - GeneID permanece a chave primária.
#   - Symbol e EnsemblGeneID são anotações, não chaves primárias.
#
# Saídas:
#   results/harmonization/gene_master.tsv
#   results/harmonization/gene_presence_audit.tsv
#
# gene_master.tsv:
#   tabela mestre para uso nas etapas posteriores.
#
# gene_presence_audit.tsv:
#   genes que NÃO estão simultaneamente presentes nos três datasets,
#   permitindo auditar explicitamente o que ficaria fora de uma futura
#   interseção.
# ==============================================================================


# ==============================================================================
# 1. PACOTES
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

ad_counts_file <- "results/counts/counts_AD.tsv"
hd_counts_file <- "results/counts/counts_HD.tsv"
pd_counts_file <- "results/counts/counts_PD.tsv"

annotation_file <- "data/raw/Human.GRCh38.p13.annot.tsv.gz"

# Arquivo AD original é usado SOMENTE para recuperar GeneSymbol original.
ad_original_file <- "data/raw/GSE53697/GSE53697_RNAseq_AD.txt.gz"

output_dir <- "results/harmonization"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. VERIFICAR ENTRADAS
# ==============================================================================

input_files <- c(
  ad_counts_file,
  hd_counts_file,
  pd_counts_file,
  annotation_file,
  ad_original_file
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {

  stop(
    paste0(
      "Arquivos de entrada não encontrados:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


# ==============================================================================
# 4. FUNÇÕES AUXILIARES
# ==============================================================================

# ------------------------------------------------------------------------------
# Padronizar GeneID como character
# ------------------------------------------------------------------------------

standardize_gene_id <- function(x) {

  x <- as.character(x)
  x <- trimws(x)

  x[
    is.na(x) |
      x == "" |
      tolower(x) %in% c(
        "na",
        "n/a",
        "nan"
      )
  ] <- NA_character_

  x
}


# ------------------------------------------------------------------------------
# Padronizar campos de anotação
# ------------------------------------------------------------------------------

clean_annotation_field <- function(x) {

  x <- as.character(x)
  x <- trimws(x)

  x[
    is.na(x) |
      x == "" |
      x == "-" |
      tolower(x) %in% c(
        "na",
        "n/a",
        "nan",
        "null"
      )
  ] <- NA_character_

  x
}


# ------------------------------------------------------------------------------
# Validar conjunto de GeneIDs
# ------------------------------------------------------------------------------

validate_gene_ids <- function(
  gene_ids,
  dataset_name
) {

  if (anyNA(gene_ids)) {

    stop(
      dataset_name,
      ": existem GeneIDs ausentes."
    )
  }

  if (anyDuplicated(gene_ids)) {

    stop(
      dataset_name,
      ": existem GeneIDs duplicados."
    )
  }

  invisible(TRUE)
}


# ==============================================================================
# 5. LER GeneIDs DAS MATRIZES PREPARADAS
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("LEITURA DOS GeneIDs\n")
cat("============================================\n")


# Lemos somente GeneID.
# As contagens não precisam ser carregadas nesta etapa.

ad_genes <- fread(
  ad_counts_file,
  select = "GeneID"
)

hd_genes <- fread(
  hd_counts_file,
  select = "GeneID"
)

pd_genes <- fread(
  pd_counts_file,
  select = "GeneID"
)


ad_ids <- standardize_gene_id(
  ad_genes$GeneID
)

hd_ids <- standardize_gene_id(
  hd_genes$GeneID
)

pd_ids <- standardize_gene_id(
  pd_genes$GeneID
)


validate_gene_ids(
  ad_ids,
  "AD"
)

validate_gene_ids(
  hd_ids,
  "HD"
)

validate_gene_ids(
  pd_ids,
  "PD"
)


cat(
  "AD:",
  length(ad_ids),
  "GeneIDs\n"
)

cat(
  "HD:",
  length(hd_ids),
  "GeneIDs\n"
)

cat(
  "PD:",
  length(pd_ids),
  "GeneIDs\n"
)


# Valores esperados a partir das etapas 00 e 02.

if (length(ad_ids) != 19185) {
  stop(
    "AD: esperados 19185 GeneIDs; encontrados ",
    length(ad_ids),
    "."
  )
}

if (length(hd_ids) != 39376) {
  stop(
    "HD: esperados 39376 GeneIDs; encontrados ",
    length(hd_ids),
    "."
  )
}

if (length(pd_ids) != 39376) {
  stop(
    "PD: esperados 39376 GeneIDs; encontrados ",
    length(pd_ids),
    "."
  )
}


# ==============================================================================
# 6. AUDITAR RELAÇÃO HD / PD
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("COMPARAÇÃO DOS UNIVERSOS HD / PD\n")
cat("============================================\n")


hd_only <- setdiff(
  hd_ids,
  pd_ids
)

pd_only <- setdiff(
  pd_ids,
  hd_ids
)


cat(
  "HD não presentes em PD:",
  length(hd_only),
  "\n"
)

cat(
  "PD não presentes em HD:",
  length(pd_only),
  "\n"
)


if (
  length(hd_only) == 0 &&
    length(pd_only) == 0
) {

  cat(
    "HD e PD possuem o mesmo conjunto de GeneIDs.\n"
  )

} else {

  warning(
    "HD e PD não possuem conjuntos idênticos de GeneIDs."
  )
}


# ==============================================================================
# 7. CONSTRUIR UNIVERSO GÊNICO
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("CONSTRUÇÃO DO UNIVERSO GÊNICO\n")
cat("============================================\n")


all_gene_ids <- sort(
  unique(
    c(
      ad_ids,
      hd_ids,
      pd_ids
    )
  )
)


gene_presence <- data.table(
  gene_id = all_gene_ids
)


gene_presence[
  ,
  in_AD := gene_id %in% ad_ids
]

gene_presence[
  ,
  in_HD := gene_id %in% hd_ids
]

gene_presence[
  ,
  in_PD := gene_id %in% pd_ids
]


# ------------------------------------------------------------------------------
# Número de datasets em que cada gene aparece
# ------------------------------------------------------------------------------

gene_presence[
  ,
  n_datasets := (
    as.integer(in_AD) +
      as.integer(in_HD) +
      as.integer(in_PD)
  )
]


# ------------------------------------------------------------------------------
# Classe explícita de presença
# ------------------------------------------------------------------------------

gene_presence[
  ,
  presence_class := fcase(

    in_AD & in_HD & in_PD,
    "AD_HD_PD",

    in_AD & in_HD & !in_PD,
    "AD_HD",

    in_AD & !in_HD & in_PD,
    "AD_PD",

    !in_AD & in_HD & in_PD,
    "HD_PD",

    in_AD & !in_HD & !in_PD,
    "AD_only",

    !in_AD & in_HD & !in_PD,
    "HD_only",

    !in_AD & !in_HD & in_PD,
    "PD_only",

    default = "unexpected"
  )
]


# ==============================================================================
# 8. RESUMO DE PRESENÇA / AUSÊNCIA
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("PRESENÇA DOS GENES ENTRE DATASETS\n")
cat("============================================\n")


cat(
  "União AD ∪ HD ∪ PD:",
  nrow(gene_presence),
  "\n"
)

cat(
  "AD ∩ HD ∩ PD:",
  sum(
    gene_presence$presence_class == "AD_HD_PD"
  ),
  "\n"
)

cat(
  "AD-only:",
  sum(
    gene_presence$presence_class == "AD_only"
  ),
  "\n"
)

cat(
  "HD+PD, ausentes em AD:",
  sum(
    gene_presence$presence_class == "HD_PD"
  ),
  "\n"
)

cat(
  "HD-only:",
  sum(
    gene_presence$presence_class == "HD_only"
  ),
  "\n"
)

cat(
  "PD-only:",
  sum(
    gene_presence$presence_class == "PD_only"
  ),
  "\n"
)

cat(
  "AD+HD, ausentes em PD:",
  sum(
    gene_presence$presence_class == "AD_HD"
  ),
  "\n"
)

cat(
  "AD+PD, ausentes em HD:",
  sum(
    gene_presence$presence_class == "AD_PD"
  ),
  "\n"
)


# Valores já observados no 00.
# Aqui funcionam como validação de continuidade do pipeline.

if (nrow(gene_presence) != 39620) {

  stop(
    "União esperada = 39620; encontrada = ",
    nrow(gene_presence),
    "."
  )
}


if (
  sum(
    gene_presence$presence_class == "AD_HD_PD"
  ) != 18941
) {

  stop(
    "Interseção tripla esperada = 18941."
  )
}


if (
  sum(
    gene_presence$presence_class == "AD_only"
  ) != 244
) {

  stop(
    "Esperados 244 GeneIDs AD-only."
  )
}


if (
  sum(
    gene_presence$presence_class == "HD_PD"
  ) != 20435
) {

  stop(
    "Esperados 20435 GeneIDs presentes em HD/PD e ausentes em AD."
  )
}


# ==============================================================================
# 9. LER ANOTAÇÃO GRCh38.p13
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("LEITURA DA ANOTAÇÃO GRCh38.p13\n")
cat("============================================\n")


annotation <- fread(
  annotation_file,
  check.names = FALSE
)


required_annotation_columns <- c(
  "GeneID",
  "Symbol",
  "GeneType",
  "EnsemblGeneID"
)


missing_annotation_columns <- setdiff(
  required_annotation_columns,
  names(annotation)
)


if (length(missing_annotation_columns) > 0) {

  stop(
    paste0(
      "Colunas ausentes na anotação: ",
      paste(
        missing_annotation_columns,
        collapse = ", "
      )
    )
  )
}


annotation <- annotation[
  ,
  .(
    gene_id = standardize_gene_id(
      GeneID
    ),

    symbol_annotation = clean_annotation_field(
      Symbol
    ),

    ensembl_id = clean_annotation_field(
      EnsemblGeneID
    ),

    gene_type = clean_annotation_field(
      GeneType
    )
  )
]


if (anyNA(annotation$gene_id)) {

  stop(
    "Existem GeneIDs ausentes na anotação."
  )
}


if (anyDuplicated(annotation$gene_id)) {

  stop(
    "Existem GeneIDs duplicados na anotação."
  )
}


cat(
  "Genes na anotação:",
  nrow(annotation),
  "\n"
)

cat(
  "Symbols disponíveis:",
  sum(
    !is.na(annotation$symbol_annotation)
  ),
  "/",
  nrow(annotation),
  "\n"
)

cat(
  "EnsemblGeneID disponíveis:",
  sum(
    !is.na(annotation$ensembl_id)
  ),
  "/",
  nrow(annotation),
  "\n"
)

cat(
  "GeneType disponíveis:",
  sum(
    !is.na(annotation$gene_type)
  ),
  "/",
  nrow(annotation),
  "\n"
)


# ==============================================================================
# 10. RECUPERAR GeneSymbol ORIGINAL DE AD
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("RECUPERAÇÃO DO GeneSymbol ORIGINAL DE AD\n")
cat("============================================\n")


ad_original_annotation <- fread(
  ad_original_file,
  select = c(
    "GeneID",
    "GeneSymbol"
  ),
  check.names = FALSE
)


ad_original_annotation <- ad_original_annotation[
  ,
  .(
    gene_id = standardize_gene_id(
      GeneID
    ),

    symbol_AD_original = clean_annotation_field(
      GeneSymbol
    )
  )
]


if (anyNA(ad_original_annotation$gene_id)) {

  stop(
    "AD original: existem GeneIDs ausentes."
  )
}


if (anyDuplicated(ad_original_annotation$gene_id)) {

  stop(
    "AD original: existem GeneIDs duplicados."
  )
}


if (nrow(ad_original_annotation) != 19185) {

  stop(
    "AD original deveria possuir 19185 genes."
  )
}


cat(
  "GeneSymbols AD disponíveis:",
  sum(
    !is.na(
      ad_original_annotation$symbol_AD_original
    )
  ),
  "/",
  nrow(ad_original_annotation),
  "\n"
)


# ==============================================================================
# 11. CONSTRUIR TABELA MESTRE
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("CONSTRUÇÃO DA TABELA MESTRE\n")
cat("============================================\n")


gene_master <- merge(
  gene_presence,
  annotation,
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)


gene_master <- merge(
  gene_master,
  ad_original_annotation,
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------------------------
# Symbol final
#
# Prioridade:
#   1. Symbol da anotação GRCh38.p13
#   2. GeneSymbol original de AD
#
# Isso permite preservar uma nomenclatura para parte dos genes AD-only
# que não existem na anotação GRCh38.p13 utilizada para HD/PD.
#
# A origem do símbolo permanece explicitamente registrada.
# ------------------------------------------------------------------------------

gene_master[
  ,
  symbol := fcoalesce(
    symbol_annotation,
    symbol_AD_original
  )
]


gene_master[
  ,
  symbol_source := fcase(

    !is.na(symbol_annotation),
    "GRCh38.p13_annotation",

    is.na(symbol_annotation) &
      !is.na(symbol_AD_original),
    "AD_original",

    default = NA_character_
  )
]


# ------------------------------------------------------------------------------
# Cobertura da anotação
# ------------------------------------------------------------------------------

gene_master[
  ,
  in_annotation := !is.na(symbol_annotation) |
    !is.na(ensembl_id) |
    !is.na(gene_type)
]


# ------------------------------------------------------------------------------
# Ordenação das colunas
# ------------------------------------------------------------------------------

setcolorder(
  gene_master,
  c(
    "gene_id",
    "symbol",
    "ensembl_id",
    "gene_type",
    "symbol_source",
    "symbol_annotation",
    "symbol_AD_original",
    "in_AD",
    "in_HD",
    "in_PD",
    "n_datasets",
    "presence_class",
    "in_annotation"
  )
)


# Ordenação determinística por GeneID numérico quando possível.

gene_master[
  ,
  gene_id_numeric__ := suppressWarnings(
    as.numeric(gene_id)
  )
]

setorder(
  gene_master,
  gene_id_numeric__,
  gene_id
)

gene_master[
  ,
  gene_id_numeric__ := NULL
]


# ==============================================================================
# 12. AUDITORIA DA ANOTAÇÃO
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("COBERTURA DA ANOTAÇÃO\n")
cat("============================================\n")


cat(
  "Universo total:",
  nrow(gene_master),
  "\n"
)

cat(
  "Com Symbol final:",
  sum(
    !is.na(gene_master$symbol)
  ),
  "/",
  nrow(gene_master),
  "\n"
)

cat(
  "Com EnsemblGeneID:",
  sum(
    !is.na(gene_master$ensembl_id)
  ),
  "/",
  nrow(gene_master),
  "\n"
)

cat(
  "Com GeneType:",
  sum(
    !is.na(gene_master$gene_type)
  ),
  "/",
  nrow(gene_master),
  "\n"
)

cat(
  "Sem correspondência na anotação GRCh38.p13:",
  sum(
    !gene_master$in_annotation
  ),
  "\n"
)


# ------------------------------------------------------------------------------
# Auditoria específica dos 244 genes AD-only
# ------------------------------------------------------------------------------

ad_only_table <- gene_master[
  presence_class == "AD_only"
]


cat("\nAD-only:\n")

cat(
  "  Total:",
  nrow(ad_only_table),
  "\n"
)

cat(
  "  Com Symbol final:",
  sum(
    !is.na(ad_only_table$symbol)
  ),
  "\n"
)

cat(
  "  Com Symbol GRCh38.p13:",
  sum(
    !is.na(ad_only_table$symbol_annotation)
  ),
  "\n"
)

cat(
  "  Com Symbol original AD:",
  sum(
    !is.na(ad_only_table$symbol_AD_original)
  ),
  "\n"
)

cat(
  "  Com EnsemblGeneID:",
  sum(
    !is.na(ad_only_table$ensembl_id)
  ),
  "\n"
)

cat(
  "  Com GeneType:",
  sum(
    !is.na(ad_only_table$gene_type)
  ),
  "\n"
)


# ------------------------------------------------------------------------------
# Auditoria dos genes HD/PD ausentes em AD
# ------------------------------------------------------------------------------

hd_pd_not_ad <- gene_master[
  presence_class == "HD_PD"
]


cat("\nHD+PD ausentes em AD:\n")

cat(
  "  Total:",
  nrow(hd_pd_not_ad),
  "\n"
)

cat(
  "  Com Symbol:",
  sum(
    !is.na(hd_pd_not_ad$symbol)
  ),
  "\n"
)

cat(
  "  Com EnsemblGeneID:",
  sum(
    !is.na(hd_pd_not_ad$ensembl_id)
  ),
  "\n"
)

cat(
  "  Com GeneType:",
  sum(
    !is.na(hd_pd_not_ad$gene_type)
  ),
  "\n"
)


# ==============================================================================
# 13. COMPARAR SYMBOL AD ORIGINAL vs ANOTAÇÃO
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("COMPARAÇÃO DE SYMBOLS AD\n")
cat("============================================\n")


ad_symbol_comparison <- gene_master[
  in_AD &
    !is.na(symbol_annotation) &
    !is.na(symbol_AD_original)
]


n_symbol_compared <- nrow(
  ad_symbol_comparison
)


n_symbol_same <- sum(
  ad_symbol_comparison$symbol_annotation ==
    ad_symbol_comparison$symbol_AD_original
)


n_symbol_different <- sum(
  ad_symbol_comparison$symbol_annotation !=
    ad_symbol_comparison$symbol_AD_original
)


cat(
  "Comparáveis:",
  n_symbol_compared,
  "\n"
)

cat(
  "Symbols idênticos:",
  n_symbol_same,
  "\n"
)

cat(
  "Symbols diferentes:",
  n_symbol_different,
  "\n"
)


# IMPORTANTE:
# diferenças de Symbol não são automaticamente corrigidas nem descartadas.
# A tabela mestre preserva ambas as versões para rastreabilidade.


# ==============================================================================
# 14. CRIAR AUDITORIA DE GENES NÃO COMUNS
# ==============================================================================

gene_presence_audit <- gene_master[
  presence_class != "AD_HD_PD"
]


# Essa tabela deve conter todos os genes que seriam afetados por uma
# eventual restrição futura ao universo comum.

expected_non_common <- (
  nrow(gene_master) -
    sum(
      gene_master$presence_class == "AD_HD_PD"
    )
)


if (
  nrow(gene_presence_audit) !=
    expected_non_common
) {

  stop(
    "Erro na construção da auditoria de genes não comuns."
  )
}


# ==============================================================================
# 15. VALIDAÇÕES FINAIS
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("VALIDAÇÕES FINAIS\n")
cat("============================================\n")


if (nrow(gene_master) != 39620) {

  stop(
    "gene_master deveria possuir 39620 genes."
  )
}


if (anyNA(gene_master$gene_id)) {

  stop(
    "gene_master contém GeneID ausente."
  )
}


if (anyDuplicated(gene_master$gene_id)) {

  stop(
    "gene_master contém GeneID duplicado."
  )
}


if (
  sum(
    gene_master$in_AD
  ) != 19185
) {

  stop(
    "Número inesperado de genes marcados como presentes em AD."
  )
}


if (
  sum(
    gene_master$in_HD
  ) != 39376
) {

  stop(
    "Número inesperado de genes marcados como presentes em HD."
  )
}


if (
  sum(
    gene_master$in_PD
  ) != 39376
) {

  stop(
    "Número inesperado de genes marcados como presentes em PD."
  )
}


if (
  sum(
    gene_master$presence_class == "AD_HD_PD"
  ) != 18941
) {

  stop(
    "Número inesperado de genes comuns aos três datasets."
  )
}


cat(
  "GeneIDs únicos:",
  nrow(gene_master),
  "/",
  nrow(gene_master),
  "\n"
)

cat(
  "Genes AD:",
  sum(gene_master$in_AD),
  "/ 19185\n"
)

cat(
  "Genes HD:",
  sum(gene_master$in_HD),
  "/ 39376\n"
)

cat(
  "Genes PD:",
  sum(gene_master$in_PD),
  "/ 39376\n"
)

cat(
  "Genes comuns:",
  sum(
    gene_master$presence_class == "AD_HD_PD"
  ),
  "/ 18941\n"
)


# ==============================================================================
# 16. GRAVAR RESULTADOS
# ==============================================================================

gene_master_file <- file.path(
  output_dir,
  "gene_master.tsv"
)

gene_presence_audit_file <- file.path(
  output_dir,
  "gene_presence_audit.tsv"
)


fwrite(
  gene_master,
  gene_master_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  gene_presence_audit,
  gene_presence_audit_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 17. FINALIZAÇÃO
# ==============================================================================

cat("\n")
cat("============================================\n")
cat("03_harmonize_genes.R CONCLUÍDO\n")
cat("============================================\n")


cat(
  "Universo total:",
  nrow(gene_master),
  "genes\n"
)

cat(
  "AD:",
  sum(gene_master$in_AD),
  "genes\n"
)

cat(
  "HD:",
  sum(gene_master$in_HD),
  "genes\n"
)

cat(
  "PD:",
  sum(gene_master$in_PD),
  "genes\n"
)

cat(
  "Comuns aos três:",
  sum(
    gene_master$presence_class == "AD_HD_PD"
  ),
  "genes\n"
)

cat(
  "AD-only:",
  sum(
    gene_master$presence_class == "AD_only"
  ),
  "genes\n"
)

cat(
  "HD+PD ausentes em AD:",
  sum(
    gene_master$presence_class == "HD_PD"
  ),
  "genes\n"
)


cat("\nArquivos gerados:\n")

cat(
  " -",
  gene_master_file,
  "\n"
)

cat(
  " -",
  gene_presence_audit_file,
  "\n"
)


cat("\nNenhum gene foi removido das matrizes de contagens.\n")
cat("Nenhuma matriz foi equalizada.\n")
cat("GeneID foi mantido como chave primária.\n")
cat("Symbol e EnsemblGeneID foram mantidos como anotações auxiliares.\n")
cat("A definição do universo final será feita somente após esta auditoria.\n")