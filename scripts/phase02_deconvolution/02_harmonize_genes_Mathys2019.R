#!/usr/bin/env Rscript

# ==============================================================================
# 02_harmonize_genes_Mathys2019.R
# ==============================================================================
#
# FASE 02 — DECONVOLUÇÃO CELULAR
# ETAPA 02 — HARMONIZAÇÃO GÊNICA DA REFERÊNCIA MATHYS 2019
#
# Objetivo
# --------
# Converter os identificadores gênicos da referência Mathys 2019, atualmente
# representados por gene symbols, para NCBI GeneID / Entrez Gene ID.
#
# Isso coloca a referência no mesmo tipo de identificador usado no restante
# do projeto e evita fazer a integração bulk × referência diretamente por
# gene symbol.
#
#
# ESTRATÉGIA DE MAPEAMENTO
# ------------------------
# O script usa org.Hs.eg.db / AnnotationDbi e aplica uma política conservadora:
#
#   1. tentativa principal:
#        SYMBOL -> ENTREZID
#
#   2. somente para símbolos sem match direto:
#        ALIAS -> ENTREZID
#
#   3. um mapeamento só é aceito automaticamente quando existe exatamente
#      UM GeneID único para o símbolo/alias.
#
# Portanto:
#
#   1 símbolo -> 1 GeneID
#       ACEITO
#
#   1 símbolo -> 0 GeneID
#       NÃO RESOLVIDO
#
#   1 símbolo -> >1 GeneID
#       AMBÍGUO; NÃO É RESOLVIDO AUTOMATICAMENTE
#
#
# DUPLICAÇÕES APÓS O MAPEAMENTO
# -----------------------------
# Símbolos diferentes podem eventualmente convergir para o mesmo GeneID.
#
# Este script:
#
#   - identifica essas colisões;
#   - registra todas;
#   - NÃO remove nenhuma linha;
#   - NÃO soma counts nesta etapa.
#
# O colapso por GeneID será feito no:
#
#   03_build_reference_Mathys2019.R
#
# pela soma das linhas que forem mapeadas para o mesmo GeneID.
#
#
# IMPORTANTE
# ----------
# Este script NÃO:
#
#   - modifica a matriz de expressão;
#   - seleciona núcleos;
#   - cria pseudobulk;
#   - calcula CPM;
#   - seleciona marcadores;
#   - intersecta ainda Mathys com ROSMAP/HD/PD;
#   - usa "multiVals = first";
#   - escolhe arbitrariamente um GeneID em mappings ambíguos.
#
#
# INPUT PRINCIPAL
# ---------------
# data/reference/mathys2019/raw/filtered_gene_row_names.txt
#
# O arquivo já foi auditado no 00:
#   - 17.926 linhas;
#   - correspondência 1:1 com as linhas da matriz;
#   - conteúdo compatível com gene symbols.
#
#
# OUTPUT PRINCIPAL
# ----------------
# data/reference/mathys2019/phase02_prepared/gene_harmonization/
#   mathys2019_gene_mapping.tsv
#
# Esse arquivo preserva matrix_row_index, permitindo que o script 03 aplique
# o mapa diretamente às linhas da matriz original sem alterar a ordem.
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

if (length(missing_packages) > 0L) {

  stop(
    "Pacotes necessários não instalados: ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\nNão altere o renv automaticamente. Instale somente os pacotes ausentes ",
    "no ambiente R atual e execute novamente."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})


# ==============================================================================
# 1. PARÂMETROS PREDEFINIDOS
# ==============================================================================

EXPECTED_MATHYS_GENE_ROWS <- 17926L

# Métodos considerados resolvidos.
RESOLVED_METHODS <- c(
  "direct_symbol",
  "rescued_unique_alias"
)


# ==============================================================================
# 2. CAMINHOS
# ==============================================================================

raw_dir <- file.path(
  "data",
  "reference",
  "mathys2019",
  "raw"
)

prepared_dir <- file.path(
  "data",
  "reference",
  "mathys2019",
  "phase02_prepared"
)

harmonized_dir <- file.path(
  prepared_dir,
  "gene_harmonization"
)

result_dir <- file.path(
  "results",
  "deconvolution",
  "mathys2019",
  "gene_harmonization"
)

dir.create(
  harmonized_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------

gene_file <- file.path(
  raw_dir,
  "filtered_gene_row_names.txt"
)

# Output do 01 usado somente para garantir que a etapa anterior foi concluída.
step01_status_file <- file.path(
  "results",
  "deconvolution",
  "mathys2019",
  "metadata",
  "01_prepare_metadata_status.tsv"
)


# ------------------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------------------

mapping_output_file <- file.path(
  harmonized_dir,
  "mathys2019_gene_mapping.tsv"
)

collapse_map_file <- file.path(
  harmonized_dir,
  "mathys2019_geneid_collapse_map.tsv"
)

resolved_rows_file <- file.path(
  harmonized_dir,
  "mathys2019_resolved_gene_rows.rds"
)

log_file <- file.path(
  result_dir,
  "02_harmonize_genes_Mathys2019.txt"
)

session_file <- file.path(
  result_dir,
  "02_harmonize_genes_Mathys2019_sessionInfo.txt"
)


# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================

clean_chr <- function(x) {

  x <- trimws(
    as.character(x)
  )

  x[
    is.na(x) |
      x == "" |
      toupper(x) == "NA"
  ] <- NA_character_

  x
}


collapse_values <- function(x) {

  x <- clean_chr(x)

  x <- sort(
    unique(
      x[
        !is.na(x)
      ]
    )
  )

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste(
    x,
    collapse = ";"
  )
}


count_unique_nonmissing <- function(x) {

  x <- clean_chr(x)

  uniqueN(
    x[
      !is.na(x)
    ]
  )
}


make_check <- function(
  check,
  observed,
  expected,
  status,
  note = NA_character_
) {

  data.table(
    check = check,
    observed = as.character(observed),
    expected = as.character(expected),
    status = status,
    note = note
  )
}


# ------------------------------------------------------------------------------
# Verificar se a etapa 01 terminou sem FAIL.
# ------------------------------------------------------------------------------

validate_previous_step <- function() {

  if (!file.exists(step01_status_file)) {

    stop(
      "Output da etapa 01 não encontrado:\n",
      step01_status_file,
      "\nExecute primeiro 01_prepare_metadata_Mathys2019.R."
    )
  }

  status01 <- fread(
    step01_status_file,
    na.strings = c(
      "NA",
      ""
    )
  )

  if (!"status" %in% names(status01)) {

    stop(
      "Arquivo de status da etapa 01 não possui a coluna 'status'."
    )
  }

  if (
    any(
      status01$status == "FAIL",
      na.rm = TRUE
    )
  ) {

    stop(
      "A etapa 01 contém check(s) FAIL. ",
      "A harmonização gênica foi interrompida."
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# Resumir mapeamento SYMBOL -> ENTREZID.
#
# AnnotationDbi::select pode retornar múltiplas linhas por chave porque uma
# chave pode estar associada a mais de uma anotação. Não escolhemos "a primeira".
# Primeiro contamos GeneIDs únicos por símbolo.
# ------------------------------------------------------------------------------

summarize_direct_mapping <- function(raw_map) {

  raw_map <- as.data.table(
    raw_map
  )

  if (!"SYMBOL" %in% names(raw_map)) {

    stop(
      "Resultado de AnnotationDbi::select por SYMBOL sem coluna SYMBOL."
    )
  }

  raw_map[
    ,
    query_symbol := clean_chr(
      SYMBOL
    )
  ]

  raw_map[
    ,
    entrez_clean := clean_chr(
      ENTREZID
    )
  ]

  raw_map[
    ,
    ensembl_clean := clean_chr(
      ENSEMBL
    )
  ]

  raw_map[
    ,
    gene_name_clean := clean_chr(
      GENENAME
    )
  ]

  raw_map[
    ,
    .(
      direct_n_entrez = count_unique_nonmissing(
        entrez_clean
      ),
      direct_entrez_ids = collapse_values(
        entrez_clean
      ),
      direct_ensembl_ids = collapse_values(
        ensembl_clean
      ),
      direct_gene_names = collapse_values(
        gene_name_clean
      )
    ),
    by = .(
      original_symbol = query_symbol
    )
  ]
}


# ------------------------------------------------------------------------------
# Resumir mapeamento ALIAS -> ENTREZID.
#
# Essa segunda camada é aplicada SOMENTE aos símbolos que não tiveram match
# direto por SYMBOL. O alias só é aceito se apontar para exatamente um GeneID.
# ------------------------------------------------------------------------------

summarize_alias_mapping <- function(raw_map) {

  raw_map <- as.data.table(
    raw_map
  )

  if (!"ALIAS" %in% names(raw_map)) {

    stop(
      "Resultado de AnnotationDbi::select por ALIAS sem coluna ALIAS."
    )
  }

  raw_map[
    ,
    query_alias := clean_chr(
      ALIAS
    )
  ]

  raw_map[
    ,
    official_symbol_clean := clean_chr(
      SYMBOL
    )
  ]

  raw_map[
    ,
    entrez_clean := clean_chr(
      ENTREZID
    )
  ]

  raw_map[
    ,
    ensembl_clean := clean_chr(
      ENSEMBL
    )
  ]

  raw_map[
    ,
    gene_name_clean := clean_chr(
      GENENAME
    )
  ]

  raw_map[
    ,
    .(
      alias_n_entrez = count_unique_nonmissing(
        entrez_clean
      ),
      alias_entrez_ids = collapse_values(
        entrez_clean
      ),
      alias_official_symbols = collapse_values(
        official_symbol_clean
      ),
      alias_ensembl_ids = collapse_values(
        ensembl_clean
      ),
      alias_gene_names = collapse_values(
        gene_name_clean
      )
    ),
    by = .(
      original_symbol = query_alias
    )
  ]
}


# ==============================================================================
# 4. FUNÇÃO PRINCIPAL
# ==============================================================================

run_harmonization <- function() {


  # ============================================================================
  # 4.1 CABEÇALHO E PRÉ-REQUISITOS
  # ============================================================================

  cat(
    "============================================================\n"
  )

  cat(
    "FASE 02 — 02_harmonize_genes_Mathys2019.R\n"
  )

  cat(
    "SYMBOL -> NCBI GeneID / ENTREZID\n"
  )

  cat(
    "============================================================\n\n"
  )

  validate_previous_step()

  cat(
    "Etapa 01 validada: nenhum check FAIL.\n"
  )

  if (!file.exists(gene_file)) {

    stop(
      "Arquivo de genes não encontrado:\n",
      gene_file
    )
  }


  # ============================================================================
  # 4.2 CARREGAR GENES MATHYS
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "1. CARREGAMENTO DOS GENES MATHYS\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  genes <- fread(
    gene_file,
    header = FALSE,
    fill = TRUE,
    na.strings = c(
      "NA",
      ""
    )
  )

  if (ncol(genes) != 1L) {

    stop(
      "Esperávamos uma única coluna em filtered_gene_row_names.txt; ",
      "foram encontradas ",
      ncol(genes),
      "."
    )
  }

  setnames(
    genes,
    "V1",
    "original_symbol"
  )

  genes[
    ,
    original_symbol := clean_chr(
      original_symbol
    )
  ]

  genes[
    ,
    matrix_row_index := seq_len(.N)
  ]

  setcolorder(
    genes,
    c(
      "matrix_row_index",
      "original_symbol"
    )
  )

  cat(
    "Linhas no arquivo de genes: ",
    nrow(genes),
    "\n",
    sep = ""
  )

  if (
    nrow(genes) !=
      EXPECTED_MATHYS_GENE_ROWS
  ) {

    stop(
      "Esperávamos ",
      EXPECTED_MATHYS_GENE_ROWS,
      " linhas com base na auditoria 00, mas encontramos ",
      nrow(genes),
      "."
    )
  }

  n_missing_symbols <- sum(
    is.na(
      genes$original_symbol
    )
  )

  n_unique_symbols <- uniqueN(
    genes$original_symbol[
      !is.na(
        genes$original_symbol
      )
    ]
  )

  n_duplicated_symbol_rows <- sum(
    duplicated(
      genes$original_symbol[
        !is.na(
          genes$original_symbol
        )
      ]
    )
  )

  cat(
    "Símbolos ausentes: ",
    n_missing_symbols,
    "\n",
    sep = ""
  )

  cat(
    "Símbolos únicos: ",
    n_unique_symbols,
    "\n",
    sep = ""
  )

  cat(
    "Linhas com símbolo duplicado: ",
    n_duplicated_symbol_rows,
    "\n",
    sep = ""
  )

  if (n_missing_symbols > 0L) {

    stop(
      "Existem símbolos gênicos ausentes. ",
      "A etapa foi interrompida."
    )
  }


  # ============================================================================
  # 4.3 AUDITAR VERSÕES DA ANOTAÇÃO
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "2. VERSÕES DA ANOTAÇÃO\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  annotation_versions <- data.table(
    resource = c(
      "R",
      "AnnotationDbi",
      "org.Hs.eg.db"
    ),
    version = c(
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

  fwrite(
    annotation_versions,
    file.path(
      result_dir,
      "02_annotation_versions.tsv"
    ),
    sep = "\t"
  )

  print(
    annotation_versions
  )

  cat(
    "\nKeytypes disponíveis em org.Hs.eg.db incluem SYMBOL: ",
    "SYMBOL" %in%
      AnnotationDbi::keytypes(
        org.Hs.eg.db::org.Hs.eg.db
      ),
    "\n",
    sep = ""
  )

  cat(
    "Keytypes disponíveis em org.Hs.eg.db incluem ALIAS: ",
    "ALIAS" %in%
      AnnotationDbi::keytypes(
        org.Hs.eg.db::org.Hs.eg.db
      ),
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.4 MAPEAMENTO DIRETO SYMBOL -> ENTREZID
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "3. MAPEAMENTO DIRETO SYMBOL -> ENTREZID\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  symbol_keys <- sort(
    unique(
      genes$original_symbol
    )
  )

  cat(
    "Símbolos submetidos ao mapeamento direto: ",
    length(symbol_keys),
    "\n",
    sep = ""
  )

  direct_raw <- suppressMessages(
    AnnotationDbi::select(
      x = org.Hs.eg.db::org.Hs.eg.db,
      keys = symbol_keys,
      keytype = "SYMBOL",
      columns = c(
        "ENTREZID",
        "ENSEMBL",
        "GENENAME"
      )
    )
  )

  direct_summary <- summarize_direct_mapping(
    direct_raw
  )

  # Garantir que TODOS os símbolos de entrada estejam representados no resumo.
  direct_summary <- merge(
    data.table(
      original_symbol = symbol_keys
    ),
    direct_summary,
    by = "original_symbol",
    all.x = TRUE,
    sort = FALSE
  )

  direct_summary[
    is.na(
      direct_n_entrez
    ),
    direct_n_entrez := 0L
  ]

  n_direct_unique <- sum(
    direct_summary$direct_n_entrez == 1L
  )

  n_direct_ambiguous <- sum(
    direct_summary$direct_n_entrez > 1L
  )

  n_direct_unmapped <- sum(
    direct_summary$direct_n_entrez == 0L
  )

  cat(
    "Mapeamento direto único: ",
    n_direct_unique,
    "\n",
    sep = ""
  )

  cat(
    "Mapeamento direto ambíguo (>1 GeneID): ",
    n_direct_ambiguous,
    "\n",
    sep = ""
  )

  cat(
    "Sem GeneID por SYMBOL: ",
    n_direct_unmapped,
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.5 RESGATE CONSERVADOR POR ALIAS
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "4. RESGATE CONSERVADOR ALIAS -> ENTREZID\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  alias_candidates <- direct_summary[
    direct_n_entrez == 0L,
    original_symbol
  ]

  cat(
    "Símbolos candidatos ao resgate por alias: ",
    length(alias_candidates),
    "\n",
    sep = ""
  )

  if (length(alias_candidates) > 0L) {

    alias_raw <- suppressMessages(
      AnnotationDbi::select(
        x = org.Hs.eg.db::org.Hs.eg.db,
        keys = alias_candidates,
        keytype = "ALIAS",
        columns = c(
          "SYMBOL",
          "ENTREZID",
          "ENSEMBL",
          "GENENAME"
        )
      )
    )

    alias_summary <- summarize_alias_mapping(
      alias_raw
    )

    alias_summary <- merge(
      data.table(
        original_symbol = alias_candidates
      ),
      alias_summary,
      by = "original_symbol",
      all.x = TRUE,
      sort = FALSE
    )

    alias_summary[
      is.na(
        alias_n_entrez
      ),
      alias_n_entrez := 0L
    ]

  } else {

    alias_summary <- data.table(
      original_symbol = character(),
      alias_n_entrez = integer(),
      alias_entrez_ids = character(),
      alias_official_symbols = character(),
      alias_ensembl_ids = character(),
      alias_gene_names = character()
    )
  }

  n_alias_unique <- sum(
    alias_summary$alias_n_entrez == 1L
  )

  n_alias_ambiguous <- sum(
    alias_summary$alias_n_entrez > 1L
  )

  n_alias_unmapped <- sum(
    alias_summary$alias_n_entrez == 0L
  )

  cat(
    "Resgatados por alias único: ",
    n_alias_unique,
    "\n",
    sep = ""
  )

  cat(
    "Alias ambíguo (>1 GeneID): ",
    n_alias_ambiguous,
    "\n",
    sep = ""
  )

  cat(
    "Ainda não resolvidos após alias: ",
    n_alias_unmapped,
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.6 CONSTRUIR MAPA FINAL POR LINHA DA MATRIZ
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "5. MAPA FINAL POR LINHA DA MATRIZ\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  mapping <- merge(
    copy(genes),
    direct_summary,
    by = "original_symbol",
    all.x = TRUE,
    sort = FALSE
  )

  if (nrow(alias_summary) > 0L) {

    mapping <- merge(
      mapping,
      alias_summary,
      by = "original_symbol",
      all.x = TRUE,
      sort = FALSE
    )

  } else {

    mapping[
      ,
      `:=`(
        alias_n_entrez = NA_integer_,
        alias_entrez_ids = NA_character_,
        alias_official_symbols = NA_character_,
        alias_ensembl_ids = NA_character_,
        alias_gene_names = NA_character_
      )
    ]
  }

  # Restaurar rigorosamente a ordem original das linhas da matriz.
  setorder(
    mapping,
    matrix_row_index
  )

  if (
    !identical(
      mapping$matrix_row_index,
      seq_len(
        nrow(genes)
      )
    )
  ) {

    stop(
      "A ordem das linhas foi alterada durante a harmonização."
    )
  }

  # Número de ocorrências de cada símbolo no input.
  mapping[
    ,
    input_symbol_n_rows := .N,
    by = original_symbol
  ]


  # --------------------------------------------------------------------------
  # Inicializar classificação.
  # --------------------------------------------------------------------------

  mapping[
    ,
    `:=`(
      gene_id = NA_character_,
      annotation_symbol = NA_character_,
      annotation_ensembl = NA_character_,
      annotation_gene_name = NA_character_,
      mapping_method = "unmapped"
    )
  ]


  # --------------------------------------------------------------------------
  # 1) Direct unique SYMBOL -> ENTREZID
  # --------------------------------------------------------------------------

  direct_idx <- which(
    mapping$direct_n_entrez == 1L
  )

  if (length(direct_idx) > 0L) {

    mapping[
      direct_idx,
      `:=`(
        gene_id = direct_entrez_ids,
        annotation_symbol = original_symbol,
        annotation_ensembl = direct_ensembl_ids,
        annotation_gene_name = direct_gene_names,
        mapping_method = "direct_symbol"
      )
    ]
  }


  # --------------------------------------------------------------------------
  # 2) Direct ambiguous
  # --------------------------------------------------------------------------

  direct_ambiguous_idx <- which(
    mapping$direct_n_entrez > 1L
  )

  if (length(direct_ambiguous_idx) > 0L) {

    mapping[
      direct_ambiguous_idx,
      mapping_method := "ambiguous_direct"
    ]
  }


  # --------------------------------------------------------------------------
  # 3) Alias unique — somente quando não houve direct mapping.
  # --------------------------------------------------------------------------

  alias_unique_idx <- which(
    mapping$direct_n_entrez == 0L &
      !is.na(
        mapping$alias_n_entrez
      ) &
      mapping$alias_n_entrez == 1L
  )

  if (length(alias_unique_idx) > 0L) {

    mapping[
      alias_unique_idx,
      `:=`(
        gene_id = alias_entrez_ids,
        annotation_symbol = alias_official_symbols,
        annotation_ensembl = alias_ensembl_ids,
        annotation_gene_name = alias_gene_names,
        mapping_method = "rescued_unique_alias"
      )
    ]
  }


  # --------------------------------------------------------------------------
  # 4) Alias ambiguous — não resolver.
  # --------------------------------------------------------------------------

  alias_ambiguous_idx <- which(
    mapping$direct_n_entrez == 0L &
      !is.na(
        mapping$alias_n_entrez
      ) &
      mapping$alias_n_entrez > 1L
  )

  if (length(alias_ambiguous_idx) > 0L) {

    mapping[
      alias_ambiguous_idx,
      mapping_method := "ambiguous_alias"
    ]
  }


  # --------------------------------------------------------------------------
  # Flag final.
  # --------------------------------------------------------------------------

  mapping[
    ,
    resolved :=
      mapping_method %in%
        RESOLVED_METHODS
  ]


  # ============================================================================
  # 4.7 VALIDAR GENE IDs RESOLVIDOS
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "6. VALIDAÇÃO DOS GeneIDs RESOLVIDOS\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  # Usar um nome diferente da coluna lógica `resolved`.
  # Isso evita colisão de escopo dentro das expressões do data.table.
  resolved_mapping <- mapping[
    resolved == TRUE
  ]

  if (nrow(resolved_mapping) == 0L) {

    stop(
      "Nenhum gene foi resolvido para GeneID."
    )
  }

  invalid_gene_id <- resolved_mapping[
    is.na(gene_id) |
      !grepl(
        "^[0-9]+$",
        gene_id
      )
  ]

  if (nrow(invalid_gene_id) > 0L) {

    stop(
      "Existem ",
      nrow(invalid_gene_id),
      " linhas marcadas como resolvidas com GeneID inválido."
    )
  }

  n_resolved_rows <- nrow(
    resolved_mapping
  )

  n_unresolved_rows <- nrow(
    mapping
  ) -
    n_resolved_rows

  n_unique_geneids <- uniqueN(
    resolved_mapping$gene_id
  )

  pct_resolved_rows <- 100 *
    n_resolved_rows /
    nrow(
      mapping
    )

  cat(
    "Linhas resolvidas: ",
    n_resolved_rows,
    " / ",
    nrow(mapping),
    " (",
    round(
      pct_resolved_rows,
      2
    ),
    "%)\n",
    sep = ""
  )

  cat(
    "Linhas não resolvidas/ambíguas: ",
    n_unresolved_rows,
    "\n",
    sep = ""
  )

  cat(
    "GeneIDs únicos resolvidos: ",
    n_unique_geneids,
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.8 COLISÕES: MÚLTIPLAS LINHAS -> MESMO GeneID
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "7. COLISÕES APÓS HARMONIZAÇÃO\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  geneid_groups <- resolved_mapping[
    ,
    .(
      n_original_rows = .N,
      n_original_symbols = uniqueN(
        original_symbol
      ),
      original_symbols = collapse_values(
        original_symbol
      ),
      mapping_methods = collapse_values(
        mapping_method
      ),
      annotation_symbols = collapse_values(
        annotation_symbol
      )
    ),
    by = gene_id
  ][
    order(
      as.numeric(
        gene_id
      )
    )
  ]

  duplicated_geneids <- geneid_groups[
    n_original_rows > 1L
  ]

  n_duplicated_geneids <- nrow(
    duplicated_geneids
  )

  n_rows_in_duplicated_geneids <- if (
    n_duplicated_geneids > 0L
  ) {
    sum(
      duplicated_geneids$n_original_rows
    )
  } else {
    0L
  }

  cat(
    "GeneIDs recebendo >1 linha original: ",
    n_duplicated_geneids,
    "\n",
    sep = ""
  )

  cat(
    "Linhas envolvidas nessas colisões: ",
    n_rows_in_duplicated_geneids,
    "\n",
    sep = ""
  )

  # Adicionar o número de linhas originais associado a cada GeneID.
  #
  # Fazemos a atribuição para TODAS as linhas de `mapping`, e não apenas para
  # o subconjunto resolvido. Para linhas sem GeneID, match() retorna NA e,
  # consequentemente, geneid_n_original_rows também fica NA.
  #
  # Isso evita incompatibilidade de comprimento entre o vetor calculado para
  # 17.926 linhas e o subconjunto de linhas resolvidas.
  mapping[
    ,
    geneid_n_original_rows :=
      as.integer(
        geneid_groups$n_original_rows[
          match(
            gene_id,
            geneid_groups$gene_id
          )
        ]
      )
  ]

  mapping[
    ,
    geneid_collision :=
      resolved == TRUE &
        !is.na(
          geneid_n_original_rows
        ) &
        geneid_n_original_rows > 1L
  ]


  # ============================================================================
  # 4.9 RESUMOS POR MÉTODO
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "8. RESUMOS DE HARMONIZAÇÃO\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  method_summary <- mapping[
    ,
    .(
      n_rows = .N,
      n_unique_input_symbols = uniqueN(
        original_symbol
      ),
      n_unique_geneids = uniqueN(
        gene_id[
          !is.na(
            gene_id
          )
        ]
      )
    ),
    by = mapping_method
  ][
    order(
      -n_rows
    )
  ]

  print(
    method_summary
  )

  harmonization_summary <- data.table(
    metric = c(
      "input_gene_rows",
      "input_unique_symbols",
      "input_duplicated_symbol_rows",
      "direct_symbol_resolved_rows",
      "alias_rescued_rows",
      "ambiguous_direct_rows",
      "ambiguous_alias_rows",
      "unmapped_rows",
      "resolved_rows_total",
      "resolved_rows_percent",
      "unique_resolved_geneids",
      "geneids_with_multiple_input_rows",
      "rows_in_geneid_collisions"
    ),
    value = c(
      nrow(
        mapping
      ),
      uniqueN(
        mapping$original_symbol
      ),
      sum(
        duplicated(
          mapping$original_symbol
        )
      ),
      sum(
        mapping$mapping_method ==
          "direct_symbol"
      ),
      sum(
        mapping$mapping_method ==
          "rescued_unique_alias"
      ),
      sum(
        mapping$mapping_method ==
          "ambiguous_direct"
      ),
      sum(
        mapping$mapping_method ==
          "ambiguous_alias"
      ),
      sum(
        mapping$mapping_method ==
          "unmapped"
      ),
      n_resolved_rows,
      round(
        pct_resolved_rows,
        6
      ),
      n_unique_geneids,
      n_duplicated_geneids,
      n_rows_in_duplicated_geneids
    )
  )

  print(
    harmonization_summary
  )


  # ============================================================================
  # 4.10 TABELAS DE AUDITORIA
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "9. SALVANDO AUDITORIAS\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  ambiguous_direct <- mapping[
    mapping_method ==
      "ambiguous_direct"
  ]

  ambiguous_alias <- mapping[
    mapping_method ==
      "ambiguous_alias"
  ]

  unmapped <- mapping[
    mapping_method ==
      "unmapped"
  ]

  alias_rescued <- mapping[
    mapping_method ==
      "rescued_unique_alias"
  ]

  fwrite(
    method_summary,
    file.path(
      result_dir,
      "02_mapping_method_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  fwrite(
    harmonization_summary,
    file.path(
      result_dir,
      "02_harmonization_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  fwrite(
    ambiguous_direct,
    file.path(
      result_dir,
      "02_ambiguous_direct_symbol.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  fwrite(
    ambiguous_alias,
    file.path(
      result_dir,
      "02_ambiguous_alias.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  fwrite(
    unmapped,
    file.path(
      result_dir,
      "02_unmapped_symbols.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  fwrite(
    alias_rescued,
    file.path(
      result_dir,
      "02_alias_rescued.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  fwrite(
    duplicated_geneids,
    file.path(
      result_dir,
      "02_duplicated_geneids.tsv"
    ),
    sep = "\t",
    na = "NA"
  )


  # ============================================================================
  # 4.11 OUTPUTS PRINCIPAIS PARA O SCRIPT 03
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "10. OUTPUTS PARA 03_build_reference_Mathys2019.R\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  # Mapa completo por linha da matriz original.
  fwrite(
    mapping,
    mapping_output_file,
    sep = "\t",
    na = "NA"
  )

  # Mapa de colapso:
  # uma linha por GeneID resolvido; o script 03 usará matrix_row_index para somar
  # counts quando n_original_rows > 1.
  collapse_map <- resolved_mapping[
    ,
    .(
      gene_id,
      matrix_row_index,
      original_symbol,
      annotation_symbol,
      mapping_method
    )
  ][
    order(
      as.numeric(
        gene_id
      ),
      matrix_row_index
    )
  ]

  fwrite(
    collapse_map,
    collapse_map_file,
    sep = "\t",
    na = "NA"
  )

  # Vetor com as linhas resolvidas da matriz original.
  # Mantém índice 1-based de R.
  resolved_rows <- resolved_mapping$matrix_row_index

  saveRDS(
    resolved_rows,
    resolved_rows_file,
    compress = FALSE
  )

  cat(
    "Mapa completo:\n  ",
    mapping_output_file,
    "\n",
    sep = ""
  )

  cat(
    "Mapa de colapso:\n  ",
    collapse_map_file,
    "\n",
    sep = ""
  )

  cat(
    "Índices de linhas resolvidas:\n  ",
    resolved_rows_file,
    "\n",
    sep = ""
  )


  # ============================================================================
  # 4.12 CHECKS FINAIS
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "11. CHECKS FINAIS\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  checks <- rbindlist(
    list(

      make_check(
        "input_gene_rows",
        nrow(mapping),
        EXPECTED_MATHYS_GENE_ROWS,
        if (
          nrow(mapping) ==
            EXPECTED_MATHYS_GENE_ROWS
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "matrix_row_index_preserved",
        identical(
          mapping$matrix_row_index,
          seq_len(
            nrow(mapping)
          )
        ),
        TRUE,
        if (
          identical(
            mapping$matrix_row_index,
            seq_len(
              nrow(mapping)
            )
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "resolved_rows_have_gene_id",
        all(
          !is.na(
            mapping$gene_id[
              mapping$resolved
            ]
          )
        ),
        TRUE,
        if (
          all(
            !is.na(
              mapping$gene_id[
                mapping$resolved
              ]
            )
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "resolved_gene_ids_numeric",
        all(
          grepl(
            "^[0-9]+$",
            mapping$gene_id[
              mapping$resolved
            ]
          )
        ),
        TRUE,
        if (
          all(
            grepl(
              "^[0-9]+$",
              mapping$gene_id[
                mapping$resolved
              ]
            )
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      ),

      make_check(
        "ambiguous_mappings_not_resolved",
        all(
          !mapping$resolved[
            mapping$mapping_method %in%
              c(
                "ambiguous_direct",
                "ambiguous_alias"
              )
          ]
        ),
        TRUE,
        if (
          all(
            !mapping$resolved[
              mapping$mapping_method %in%
                c(
                  "ambiguous_direct",
                  "ambiguous_alias"
                )
            ]
          )
        ) {
          "PASS"
        } else {
          "FAIL"
        },
        "Nenhum mapping ambíguo pode receber GeneID automaticamente."
      ),

      make_check(
        "mapping_rows_equal_input_rows",
        nrow(mapping) ==
          nrow(genes),
        TRUE,
        if (
          nrow(mapping) ==
            nrow(genes)
        ) {
          "PASS"
        } else {
          "FAIL"
        }
      )
    ),
    fill = TRUE
  )

  fwrite(
    checks,
    file.path(
      result_dir,
      "02_harmonization_status.tsv"
    ),
    sep = "\t",
    na = "NA"
  )

  print(
    checks
  )

  if (
    any(
      checks$status == "FAIL",
      na.rm = TRUE
    )
  ) {

    stop(
      "Há check(s) FAIL após a harmonização."
    )
  }


  # ============================================================================
  # 4.13 SESSION INFO
  # ============================================================================

  sink(
    session_file
  )

  print(
    sessionInfo()
  )

  sink()


  # ============================================================================
  # 4.14 CONCLUSÃO
  # ============================================================================

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "12. CONCLUSÃO\n"
  )

  cat(
    "------------------------------------------------------------\n"
  )

  cat(
    "Genes de entrada: ",
    nrow(mapping),
    "\n",
    sep = ""
  )

  cat(
    "Linhas resolvidas para GeneID: ",
    n_resolved_rows,
    " (",
    round(
      pct_resolved_rows,
      2
    ),
    "%)\n",
    sep = ""
  )

  cat(
    "GeneIDs únicos: ",
    n_unique_geneids,
    "\n",
    sep = ""
  )

  cat(
    "Resgatados por alias único: ",
    sum(
      mapping$mapping_method ==
        "rescued_unique_alias"
    ),
    "\n",
    sep = ""
  )

  cat(
    "Mappings ambíguos não aceitos: ",
    sum(
      mapping$mapping_method %in%
        c(
          "ambiguous_direct",
          "ambiguous_alias"
        )
    ),
    "\n",
    sep = ""
  )

  cat(
    "Símbolos não resolvidos: ",
    sum(
      mapping$mapping_method ==
        "unmapped"
    ),
    "\n",
    sep = ""
  )

  cat(
    "GeneIDs com múltiplas linhas originais: ",
    n_duplicated_geneids,
    "\n",
    sep = ""
  )

  cat(
    "\nNenhuma matriz de expressão foi modificada.\n"
  )

  cat(
    "Nenhuma linha duplicada foi colapsada nesta etapa.\n"
  )

  cat(
    "Nenhum GeneID ambíguo foi escolhido automaticamente.\n"
  )

  cat(
    "\nPróxima etapa:\n"
  )

  cat(
    "  03_build_reference_Mathys2019.R\n"
  )

  invisible(TRUE)
}


# ==============================================================================
# 5. EXECUÇÃO COM LOG
# ==============================================================================

log_connection <- file(
  log_file,
  open = "wt"
)

sink(
  log_connection,
  split = TRUE
)

execution_error <- NULL

tryCatch(
  {
    run_harmonization()
  },
  error = function(e) {

    execution_error <<- e

    cat(
      "\n============================================================\n"
    )

    cat(
      "ERRO NA HARMONIZAÇÃO GÊNICA\n"
    )

    cat(
      "============================================================\n"
    )

    cat(
      conditionMessage(
        e
      ),
      "\n"
    )
  }
)

sink()

close(
  log_connection
)

if (!is.null(execution_error)) {

  stop(
    conditionMessage(
      execution_error
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 6. MENSAGEM FINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "02_harmonize_genes_Mathys2019.R CONCLUÍDO\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Log:\n  ",
  log_file,
  "\n\n",
  sep = ""
)

cat(
  "Resumo:\n  ",
  file.path(
    result_dir,
    "02_harmonization_summary.tsv"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Status:\n  ",
  file.path(
    result_dir,
    "02_harmonization_status.tsv"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Mapa final:\n  ",
  mapping_output_file,
  "\n\n",
  sep = ""
)

cat(
  "Colisões por GeneID:\n  ",
  file.path(
    result_dir,
    "02_duplicated_geneids.tsv"
  ),
  "\n",
  sep = ""
)
