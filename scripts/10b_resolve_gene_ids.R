#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 08b_resolve_gene_ids.R
#
# Segunda camada de harmonização gênica:
#
# 1. preserva matches diretos por símbolo;
# 2. tenta resgatar genes não resolvidos usando:
#       Ensembl gene_id -> annotation -> Symbol;
# 3. verifica se o símbolo anotado existe na referência Mathys;
# 4. audita ambiguidades e duplicações;
# 5. NÃO modifica matrizes de expressão;
# 6. NÃO colapsa genes duplicados nesta etapa.
# ============================================================


# ============================================================
# 1. CAMINHOS
# ============================================================

mathys_file <- paste0(
  "data/reference/mathys2019/prepared/",
  "mathys2019_genes.tsv"
)

annotation_file <- paste0(
  "data/raw/",
  "Human.GRCh38.p13.annot.tsv.gz"
)

bulk_files <- c(
  AD = "data/processed/GSE53697_expression_prepared.csv",
  HD = "data/processed/GSE64810_expression_prepared.csv",
  PD = "data/processed/GSE68719_expression_prepared.csv"
)

out_dir <- paste0(
  "data/reference/mathys2019/",
  "harmonization/resolved"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. VERIFICAR INPUTS
# ============================================================

cat("\n============================================\n")
cat("1. VERIFICAÇÃO DOS INPUTS\n")
cat("============================================\n\n")

required_files <- c(
  mathys_file,
  annotation_file,
  bulk_files
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  cat("Arquivos ausentes:\n")
  print(missing_files)

  stop("Inputs ausentes.")
}

cat("Todos os inputs foram encontrados.\n")


# ============================================================
# 3. REFERÊNCIA MATHYS
# ============================================================

cat("\n============================================\n")
cat("2. REFERÊNCIA MATHYS\n")
cat("============================================\n\n")

mathys <- fread(mathys_file)

if (!"gene_symbol" %in% names(mathys)) {
  stop("gene_symbol ausente na referência Mathys.")
}

mathys_symbols <- unique(
  mathys$gene_symbol[
    !is.na(mathys$gene_symbol) &
    trimws(mathys$gene_symbol) != ""
  ]
)

cat(
  "Genes Mathys válidos:",
  length(mathys_symbols),
  "\n"
)


# ============================================================
# 4. CARREGAR ANOTAÇÃO GRCh38.p13
# ============================================================

cat("\n============================================\n")
cat("3. ANOTAÇÃO GRCh38.p13\n")
cat("============================================\n\n")

annot <- fread(annotation_file)

cat(
  "Linhas da anotação:",
  nrow(annot),
  "\n"
)

cat(
  "Colunas disponíveis:\n"
)

print(names(annot))


required_annot_cols <- c(
  "Symbol",
  "EnsemblGeneID"
)

missing_annot_cols <- setdiff(
  required_annot_cols,
  names(annot)
)

if (length(missing_annot_cols) > 0) {

  stop(
    paste(
      "Colunas ausentes na anotação:",
      paste(
        missing_annot_cols,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 5. LIMPAR MAPEAMENTO ENSEMBL -> SYMBOL
# ============================================================

cat("\n============================================\n")
cat("4. PREPARANDO MAPA ENSEMBL -> SYMBOL\n")
cat("============================================\n\n")

annot_map <- annot[
  !is.na(EnsemblGeneID) &
  trimws(EnsemblGeneID) != "" &
  !is.na(Symbol) &
  trimws(Symbol) != "",
  .(
    gene_id = trimws(EnsemblGeneID),
    annotation_symbol = trimws(Symbol)
  )
]

# Remover versão do Ensembl caso exista.
annot_map[
  ,
  gene_id := sub(
    "\\.[0-9]+$",
    "",
    gene_id
  )
]

annot_map <- unique(annot_map)

cat(
  "Pares Ensembl-Symbol únicos:",
  nrow(annot_map),
  "\n"
)


# ============================================================
# 6. AUDITAR AMBIGUIDADES DA ANOTAÇÃO
# ============================================================

ensembl_ambiguity <- annot_map[
  ,
  .(
    n_symbols = uniqueN(annotation_symbol),
    symbols = paste(
      sort(unique(annotation_symbol)),
      collapse = ";"
    )
  ),
  by = gene_id
][
  n_symbols > 1
]

cat(
  "Ensembl IDs associados a >1 símbolo:",
  nrow(ensembl_ambiguity),
  "\n"
)

ambiguous_ensembl <- ensembl_ambiguity$gene_id

# Para o resgate automático, usar apenas relações 1:1.
annot_map_unique <- annot_map[
  !gene_id %in% ambiguous_ensembl
]

if (
  anyDuplicated(
    annot_map_unique$gene_id
  ) > 0
) {

  stop(
    "Ainda existem Ensembl IDs duplicados após remoção das ambiguidades."
  )
}

cat(
  "Mapeamentos Ensembl 1:1 utilizáveis:",
  nrow(annot_map_unique),
  "\n"
)


# ============================================================
# 7. FUNÇÃO DE HARMONIZAÇÃO
# ============================================================

cat("\n============================================\n")
cat("5. HARMONIZAÇÃO POR DATASET\n")
cat("============================================\n\n")

resolve_dataset <- function(
  file,
  dataset,
  mathys_symbols,
  annot_map_unique,
  ambiguous_ensembl
) {

  cat(
    "\n---",
    dataset,
    "---\n"
  )

  x <- fread(file)

  required_cols <- c(
    "gene_id",
    "original_id",
    "symbol"
  )

  missing_cols <- setdiff(
    required_cols,
    names(x)
  )

  if (length(missing_cols) > 0) {

    stop(
      paste(
        dataset,
        "sem colunas:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  # Trabalhar apenas com anotação.
  x <- x[
    ,
    .(
      gene_id,
      original_id,
      original_symbol = symbol
    )
  ]

  # Padronizar Ensembl removendo versão.
  x[
    ,
    gene_id_clean := sub(
      "\\.[0-9]+$",
      "",
      trimws(gene_id)
    )
  ]

  # ----------------------------------------------------------
  # Match direto por símbolo
  # ----------------------------------------------------------

  x[
    ,
    direct_match :=
      !is.na(original_symbol) &
      trimws(original_symbol) != "" &
      original_symbol %in% mathys_symbols
  ]


  # ----------------------------------------------------------
  # Buscar símbolo pela anotação Ensembl
  # ----------------------------------------------------------

  idx <- match(
    x$gene_id_clean,
    annot_map_unique$gene_id
  )

  x[
    ,
    annotation_symbol :=
      annot_map_unique$annotation_symbol[idx]
  ]


  # ----------------------------------------------------------
  # Verificar se símbolo recuperado existe em Mathys
  # ----------------------------------------------------------

  x[
    ,
    ensembl_mathys_match :=
      !is.na(annotation_symbol) &
      annotation_symbol %in% mathys_symbols
  ]


  # ----------------------------------------------------------
  # Identificar Ensembl ambíguo
  # ----------------------------------------------------------

  x[
    ,
    ensembl_ambiguous :=
      gene_id_clean %in%
      ambiguous_ensembl
  ]


  # ----------------------------------------------------------
  # Definir símbolo harmonizado
  #
  # Prioridade:
  # 1. símbolo original com match direto
  # 2. símbolo recuperado por Ensembl com match Mathys
  # 3. NA
  # ----------------------------------------------------------

  x[
    ,
    harmonized_symbol := fifelse(
      direct_match,
      original_symbol,
      fifelse(
        ensembl_mathys_match,
        annotation_symbol,
        NA_character_
      )
    )
  ]


  # ----------------------------------------------------------
  # Método de resolução
  # ----------------------------------------------------------

  x[
    ,
    mapping_method := fifelse(
      direct_match,
      "direct_symbol",
      fifelse(
        ensembl_mathys_match,
        "rescued_by_ensembl",
        fifelse(
          ensembl_ambiguous,
          "ambiguous_ensembl",
          "unresolved"
        )
      )
    )
  ]


  # ----------------------------------------------------------
  # Estatísticas
  # ----------------------------------------------------------

  n_direct <- sum(
    x$mapping_method ==
    "direct_symbol"
  )

  n_rescued <- sum(
    x$mapping_method ==
    "rescued_by_ensembl"
  )

  n_ambiguous <- sum(
    x$mapping_method ==
    "ambiguous_ensembl"
  )

  n_unresolved <- sum(
    x$mapping_method ==
    "unresolved"
  )

  cat(
    "Linhas:",
    nrow(x),
    "\n"
  )

  cat(
    "Direct symbol:",
    n_direct,
    "\n"
  )

  cat(
    "Resgatados por Ensembl:",
    n_rescued,
    "\n"
  )

  cat(
    "Ensembl ambíguo:",
    n_ambiguous,
    "\n"
  )

  cat(
    "Não resolvidos:",
    n_unresolved,
    "\n"
  )

  cat(
    "Total compatível Mathys:",
    n_direct + n_rescued,
    "\n"
  )


  # ----------------------------------------------------------
  # Símbolos harmonizados duplicados
  # ----------------------------------------------------------

  resolved <- x[
    !is.na(harmonized_symbol)
  ]

  duplicate_symbols <- resolved[
    ,
    .N,
    by = harmonized_symbol
  ][
    N > 1
  ]

  cat(
    "Símbolos harmonizados duplicados:",
    nrow(duplicate_symbols),
    "\n"
  )


  # ----------------------------------------------------------
  # Resumo
  # ----------------------------------------------------------

  summary <- data.table(
    dataset = dataset,
    n_bulk_rows = nrow(x),
    direct_symbol = n_direct,
    rescued_by_ensembl = n_rescued,
    ambiguous_ensembl = n_ambiguous,
    unresolved = n_unresolved,
    compatible_rows =
      n_direct + n_rescued,
    unique_harmonized_symbols =
      uniqueN(
        x$harmonized_symbol[
          !is.na(x$harmonized_symbol)
        ]
      ),
    duplicated_harmonized_symbols =
      nrow(duplicate_symbols)
  )

  list(
    mapping = x,
    summary = summary,
    duplicates = duplicate_symbols
  )
}


# ============================================================
# 8. EXECUTAR PARA AD / HD / PD
# ============================================================

resolved <- list()

for (dataset in names(bulk_files)) {

  resolved[[dataset]] <- resolve_dataset(
    file = bulk_files[[dataset]],
    dataset = dataset,
    mathys_symbols = mathys_symbols,
    annot_map_unique = annot_map_unique,
    ambiguous_ensembl = ambiguous_ensembl
  )
}


# ============================================================
# 9. RESUMO GLOBAL
# ============================================================

cat("\n============================================\n")
cat("6. RESUMO GLOBAL\n")
cat("============================================\n\n")

resolution_summary <- rbindlist(
  lapply(
    resolved,
    function(z) z$summary
  )
)

print(resolution_summary)


# ============================================================
# 10. UNIVERSOS FINAIS APÓS RESGATE
# ============================================================

cat("\n============================================\n")
cat("7. UNIVERSOS APÓS RESGATE\n")
cat("============================================\n\n")

final_symbols <- lapply(
  resolved,
  function(z) {

    unique(
      z$mapping$harmonized_symbol[
        !is.na(
          z$mapping$harmonized_symbol
        )
      ]
    )
  }
)

for (dataset in names(final_symbols)) {

  cat(
    "Mathys ∩",
    dataset,
    ":",
    length(final_symbols[[dataset]]),
    "símbolos únicos\n"
  )
}


common_all <- Reduce(
  intersect,
  c(
    list(mathys_symbols),
    final_symbols
  )
)

cat(
  "\nMathys ∩ AD ∩ HD ∩ PD:",
  length(common_all),
  "genes\n"
)

cat(
  "Percentual da referência Mathys:",
  round(
    100 *
    length(common_all) /
    length(mathys_symbols),
    2
  ),
  "%\n"
)


# ============================================================
# 11. STATUS GLOBAL DOS GENES MATHYS
# ============================================================

cat("\n============================================\n")
cat("8. STATUS GLOBAL MATHYS\n")
cat("============================================\n\n")

mathys_status <- data.table(
  mathys_gene = mathys_symbols
)

for (dataset in names(final_symbols)) {

  mathys_status[
    ,
    (paste0("match_", dataset)) :=
      mathys_gene %in%
      final_symbols[[dataset]]
  ]
}

mathys_status[
  ,
  match_all :=
    match_AD &
    match_HD &
    match_PD
]

mathys_status[
  ,
  n_datasets :=
    as.integer(match_AD) +
    as.integer(match_HD) +
    as.integer(match_PD)
]

status_summary <- mathys_status[
  ,
  .N,
  by = n_datasets
][order(n_datasets)]

print(status_summary)


# ============================================================
# 12. SALVAR OUTPUTS
# ============================================================

cat("\n============================================\n")
cat("9. SALVANDO RESULTADOS\n")
cat("============================================\n\n")

fwrite(
  resolution_summary,
  file.path(
    out_dir,
    "resolution_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  ensembl_ambiguity,
  file.path(
    out_dir,
    "annotation_ensembl_ambiguities.tsv"
  ),
  sep = "\t"
)

for (dataset in names(resolved)) {

  fwrite(
    resolved[[dataset]]$mapping,
    file.path(
      out_dir,
      paste0(
        dataset,
        "_resolved_gene_mapping.tsv"
      )
    ),
    sep = "\t"
  )

  fwrite(
    resolved[[dataset]]$duplicates,
    file.path(
      out_dir,
      paste0(
        dataset,
        "_duplicated_harmonized_symbols.tsv"
      )
    ),
    sep = "\t"
  )

  fwrite(
    data.table(
      gene_symbol =
        final_symbols[[dataset]]
    ),
    file.path(
      out_dir,
      paste0(
        dataset,
        "_mathys_compatible_symbols.tsv"
      )
    ),
    sep = "\t"
  )
}

fwrite(
  data.table(
    gene_symbol = common_all
  ),
  file.path(
    out_dir,
    "common_mathys_AD_HD_PD_symbols.tsv"
  ),
  sep = "\t"
)

fwrite(
  mathys_status,
  file.path(
    out_dir,
    "mathys_final_match_status.tsv"
  ),
  sep = "\t"
)

fwrite(
  status_summary,
  file.path(
    out_dir,
    "mathys_final_match_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 13. VALIDAÇÃO FINAL
# ============================================================

cat("\n============================================\n")
cat("10. VALIDAÇÃO FINAL\n")
cat("============================================\n\n")

cat(
  "Mathys:",
  length(mathys_symbols),
  "genes\n"
)

for (dataset in names(final_symbols)) {

  cat(
    "Mathys ∩",
    dataset,
    ":",
    length(final_symbols[[dataset]]),
    "\n"
  )
}

cat(
  "Mathys ∩ AD ∩ HD ∩ PD:",
  length(common_all),
  "\n"
)

cat(
  "\nOutputs:\n",
  out_dir,
  "\n",
  sep = ""
)

cat(
  "\nNenhuma matriz de expressão foi modificada.\n"
)

cat(
  "Nenhum gene duplicado foi colapsado.\n"
)

cat(
  "Nenhum alias histórico foi resolvido nesta etapa.\n"
)

cat(
  "Resgate por Ensembl concluído.\n"
)
