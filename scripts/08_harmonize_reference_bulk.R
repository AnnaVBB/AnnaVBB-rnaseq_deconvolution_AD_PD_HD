#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 08_harmonize_reference_bulk.R
#
# Auditoria e harmonização dos identificadores gênicos entre:
# - referência snRNA-seq Mathys et al. 2019
# - GSE53697 (AD)
# - GSE64810 (HD)
# - GSE68719 (PD)
#
# Nesta etapa:
# - nenhuma expressão é transformada;
# - nenhum gene é removido das matrizes originais;
# - nenhum símbolo é atualizado por alias;
# - nenhum universo é imposto aos algoritmos;
# - apenas matches diretos por gene symbol são avaliados.
# ============================================================


# ============================================================
# 1. CAMINHOS
# ============================================================

mathys_file <- paste0(
  "data/reference/mathys2019/prepared/",
  "mathys2019_genes.tsv"
)

bulk_files <- c(
  AD = "data/processed/GSE53697_expression_prepared.csv",
  HD = "data/processed/GSE64810_expression_prepared.csv",
  PD = "data/processed/GSE68719_expression_prepared.csv"
)

out_dir <- "data/reference/mathys2019/harmonization"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. VERIFICAR ARQUIVOS
# ============================================================

cat("\n============================================\n")
cat("1. VERIFICAÇÃO DOS ARQUIVOS\n")
cat("============================================\n\n")

required_files <- c(
  mathys_file,
  bulk_files
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  cat("Arquivos ausentes:\n")
  print(missing_files)

  stop(
    "Harmonização interrompida: arquivos ausentes."
  )
}

cat("Todos os arquivos foram encontrados.\n")


# ============================================================
# 3. CARREGAR GENES MATHYS
# ============================================================

cat("\n============================================\n")
cat("2. REFERÊNCIA MATHYS\n")
cat("============================================\n\n")

mathys <- fread(mathys_file)

if (!"gene_symbol" %in% names(mathys)) {
  stop("Coluna gene_symbol ausente na referência Mathys.")
}

mathys[, mathys_gene_original := gene_symbol]

cat(
  "Genes Mathys:",
  nrow(mathys),
  "\n"
)

cat(
  "Símbolos únicos:",
  uniqueN(mathys$gene_symbol),
  "\n"
)

cat(
  "Símbolos ausentes:",
  sum(
    is.na(mathys$gene_symbol) |
    mathys$gene_symbol == ""
  ),
  "\n"
)

cat(
  "Símbolos duplicados:",
  sum(duplicated(mathys$gene_symbol)),
  "\n"
)


# ============================================================
# 4. CARREGAR ANOTAÇÕES BULK
# ============================================================

cat("\n============================================\n")
cat("3. DATASETS BULK\n")
cat("============================================\n\n")

bulk <- list()

for (dataset in names(bulk_files)) {

  cat("Carregando", dataset, "...\n")

  x <- fread(
    bulk_files[[dataset]]
  )

  required_columns <- c(
    "gene_id",
    "original_id",
    "symbol"
  )

  missing_columns <- setdiff(
    required_columns,
    names(x)
  )

  if (length(missing_columns) > 0) {

    stop(
      paste(
        dataset,
        "não possui as colunas:",
        paste(missing_columns, collapse = ", ")
      )
    )
  }

  # Nesta etapa precisamos apenas das anotações.
  # Não mantemos as dezenas de colunas de expressão em memória.

  bulk[[dataset]] <- x[
    ,
    .(
      gene_id,
      original_id,
      symbol
    )
  ]

  rm(x)
  gc()

  cat(
    dataset,
    ":",
    nrow(bulk[[dataset]]),
    "genes\n"
  )
}


# ============================================================
# 5. AUDITORIA DOS IDs BULK
# ============================================================

cat("\n============================================\n")
cat("4. AUDITORIA DOS IDENTIFICADORES BULK\n")
cat("============================================\n\n")

audit_bulk <- rbindlist(
  lapply(
    names(bulk),
    function(dataset) {

      x <- bulk[[dataset]]

      valid_symbol <- (
        !is.na(x$symbol) &
        trimws(x$symbol) != ""
      )

      data.table(
        dataset = dataset,

        n_rows = nrow(x),

        n_gene_id_present = sum(
          !is.na(x$gene_id) &
          trimws(x$gene_id) != ""
        ),

        n_gene_id_unique = uniqueN(
          x$gene_id[
            !is.na(x$gene_id) &
            trimws(x$gene_id) != ""
          ]
        ),

        n_symbol_present = sum(valid_symbol),

        n_symbol_missing = sum(!valid_symbol),

        n_symbol_unique = uniqueN(
          x$symbol[valid_symbol]
        ),

        n_symbol_duplicated_rows = sum(
          duplicated(
            x$symbol[valid_symbol]
          )
        )
      )
    }
  ),
  fill = TRUE
)

print(audit_bulk)


# ============================================================
# 6. PREPARAR CONJUNTOS DE SÍMBOLOS
# ============================================================

cat("\n============================================\n")
cat("5. PREPARAÇÃO DOS SÍMBOLOS\n")
cat("============================================\n\n")

valid_symbols <- function(x) {

  unique(
    x[
      !is.na(x) &
      trimws(x) != ""
    ]
  )
}

mathys_symbols <- valid_symbols(
  mathys$gene_symbol
)

bulk_symbols <- lapply(
  bulk,
  function(x) {
    valid_symbols(x$symbol)
  }
)

cat(
  "Mathys:",
  length(mathys_symbols),
  "símbolos válidos\n"
)

for (dataset in names(bulk_symbols)) {

  cat(
    dataset,
    ":",
    length(bulk_symbols[[dataset]]),
    "símbolos válidos\n"
  )
}


# ============================================================
# 7. INTERSEÇÕES MATHYS × BULK
# ============================================================

cat("\n============================================\n")
cat("6. INTERSEÇÕES MATHYS x BULK\n")
cat("============================================\n\n")

direct_matches <- lapply(
  bulk_symbols,
  function(x) {
    intersect(
      mathys_symbols,
      x
    )
  }
)

for (dataset in names(direct_matches)) {

  n_match <- length(
    direct_matches[[dataset]]
  )

  pct_mathys <- 100 *
    n_match /
    length(mathys_symbols)

  pct_bulk <- 100 *
    n_match /
    length(bulk_symbols[[dataset]])

  cat(
    dataset,
    ":\n"
  )

  cat(
    "  Match direto:",
    n_match,
    "\n"
  )

  cat(
    "  % da referência Mathys:",
    round(pct_mathys, 2),
    "%\n"
  )

  cat(
    "  % do universo bulk:",
    round(pct_bulk, 2),
    "%\n\n"
  )
}


# ============================================================
# 8. UNIVERSO COMUM AOS 4 CONJUNTOS
# ============================================================

cat("\n============================================\n")
cat("7. UNIVERSO COMUM\n")
cat("============================================\n\n")

common_all <- Reduce(
  intersect,
  c(
    list(mathys_symbols),
    bulk_symbols
  )
)

cat(
  "Mathys ∩ AD ∩ HD ∩ PD:",
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
# 9. GENES MATHYS NÃO ENCONTRADOS
# ============================================================

cat("\n============================================\n")
cat("8. GENES MATHYS NÃO ENCONTRADOS\n")
cat("============================================\n\n")

unmatched_mathys <- list()

for (dataset in names(bulk_symbols)) {

  unmatched <- setdiff(
    mathys_symbols,
    bulk_symbols[[dataset]]
  )

  unmatched_mathys[[dataset]] <- unmatched

  cat(
    dataset,
    ":",
    length(unmatched),
    "genes Mathys sem match direto\n"
  )
}


# ============================================================
# 10. STATUS DE CADA GENE MATHYS
# ============================================================

cat("\n============================================\n")
cat("9. TABELA DE STATUS MATHYS\n")
cat("============================================\n\n")

mathys_status <- data.table(
  mathys_gene_original = mathys_symbols
)

mathys_status[
  ,
  match_AD :=
    mathys_gene_original %in%
    bulk_symbols$AD
]

mathys_status[
  ,
  match_HD :=
    mathys_gene_original %in%
    bulk_symbols$HD
]

mathys_status[
  ,
  match_PD :=
    mathys_gene_original %in%
    bulk_symbols$PD
]

mathys_status[
  ,
  match_all :=
    match_AD &
    match_HD &
    match_PD
]

mathys_status[
  ,
  n_bulk_matches :=
    as.integer(match_AD) +
    as.integer(match_HD) +
    as.integer(match_PD)
]

mathys_status[
  ,
  match_pattern := fifelse(
    match_all,
    "AD_HD_PD",
    fifelse(
      match_AD & match_HD,
      "AD_HD",
      fifelse(
        match_AD & match_PD,
        "AD_PD",
        fifelse(
          match_HD & match_PD,
          "HD_PD",
          fifelse(
            match_AD,
            "AD_only",
            fifelse(
              match_HD,
              "HD_only",
              fifelse(
                match_PD,
                "PD_only",
                "None"
              )
            )
          )
        )
      )
    )
  )
]

pattern_summary <- mathys_status[
  ,
  .N,
  by = match_pattern
][order(-N)]

cat("Padrões de compatibilidade:\n")

print(pattern_summary)


# ============================================================
# 11. TABELAS DE MAPEAMENTO POR DATASET
# ============================================================

cat("\n============================================\n")
cat("10. TABELAS DE MAPEAMENTO\n")
cat("============================================\n\n")

mapping_tables <- list()

for (dataset in names(bulk)) {

  x <- copy(
    bulk[[dataset]]
  )

  x[
    ,
    direct_mathys_match :=
      !is.na(symbol) &
      symbol %in% mathys_symbols
  ]

  mapping_tables[[dataset]] <- x

  cat(
    dataset,
    "- linhas com match Mathys:",
    sum(x$direct_mathys_match),
    "\n"
  )
}


# ============================================================
# 12. RESUMO DAS INTERSEÇÕES
# ============================================================

intersection_summary <- rbindlist(
  lapply(
    names(bulk_symbols),
    function(dataset) {

      matches <- direct_matches[[dataset]]

      data.table(
        dataset = dataset,

        mathys_genes =
          length(mathys_symbols),

        bulk_symbols =
          length(bulk_symbols[[dataset]]),

        direct_matches =
          length(matches),

        mathys_unmatched =
          length(
            unmatched_mathys[[dataset]]
          ),

        pct_mathys_matched =
          100 *
          length(matches) /
          length(mathys_symbols),

        pct_bulk_matched =
          100 *
          length(matches) /
          length(bulk_symbols[[dataset]])
      )
    }
  )
)


# ============================================================
# 13. SALVAR RESULTADOS
# ============================================================

cat("\n============================================\n")
cat("11. SALVANDO AUDITORIA\n")
cat("============================================\n\n")

fwrite(
  audit_bulk,
  file.path(
    out_dir,
    "bulk_identifier_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  intersection_summary,
  file.path(
    out_dir,
    "mathys_bulk_intersection_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  mathys_status,
  file.path(
    out_dir,
    "mathys_gene_match_status.tsv"
  ),
  sep = "\t"
)

fwrite(
  pattern_summary,
  file.path(
    out_dir,
    "mathys_match_pattern_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  data.table(
    gene_symbol = common_all
  ),
  file.path(
    out_dir,
    "mathys_AD_HD_PD_common_symbols.tsv"
  ),
  sep = "\t"
)

for (dataset in names(mapping_tables)) {

  fwrite(
    mapping_tables[[dataset]],
    file.path(
      out_dir,
      paste0(
        dataset,
        "_mathys_mapping.tsv"
      )
    ),
    sep = "\t"
  )

  fwrite(
    data.table(
      mathys_gene_original =
        unmatched_mathys[[dataset]]
    ),
    file.path(
      out_dir,
      paste0(
        "mathys_unmatched_",
        dataset,
        ".tsv"
      )
    ),
    sep = "\t"
  )
}


# ============================================================
# 14. VALIDAÇÃO FINAL
# ============================================================

cat("\n============================================\n")
cat("12. VALIDAÇÃO FINAL\n")
cat("============================================\n\n")

cat(
  "Genes Mathys:",
  length(mathys_symbols),
  "\n"
)

for (dataset in names(direct_matches)) {

  cat(
    "Mathys ∩",
    dataset,
    ":",
    length(direct_matches[[dataset]]),
    "\n"
  )
}

cat(
  "Mathys ∩ AD ∩ HD ∩ PD:",
  length(common_all),
  "\n"
)

cat(
  "\nResultados salvos em:\n",
  out_dir,
  "\n",
  sep = ""
)

cat(
  "\nNenhuma matriz de expressão foi modificada.\n"
)

cat(
  "Nenhum alias foi resolvido nesta etapa.\n"
)

cat(
  "Harmonização direta por símbolos concluída.\n"
)
