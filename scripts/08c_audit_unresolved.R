#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 08c_audit_unresolved.R
#
# Auditoria final dos genes não resolvidos após o 08b.
# Nenhuma matriz é modificada.
# ============================================================

in_dir <- "data/reference/mathys2019/harmonization/resolved"

out_dir <- file.path(
  in_dir,
  "audit"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

datasets <- c("AD", "HD", "PD")

all_summary <- list()


for (dataset in datasets) {

  cat("\n============================================\n")
  cat(dataset, "\n")
  cat("============================================\n\n")

  file <- file.path(
    in_dir,
    paste0(
      dataset,
      "_resolved_gene_mapping.tsv"
    )
  )

  x <- fread(file)


  # ----------------------------------------------------------
  # Separar categorias
  # ----------------------------------------------------------

  unresolved <- x[
    mapping_method == "unresolved"
  ]

  unresolved[
    ,
    has_annotation_symbol :=
      !is.na(annotation_symbol) &
      trimws(annotation_symbol) != ""
  ]

  unresolved[
    ,
    annotation_not_in_mathys :=
      has_annotation_symbol &
      !ensembl_mathys_match
  ]

  no_annotation <- unresolved[
    has_annotation_symbol == FALSE
  ]

  annotated_not_mathys <- unresolved[
    annotation_not_in_mathys == TRUE
  ]


  # ----------------------------------------------------------
  # Duplicações entre símbolos harmonizados
  # ----------------------------------------------------------

  resolved <- x[
  !is.na(harmonized_symbol) &
  trimws(harmonized_symbol) != "" &
  harmonized_symbol != "NA"
  ]

  duplicates <- resolved[
    ,
    .N,
    by = harmonized_symbol
  ][
    N > 1
  ][
    order(-N, harmonized_symbol)
  ]


  # ----------------------------------------------------------
  # Estatísticas
  # ----------------------------------------------------------

  summary <- data.table(
    dataset = dataset,

    total_rows = nrow(x),

    direct_symbol = sum(
      x$mapping_method ==
        "direct_symbol"
    ),

    rescued_by_ensembl = sum(
      x$mapping_method ==
        "rescued_by_ensembl"
    ),

    ambiguous_ensembl = sum(
      x$mapping_method ==
        "ambiguous_ensembl"
    ),

    unresolved = nrow(unresolved),

    unresolved_with_annotation_symbol =
      sum(
        unresolved$has_annotation_symbol
      ),

    unresolved_annotation_not_in_mathys =
      sum(
        unresolved$annotation_not_in_mathys == TRUE
      ),

    unresolved_without_annotation_symbol =
      nrow(no_annotation),

    unique_harmonized_symbols =
      uniqueN(
        resolved$harmonized_symbol
      ),

    duplicated_harmonized_symbols =
      nrow(duplicates)
  )

  print(summary)


  cat(
    "\nUnresolved com símbolo anotado,",
    "mas fora da Mathys:",
    nrow(annotated_not_mathys),
    "\n"
  )

  cat(
    "Unresolved sem símbolo recuperável:",
    nrow(no_annotation),
    "\n"
  )

  cat(
    "Símbolos harmonizados duplicados:",
    nrow(duplicates),
    "\n"
  )


  # ----------------------------------------------------------
  # Salvar auditorias
  # ----------------------------------------------------------

  fwrite(
    unresolved,
    file.path(
      out_dir,
      paste0(
        dataset,
        "_unresolved_detailed.tsv"
      )
    ),
    sep = "\t"
  )

  fwrite(
    annotated_not_mathys,
    file.path(
      out_dir,
      paste0(
        dataset,
        "_annotated_but_not_mathys.tsv"
      )
    ),
    sep = "\t"
  )

  fwrite(
    no_annotation,
    file.path(
      out_dir,
      paste0(
        dataset,
        "_without_annotation_symbol.tsv"
      )
    ),
    sep = "\t"
  )

  fwrite(
    duplicates,
    file.path(
      out_dir,
      paste0(
        dataset,
        "_duplicated_symbols.tsv"
      )
    ),
    sep = "\t"
  )

  all_summary[[dataset]] <- summary
}


# ============================================================
# RESUMO FINAL
# ============================================================

cat("\n============================================\n")
cat("RESUMO FINAL\n")
cat("============================================\n\n")

final_summary <- rbindlist(
  all_summary
)

print(final_summary)

fwrite(
  final_summary,
  file.path(
    out_dir,
    "unresolved_audit_summary.tsv"
  ),
  sep = "\t"
)

cat(
  "\nAuditoria concluída.\n",
  "Nenhuma matriz de expressão foi modificada.\n",
  sep = ""
)
