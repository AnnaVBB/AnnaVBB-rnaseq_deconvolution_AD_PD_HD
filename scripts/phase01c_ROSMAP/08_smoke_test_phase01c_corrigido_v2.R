# 08_smoke_test_phase01c_corrigido_v2.R
# Uso: Rscript 08_smoke_test_phase01c_corrigido_v2.R "/caminho/absoluto/08_integrate_AD_HD_PD_phase01c_corrigido.R"
# Requer data.table e ggplot2 no ambiente; NAO usa nem altera o repositorio real.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Forneca o caminho do script corrigido.")
main_script <- normalizePath(args[1L], mustWork = TRUE)
if (!requireNamespace("data.table", quietly = TRUE)) stop("Instale data.table")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Instale ggplot2")

wd_original <- getwd()
test_root <- tempfile("TCC_test_script_08_")
dir.create(test_root, recursive = TRUE)
setwd(test_root)

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(df, path, sep = "\t", na = "NA")
}

# Gerar dados DESeq2 sem ambiguidade sobre grupos/limiares.
set.seed(1248)
ids <- as.character(1001:1080)
make_de <- function(n_deg, shift = 0, id_col = "GeneID") {
  effects <- rnorm(length(ids), 0, 0.18) + shift
  effects[seq_len(n_deg)] <- ifelse(seq_len(n_deg) %% 3L == 0L, -1.2, 1.2)
  pval <- rep(0.7, length(ids))
  padj <- rep(0.9, length(ids))
  pval[seq_len(n_deg)] <- 0.00001
  padj[seq_len(n_deg)] <- 0.001
  x <- data.frame(GeneID = ids, symbol = paste0("GENE", ids),
                  baseMean = 500, log2FoldChange = effects,
                  lfcSE = 0.1, stat = effects/0.1,
                  pvalue = pval, padj = padj,
                  stringsAsFactors = FALSE)
  names(x)[names(x) == "GeneID"] <- id_col
  x
}
write_tsv(make_de(8), "results/candidate_ROSMAP/deseq2/ROSMAP_deseq2_results.tsv")
write_tsv(make_de(6, 0.02),
          "results/candidate_ROSMAP/deseq2/ROSMAP_deseq2_results_adjusted.tsv")
write_tsv(make_de(20),
          "results/candidate_GSE125583/deseq2/GSE125583_deseq2_results.tsv")
write_tsv(make_de(18, -0.01),
          "results/candidate_GSE125583/deseq2/GSE125583_deseq2_adjusted_results.tsv")
write_tsv(make_de(30, 0, "gene_id"), "results/deseq2/HD_deseq2_results.tsv")
write_tsv(make_de(25, 0, "gene_id"), "results/deseq2/PD_deseq2_results.tsv")

# Iniciar sem proporcoes para testar o estado 'PENDENTE' com MuSiC2 ausente.
cat("Teste isolado em: ", test_root, "\n", sep = "")
output <- suppressWarnings(system2("Rscript", args = shQuote(main_script),
                                   stdout = TRUE, stderr = TRUE))
status <- attr(output, "status")
if (is.null(status)) status <- 0L
if (status != 0L) {
  cat(paste(output, collapse = "\n"))
  stop("Script 08 falhou; status Rscript = ", status)
}

comparison <- "results/phase01c_comparacao_coortes_AD"
stopifnot(file.exists(file.path(comparison, "08_criterios_multicriterio_AD.tsv")))
stopifnot(file.exists(file.path(comparison, "08_two_candidate_integration_summary.tsv")))
criteria <- data.table::fread(file.path(comparison, "08_criterios_multicriterio_AD.tsv"))
stopifnot(nrow(criteria) == 2L)
stopifnot(setequal(criteria$dataset, c("ROSMAP", "GSE125583")))
stopifnot(all(criteria$decision_status == "PENDENTE_DOS_TRES_METODOS"))
stopifnot(file.exists(file.path(comparison, "08_method_availability.tsv")))
stopifnot(file.exists(file.path(comparison, "08_DE_model_stability.tsv")))

for (cohort in c("ROSMAP", "GSE125583")) {
  root <- file.path("results", paste0("candidate_", cohort),
                    "integration_candidate_AD_HD_PD")
  for (filename in c("venn_final_primary_AD_HD_PD.png",
                     "venn_final_strict_AD_HD_PD.png",
                     "venn_final_article_like_nominal_AD_HD_PD.png",
                     "log2FC_AD_vs_HD_common_universe.png",
                     "log2FC_AD_vs_PD_common_universe.png",
                     "log2FC_HD_vs_PD_common_universe.png")) {
    stopifnot(file.exists(file.path(root, "plots", filename)))
  }
  tab <- data.table::fread(file.path(root, "tables",
                                     "final_AD_HD_PD_integration_summary.tsv"))
  vals <- setNames(tab$value, tab$metric)
  stopifnot(vals[["common_testable_universe_geneids"]] == 80)
  # Resultado conhecido: 8 DEGs ROSMAP vs 20 GSE125583 e genes 1001..1080
  expected_triple <- if (cohort == "ROSMAP") 8L else 20L
  stopifnot(vals[["AD_HD_PD_primary_shared"]] == expected_triple)
  stopifnot(vals[["AD_primary"]] == expected_triple)
  stopifnot(vals[["HD_primary"]] == 30L)
  stopifnot(vals[["PD_primary"]] == 25L)
  stopifnot(file.exists(file.path(root, "tables", "final_primary_shared_AD_HD_PD.tsv")))
  stopifnot(file.exists(file.path(root, paste0("validation_",
                    if (cohort == "ROSMAP") "GSE125583" else "ROSMAP"),
                    if (cohort == "ROSMAP") "GSE125583_validation_summary.tsv"
                    else "ROSMAP_validation_summary.tsv")))
}

cat("PASSOU: execucao completa com dados sinteticos e metodos pendentes.\n")
cat("Arquivos em ", test_root, "\n", sep = "")
setwd(wd_original)
