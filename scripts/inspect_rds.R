files <- list.files("data/processed", pattern = "\\.rds$", full.names = TRUE)
for (f in files) {
  cat("\n==================================================\n")
  cat("Arquivo:", f, "\n")
  obj <- readRDS(f)
  cat("Formato do Objeto:", class(obj), "\n")
  cat("Genes retidos:", nrow(obj$counts), "\n")
  cat("Amostras retidas:", ncol(obj$counts), "\n")
  cat("Primeiros 5 genes:", paste(head(rownames(obj$counts), 5), collapse = ", "), "\n")
}
