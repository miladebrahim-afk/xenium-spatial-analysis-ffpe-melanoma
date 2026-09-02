# CellChat is separate because it is computationally heavy and version-sensitive.
# Run after the publication-reference meta-niche objects have been generated.

run_script <- function(path) {
  cat("\n============================================================\n")
  cat("RUNNING: ", path, "\n", sep = "")
  cat("============================================================\n")
  source(path, echo = FALSE, chdir = FALSE)
  invisible(gc())
}

run_script("scripts/04_cellchat/01_global_cellchat_nf1.R")
run_script("scripts/04_cellchat/02_meta_niche_cellchat_nf1.R")
run_script("scripts/04_cellchat/03_meta_niche_cellchat_figures.R")
run_script("scripts/04_cellchat/04_global_cellchat_figures.R")

cat("\nCellChat workflow complete. Outputs are under results/cellchat/.\n")
