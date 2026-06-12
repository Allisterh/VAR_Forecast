#!/usr/bin/env Rscript
# run_all.R -- clone-to-result convenience wrapper:
#   Rscript run_all.R
# Restores pinned dependencies if renv is present, then runs the pipeline.

if (file.exists("renv.lock") && requireNamespace("renv", quietly = TRUE)) {
  message("Restoring renv library (first run may take a while)...")
  options(renv.config.install.verbose = FALSE)
  try(renv::restore(prompt = FALSE), silent = TRUE)
}

if (!requireNamespace("targets", quietly = TRUE))
  stop("the 'targets' package is required: install.packages('targets')")

targets::tar_make()

cat("\nDone. Outputs:\n",
    "  output/tables/     evaluation + diagnostics tables\n",
    "  output/figures/    fan charts, PIT, weights, score plots\n",
    "  output/forecasts/  forecast tables\n",
    "  reports/report.html (if quarto available)\n")
