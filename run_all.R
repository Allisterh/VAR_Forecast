#!/usr/bin/env Rscript
# run_all.R -- clone-to-result convenience wrapper:
#   Rscript run_all.R
# Restores pinned dependencies if renv is present, then runs the pipeline.
#
# The DEFAULT is real data, which needs FRED_API_KEY (foreign block) in a
# gitignored .Renviron. For a fully offline run set data.source: synthetic in
# config/config.yml.

if (file.exists(".Renviron")) readRenviron(".Renviron")

if (file.exists("renv.lock") && requireNamespace("renv", quietly = TRUE)) {
  message("Restoring renv library (first run may take a while)...")
  options(renv.config.install.verbose = FALSE)
  try(renv::restore(prompt = FALSE), silent = TRUE)
}

if (!requireNamespace("targets", quietly = TRUE))
  stop("the 'targets' package is required: install.packages('targets')")

cfg_src <- tryCatch(yaml::read_yaml("config/config.yml")$default$data$source,
                    error = function(e) NA)
if (identical(cfg_src, "real") && !nzchar(Sys.getenv("FRED_API_KEY")))
  message("NOTE: data.source is 'real' but FRED_API_KEY is not set. ",
          "Put it in .Renviron, or set data.source: synthetic for offline.")

targets::tar_make()

cat("\nDone. Outputs:\n",
    "  output/tables/     evaluation + diagnostics tables\n",
    "  output/figures/    fan charts, PIT, weights, score plots\n",
    "  output/forecasts/  forecast tables\n",
    "  reports/report.html (if quarto available)\n")
