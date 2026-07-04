#!/usr/bin/env Rscript
# covid_robustness_2021q3.R -- robustness check: does adding the 2021Q3 Delta
# lockdown quarter (AU GDP ~ -1.9%, a ~4-sigma event) to the treated COVID
# quarters change the evaluation?  (README.md D17/D22.)
#
# Runs the recursive OOS with covid.quarters extended to include 2021-09-01
# for the members the treatment can affect (constant-volatility engines +
# rw/ar4/ucmean; the SV members and ucsv use t-errors, not LP weights, and are
# unaffected by construction), then compares mean scores against the baseline
# run and writes reports/covid_2021q3_robustness.csv.
#
# The modified config gets its own cache directory automatically (the OOS
# cache is keyed by a hash of the estimation-relevant config), so this run
# never contaminates the production cache. Compute: ~ the cost of one main
# pipeline run for the affected members. Usage:
#   Rscript scripts/covid_robustness_2021q3.R

if (file.exists(".Renviron")) readRenviron(".Renviron")
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
setup_logging()

cfg <- load_config("config/config.yml")
baseline <- "output/tables/scores_by_horizon.csv"
if (!file.exists(baseline))
  stop("baseline scores not found -- run the main pipeline first")

# extended treatment window
cfg2 <- cfg
cfg2$covid$quarters <- c(cfg$covid$quarters, "2021-09-01")

spec <- build_transform_spec(cfg2)
raw  <- get_raw_data(cfg2, spec)
td   <- transform_data(raw, spec)
check_data(td, spec)

# members the LP treatment can affect
affected <- function(m) !(m$engine %in% c("sv", "ucsv"))
members <- Filter(affected, all_members(cfg2))
log_info("2021Q3 robustness: {length(members)} affected members")

future::plan(future::multisession, workers = cfg2$parallel$workers)
oos <- setNames(
  lapply(members, function(m)
    timed(paste0("OOS(2021Q3) ", m$name), run_oos_member(m, td, spec, cfg2))),
  vapply(members, `[[`, "", "name"))
future::plan(future::sequential)

scores2 <- do.call(rbind, lapply(names(oos), function(m)
  score_member(m, oos[[m]], td, spec, cfg2)))
sm2 <- summarise_scores(scores2)

sm1 <- read.csv(baseline, stringsAsFactors = FALSE)
sm1 <- sm1[sm1$member %in% names(oos), ]
cmp <- merge(sm1, sm2, by = c("member", "variable", "measure", "h"),
             suffixes = c("_base", "_2021q3"))
cmp$crps_ratio <- cmp$crps_2021q3 / cmp$crps_base
cmp$logdens_diff <- cmp$logdens_2021q3 - cmp$logdens_base
dir.create("reports", showWarnings = FALSE)
write.csv(cmp, "reports/covid_2021q3_robustness.csv", row.names = FALSE)

# headline summary to the console
agg <- aggregate(cbind(crps_ratio, logdens_diff) ~ member,
                 cmp[cmp$measure == "q", ], mean)
print(agg, row.names = FALSE)
log_info(paste0(
  "2021Q3 robustness written to reports/covid_2021q3_robustness.csv; ",
  "crps_ratio ~ 1 and logdens_diff ~ 0 mean the extra treated quarter ",
  "does not matter"))
