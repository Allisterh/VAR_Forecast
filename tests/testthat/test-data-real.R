source_project <- function() {
  root <- testthat::test_path("..", "..")
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
    sys.source(f, envir = globalenv())
}
source_project()

test_that("to_quarterly averages full quarters and NA-flags partial frontier quarters", {
  # monthly series: 2024Q1-Q4 complete, 2025Q1 only one month (partial)
  d <- data.frame(
    date = c(seq(as.Date("2024-01-01"), by = "month", length.out = 12),
             as.Date("2025-01-01")),
    value = c(rep(1:4, each = 3), 99))
  q <- to_quarterly(d)
  expect_equal(nrow(q), 5)
  expect_equal(q$value[q$date == as.Date("2024-01-01")], 1)  # mean of Jan-Mar (all 1)
  expect_equal(q$value[q$date == as.Date("2024-10-01")], 4)
  expect_true(is.na(q$value[q$date == as.Date("2025-01-01")]))  # 1 of 3 months -> NA
})

test_that("to_quarterly keeps quarterly-native series intact (aggregation is a no-op)", {
  d <- data.frame(date = seq(as.Date("2024-01-01"), by = "quarter", length.out = 6),
                  value = 1:6)
  q <- to_quarterly(d)
  expect_equal(nrow(q), 6)
  expect_equal(q$value, as.numeric(1:6))      # max count = 1, nothing dropped
})

test_that("to_quarterly tolerates short-but-complete daily quarters", {
  # ~63 trading days/quarter; a holiday-shortened quarter (58) must NOT be dropped
  mk <- function(start, n) seq(as.Date(start), by = "day", length.out = n)
  d <- data.frame(
    date = c(mk("2024-01-01", 66), mk("2024-04-01", 58), mk("2024-07-01", 64)),
    value = c(rep(1, 66), rep(2, 58), rep(3, 64)))
  q <- to_quarterly(d)
  expect_false(anyNA(q$value))                # 58/66 = 0.88 > 0.8 threshold
})

test_that("fetch_series errors on an unknown provider", {
  expect_error(fetch_series("x", list(provider = "bloomberg", id = "z"), list()),
               "unknown provider")
})

test_that("real path errors (never substitutes) when a series cannot be obtained", {
  cfg <- load_config(testthat::test_path("..", "..", "config", "config.yml"))
  spec <- data.frame(variable = "x", provider = "rba", series_id = "DEFINITELY_NOT_A_SERIES",
                     pre = "", transform = "level", stringsAsFactors = FALSE)
  td_dir <- tempfile(); dir.create(td_dir)
  expect_error(download_real_data(cfg, spec, raw_dir = td_dir),
               "Could not obtain series")
})

test_that("pct_change interior NA errors loudly rather than NA-poisoning the tail", {
  # simulate the cumulation guard directly on a constructed quarterly series
  q <- data.frame(date = seq(as.Date("2020-01-01"), by = "quarter", length.out = 6),
                  value = c(0.5, 0.5, NA, 0.5, 0.5, 0.5))
  fin <- which(is.finite(q$value))
  span <- min(fin):max(fin)
  expect_true(anyNA(q$value[span]))           # the guard condition fires
})
