source_project <- function() {
  root <- testthat::test_path("..", "..")
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
    sys.source(f, envir = globalenv())
}
source_project()

# ---- FIX 1: .assert_nonempty guard shared by every fetch_series provider ----

test_that(".assert_nonempty passes nonempty data.frames through unchanged", {
  df <- data.frame(date = as.Date("2024-01-01"), value = 1)
  expect_identical(.assert_nonempty(df, "RBA", "SERIES1"), df)
})

test_that(".assert_nonempty errors, naming provider and id, for every provider", {
  empty <- data.frame(date = as.Date(character()), value = numeric())
  expect_error(.assert_nonempty(empty, "RBA", "FIRMMCRTD"),
               "RBA series 'FIRMMCRTD' returned no observations")
  expect_error(.assert_nonempty(empty, "ABS", "A2304200A"),
               "ABS series 'A2304200A' returned no observations")
  expect_error(.assert_nonempty(empty, "FRED", "GDPC1"),
               "FRED series 'GDPC1' returned no observations")
  expect_error(.assert_nonempty(empty, "dbnomics", "SOME/ID"),
               "dbnomics series 'SOME/ID' returned no observations")
})

# ---- FIX 2: check_data loglevel_range ----

mk_spec <- function() {
  data.frame(variable = c("rtwi", "cash_rate", "gdp_growth"),
             transform = c("loglevel", "level", "dlog"),
             stringsAsFactors = FALSE)
}

mk_panel <- function(rtwi_val) {
  data.frame(date = seq(as.Date("2000-01-01"), by = "quarter", length.out = 90),
             rtwi = rep(rtwi_val, 90),
             cash_rate = rep(3.5, 90),
             gdp_growth = rep(0.5, 90))
}

test_that("check_data passes with a plausible loglevel value (100*log index)", {
  spec <- mk_spec()
  td <- mk_panel(460)  # 100*log(index); plausible
  expect_true(isTRUE(check_data(td, spec)$loglevel_range))
})

test_that("check_data flags loglevel_range for a raw (unlogged) index", {
  spec <- mk_spec()
  td <- mk_panel(65)   # someone forgot the 100*log transform
  expect_error(check_data(td, spec), "loglevel_range")
})

test_that("check_data flags loglevel_range for a doubled 100x scaling", {
  spec <- mk_spec()
  td <- mk_panel(46000)  # e.g. 100*log applied twice, or an extra 100x
  expect_error(check_data(td, spec), "loglevel_range")
})

# ---- FIX 3: transform_data trailing/interior NA handling unchanged ----

test_that("transform_data trims to the balanced sample on a trailing NA", {
  spec <- data.frame(variable = c("a", "b"), transform = c("level", "level"),
                     stringsAsFactors = FALSE)
  raw <- list(data = data.frame(
    date = seq(as.Date("2020-01-01"), by = "quarter", length.out = 5),
    a = c(1, 2, 3, 4, 5),
    b = c(1, 2, 3, 4, NA)))   # trailing NA (publication lag at the frontier)
  out <- transform_data(raw, spec)
  expect_equal(nrow(out), 4)
  expect_false(anyNA(out[, -1]))
})

test_that("transform_data still stops with 'interior missing values' on an interior gap", {
  spec <- data.frame(variable = c("a", "b"), transform = c("level", "level"),
                     stringsAsFactors = FALSE)
  raw <- list(data = data.frame(
    date = seq(as.Date("2020-01-01"), by = "quarter", length.out = 5),
    a = c(1, 2, 3, 4, 5),
    b = c(1, 2, NA, 4, 5)))   # interior NA, bounded by finite rows
  expect_error(transform_data(raw, spec), "interior missing values")
})
