source_project <- function() {
  root <- testthat::test_path("..", "..")
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
    sys.source(f, envir = globalenv())
}
source_project()

# small "add" harness matching write_model_scorecard()'s closure, so we can
# capture and inspect the markdown lines .combo_vs_best()/.winner_matrix()
# emit.
.capture <- function(fn, ...) {
  L <- c()
  add <- function(...) L <<- c(L, paste0(...))
  fn(add, ...)
  L
}

# ---- .combo_vs_best (FIX 1) --------------------------------------------------

.combo_view <- function() {
  v <- expand.grid(member = c("A", "B", "C"), h = 1:3, stringsAsFactors = FALSE)
  v$crps <- c(1, 1.5, 2,      # h=1 (near)
              2, 2.5, 3,      # h=2 (medium)
              3, 3.5, 4)      # h=3 (far)
  v
}
.combo_buckets <- list(`near (1-4)` = 1, `medium (5-8)` = 2, `far (9-12)` = 3)

test_that(".combo_vs_best builds a table with all-finite input, best bolded", {
  v <- .combo_view()
  L <- NULL
  expect_error(L <- .capture(.combo_vs_best, v, .combo_buckets), NA)
  a_row <- grep("^\\| A \\|", L, value = TRUE)
  expect_length(a_row, 1)
  expect_match(a_row, "\\*\\*1\\.000\\*\\*")  # A is best (lowest) in every column
  expect_false(any(grepl("—", L)))
})

test_that(".combo_vs_best survives a single NA cell: no crash, cell -> em-dash, others unaffected", {
  v <- .combo_view()
  v$crps[v$member == "B" & v$h == 1] <- NA_real_
  L <- NULL
  expect_error(L <- .capture(.combo_vs_best, v, .combo_buckets), NA)
  b_row <- grep("^\\| B \\|", L, value = TRUE)
  expect_length(b_row, 1)
  cells <- strsplit(trimws(b_row), "\\|")[[1]]
  cells <- trimws(cells[cells != ""])
  expect_equal(cells[1], "B")
  expect_equal(cells[2], "—")             # near: the wiped bucket
  expect_false(grepl("—", cells[3]))      # medium: unaffected
  expect_false(grepl("—", cells[4]))      # far: unaffected
  # A's row is untouched and still bolds the near column
  a_row <- grep("^\\| A \\|", L, value = TRUE)
  expect_match(a_row, "\\*\\*1\\.000\\*\\*")
})

test_that(".combo_vs_best survives a single NaN cell the same way as NA", {
  v <- .combo_view()
  v$crps[v$member == "C" & v$h == 3] <- NaN
  L <- NULL
  expect_error(L <- .capture(.combo_vs_best, v, .combo_buckets), NA)
  c_row <- grep("^\\| C \\|", L, value = TRUE)
  expect_length(c_row, 1)
  cells <- strsplit(trimws(c_row), "\\|")[[1]]
  cells <- trimws(cells[cells != ""])
  expect_equal(cells[4], "—")             # far: the NaN'd bucket
  expect_false(grepl("—", cells[2]))
  expect_false(grepl("—", cells[3]))
})

test_that(".combo_vs_best handles an entirely-NA member column without crashing", {
  v <- .combo_view()
  v$crps[v$member == "B"] <- NA_real_
  L <- NULL
  expect_error(L <- .capture(.combo_vs_best, v, .combo_buckets), NA)
  b_row <- grep("^\\| B \\|", L, value = TRUE)
  cells <- strsplit(trimws(b_row), "\\|")[[1]]
  cells <- trimws(cells[cells != ""])
  expect_true(all(cells[2:4] == "—"))
})

# ---- .winner_matrix (FIX 2) --------------------------------------------------

test_that(".winner_matrix renders em-dash (not a crash) when a cell has zero finite crps", {
  scores <- data.frame(
    variable = c("gdp", "gdp", "gdp", "gdp"),
    measure  = "level",
    member   = c("m1", "m2", "m1", "m2"),
    h        = c(1, 1, 3, 4),
    crps     = c(NA_real_, NaN, 0.4, 0.5))
  buckets <- list(`near (1-4)` = 1:2, `far (9-12)` = 3:4)
  out <- NULL
  expect_error(out <- .winner_matrix(scores, "gdp", buckets), NA)
  rows <- strsplit(out, "\n")[[1]]
  gdp_row <- grep("^\\| gdp \\|", rows, value = TRUE)
  expect_length(gdp_row, 1)
  cells <- strsplit(trimws(gdp_row), "\\|")[[1]]
  cells <- trimws(cells[cells != ""])
  expect_equal(cells[2], "—")            # near: both rows NA/NaN -> no finite crps
  expect_match(cells[3], "^m1 \\(0\\.400\\)$")  # far: normal winner cell
})

test_that(".winner_matrix renders em-dash when a (variable, bucket) has zero rows at all", {
  scores <- data.frame(variable = "unemp", measure = "level", member = "m1",
                       h = 1, crps = 0.2, stringsAsFactors = FALSE)
  buckets <- list(`near (1-4)` = 1:4, `far (9-12)` = 9:12)  # far has no rows
  out <- NULL
  expect_error(out <- .winner_matrix(scores, "unemp", buckets), NA)
  rows <- strsplit(out, "\n")[[1]]
  row <- grep("^\\| unemp \\|", rows, value = TRUE)
  cells <- strsplit(trimws(row), "\\|")[[1]]
  cells <- trimws(cells[cells != ""])
  expect_equal(cells[3], "—")
})

# ---- read_forecast_profiles (FIX 3) -----------------------------------------

.write_csv <- function(df) {
  f <- tempfile(fileext = ".csv")
  write.csv(df, f, row.names = FALSE)
  f
}

test_that("read_forecast_profiles errors loudly on a non-numeric value cell, naming variable + raw value", {
  df <- data.frame(variable = c("f_act", "f_act"),
                   date = c("2026-04-01", "2026-07-01"),
                   value = c("TBD", "0.50"), stringsAsFactors = FALSE)
  f <- .write_csv(df)
  expect_error(
    read_forecast_profiles(f, as.Date("2026-01-01"), 12),
    regexp = "f_act.*TBD")
})

test_that("read_forecast_profiles errors loudly on duplicate (variable, h) rows", {
  df <- data.frame(variable = c("f_act", "f_act"), h = c(1L, 1L),
                   value = c(0.5, 0.6))
  f <- .write_csv(df)
  expect_error(
    read_forecast_profiles(f, as.Date("2026-01-01"), 12),
    regexp = "duplicate")
})

test_that("read_forecast_profiles errors loudly when two dates map to the same quarter", {
  df <- data.frame(variable = c("cpi_inflation", "cpi_inflation"),
                   date = c("2026-04-01", "2026-05-15"),
                   value = c(0.5, 0.6), stringsAsFactors = FALSE)
  f <- .write_csv(df)
  expect_error(
    read_forecast_profiles(f, as.Date("2026-01-01"), 12),
    regexp = "duplicate")
})

test_that("read_forecast_profiles parses a clean CSV with a trailing note column", {
  df <- data.frame(
    variable = c("f_act", "f_act"),
    date = c("2026-04-01", "2026-07-01"),
    value = c(0.45, 0.50),
    note = c("illustrative path", ""), stringsAsFactors = FALSE)
  f <- .write_csv(df)
  p <- NULL
  expect_error(p <- read_forecast_profiles(f, as.Date("2026-01-01"), 12), NA)
  expect_true(is.list(p))
  expect_true("f_act" %in% names(p))
  expect_length(p$f_act, 12)
  expect_equal(p$f_act[1], 0.45)
  expect_equal(p$f_act[2], 0.50)
  expect_true(all(is.na(p$f_act[3:12])))
})

test_that("the bundled example_forecast_profiles.csv itself parses cleanly", {
  path <- testthat::test_path("..", "..", "data", "example_forecast_profiles.csv")
  skip_if_not(file.exists(path), "bundled example profiles file not found")
  p <- NULL
  expect_error(p <- read_forecast_profiles(path, as.Date("2026-01-01"), 12), NA)
  expect_true(is.list(p))
  expect_true(all(c("f_act", "f_rate", "cpi_inflation", "unemp_rate") %in% names(p)))
  expect_true(all(vapply(p, is.numeric, logical(1))))
  expect_true(all(vapply(p, length, integer(1)) == 12))
})
