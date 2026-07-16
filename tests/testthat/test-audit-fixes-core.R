source_project <- function() {
  root <- testthat::test_path("..", "..")
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
    sys.source(f, envir = globalenv())
}
source_project()

# ---- FIX 1a: ess_basic (dependency-free ESS) --------------------------------

test_that("ess_basic is near n for iid chains", {
  set.seed(101)
  x <- rnorm(1000)
  es <- ess_basic(x)
  expect_true(is.finite(es))
  expect_gte(es, 0.5 * length(x))
  expect_lte(es, length(x))
})

test_that("ess_basic is far below n for a strongly autocorrelated AR(1) chain", {
  set.seed(102)
  e <- rnorm(1000)
  x <- as.numeric(filter(e, 0.95, method = "recursive"))
  es <- ess_basic(x)
  expect_true(is.finite(es))
  expect_lte(es, 0.25 * length(x))
  expect_gte(es, 1)
})

test_that("ess_basic agrees with coda::effectiveSize within a factor of 2", {
  testthat::skip_if_not_installed("coda")
  set.seed(103)
  x_iid <- rnorm(1000)
  e <- rnorm(1000)
  x_ar <- as.numeric(filter(e, 0.95, method = "recursive"))
  for (x in list(x_iid, x_ar)) {
    a <- ess_basic(x)
    b <- as.numeric(coda::effectiveSize(x))
    ratio <- a / b
    expect_true(ratio > 0.5 && ratio < 2,
               info = paste("ess_basic =", a, "coda =", b))
  }
})

test_that("ess_basic is deterministic", {
  set.seed(104)
  x <- as.numeric(filter(rnorm(500), 0.7, method = "recursive"))
  expect_identical(ess_basic(x), ess_basic(x))
})

test_that("ess_basic returns NA for short or degenerate chains", {
  expect_true(is.na(ess_basic(rnorm(5))))
  expect_true(is.na(ess_basic(rep(1, 100))))
})

# ---- FIX 1b/c/d: THE INVARIANT -- coda never changes posterior draws --------
# fit_sv's adaptive-thinning retry must fire (or not) identically whether or
# not coda is installed, because the retry decision is gated on ess_basic()
# (built-in), not safe_ess() (coda). Verified on a dataset that does NOT
# trigger the retry and one that DOES (mirrors the audit verifier's trigger:
# small T, near-unit-root domestic block + a late outlier).

.core_spec_m <- function() {
  data.frame(variable = c("f", "d"), block = c("foreign", "domestic"),
            delta = c(1, 1), transform = "level", stringsAsFactors = FALSE)
}

# T = 30, p = 4, near-unit-root domestic block with a -8 spike 3 quarters
# from the end: at fit_seed 999 this reliably drives at least one equation's
# own-lag ESS below the retry gate (verified empirically).
.mk_sv_data_hard <- function(seed, Tn = 30) {
  set.seed(seed)
  e1 <- rnorm(Tn); e2 <- rnorm(Tn)
  y1 <- as.numeric(filter(e1, 0.97, method = "recursive"))
  y2 <- as.numeric(filter(e2, 0.95, method = "recursive")) + 0.3 * y1
  y2[Tn - 2] <- y2[Tn - 2] - 8
  cbind(f = y1, d = y2)
}

# well-behaved, well-mixing data: no retry expected.
.mk_sv_data_easy <- function(seed, Tn = 80) {
  set.seed(seed)
  e1 <- rnorm(Tn); e2 <- rnorm(Tn)
  cbind(f = e1, d = 0.2 * e1 + e2)
}

.run_fit_sv_capture <- function(y, member, cfg, disable_coda, fit_seed) {
  # force the base-R message() fallback (soe.no_logger) so the retry's
  # log_debug() is interceptable via withCallingHandlers below, regardless
  # of whether the `logger` package (which writes straight to the console,
  # bypassing R's message condition system) happens to be installed.
  opts <- options(soe.disable_coda = disable_coda, soe.no_logger = TRUE)
  on.exit(options(opts))
  log_threshold("DEBUG")   # the fallback's own threshold, set after the
                           # soe.no_logger option above takes effect
  msgs <- character(0)
  res <- withCallingHandlers(
    { set.seed(fit_seed); fit_sv(y, member, .core_spec_m(), cfg,
                                 list(lambda = 0.2)) },
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  list(post = res, retried = any(grepl("retrying", msgs)))
}

test_that("fit_sv draws are bit-identical with/without coda (no retry triggered)", {
  testthat::skip_if_not_installed("stochvol")
  on.exit(log_threshold("INFO"), add = TRUE)

  y <- .mk_sv_data_easy(1, 80)
  member <- list(name = "sv", lags = 1)
  cfg <- list(mcmc = list(ndraw = 150, nburn = 50), covid = list(sv_t_errors = TRUE))

  a <- .run_fit_sv_capture(y, member, cfg, disable_coda = FALSE, fit_seed = 42)
  b <- .run_fit_sv_capture(y, member, cfg, disable_coda = TRUE,  fit_seed = 42)

  expect_false(a$retried)   # sanity: this dataset is the "no retry" case
  expect_false(b$retried)
  expect_identical(lapply(a$post$eqs, `[[`, "beta"), lapply(b$post$eqs, `[[`, "beta"))
  expect_identical(lapply(a$post$eqs, `[[`, "hT"),   lapply(b$post$eqs, `[[`, "hT"))
})

test_that("fit_sv draws are bit-identical with/without coda (retry triggered)", {
  testthat::skip_if_not_installed("stochvol")
  on.exit(log_threshold("INFO"), add = TRUE)

  y <- .mk_sv_data_hard(1, 30)
  member <- list(name = "sv", lags = 4)
  cfg <- list(mcmc = list(ndraw = 150, nburn = 150), covid = list(sv_t_errors = TRUE))

  a <- .run_fit_sv_capture(y, member, cfg, disable_coda = FALSE, fit_seed = 999)
  b <- .run_fit_sv_capture(y, member, cfg, disable_coda = TRUE,  fit_seed = 999)

  expect_true(a$retried)    # sanity: this dataset/seed DOES trigger a retry
  expect_true(b$retried)    # -- and identically so without coda
  expect_identical(lapply(a$post$eqs, `[[`, "beta"), lapply(b$post$eqs, `[[`, "beta"))
  expect_identical(lapply(a$post$eqs, `[[`, "hT"),   lapply(b$post$eqs, `[[`, "hT"))
})

# ---- FIX 4: config_hash stops on a missing estimation file ------------------

test_that("config_hash stops with an informative error when a file is missing", {
  tmp_r <- file.path(tempdir(), "audit_core_missing_r")
  unlink(tmp_r, recursive = TRUE)
  dir.create(tmp_r)
  on.exit(unlink(tmp_r, recursive = TRUE), add = TRUE)
  root <- testthat::test_path("..", "..")
  src_files <- c("utils", "data_sources", "transforms", "priors", "engines",
                "forecast", "benchmarks", "evaluate", "covid")
  for (nm in src_files) {
    src <- file.path(root, "R", paste0(nm, ".R"))
    if (nm != "engines") file.copy(src, file.path(tmp_r, paste0(nm, ".R")))
  }
  cfg <- list(master_seed = 1, data = NULL, synthetic = NULL, variables = NULL,
             horizons = NULL, mcmc = NULL, glp = NULL, covid = NULL,
             suite = NULL, benchmarks = NULL, evaluation = NULL)
  expect_error(config_hash(cfg, r_dir = tmp_r), "engines\\.R")
})

# ---- FIX 2: conj_br stability on the IMPLIED JOINT reduced form -------------

test_that("conj_br stability share reflects the joint system, not just the foreign block", {
  # Hand-built two-block posterior (nf = 2 foreign, nd = 2 domestic, p = 1):
  # the foreign marginal VAR alone is stable, but the domestic conditional's
  # own-lag (1.06) plus strong feedback through the contemporaneous-foreign
  # loading makes the IMPLIED JOINT system explosive.
  nf <- 2; nd <- 2; M <- 4; p <- 1; n <- 50
  Kf <- 1 + nf * p
  Bf <- array(0, c(n, Kf, nf))
  Bf_ <- matrix(c(0, 0, 0.3, 0.05, 0.02, 0.25), Kf, nf, byrow = TRUE)
  for (d in seq_len(n)) Bf[d, , ] <- Bf_
  Sf <- array(0, c(n, nf, nf))
  for (d in seq_len(n)) Sf[d, , ] <- diag(nf) * 0.5

  Kd <- 1 + M * p + nf
  Bd <- array(0, c(n, Kd, nd))
  Bd_ <- matrix(0, Kd, nd)
  Bd_[4, 1] <- 1.06   # domestic-1 own lag: explosive on its own
  Bd_[5, 2] <- 0.2    # domestic-2 own lag: mild
  Bd_[6, 1] <- 0.6    # feedback: domestic-1 on contemporaneous foreign-1
  Bd_[7, 2] <- 0.1
  for (d in seq_len(n)) Bd[d, , ] <- Bd_
  Sd <- array(0, c(n, nd, nd))
  for (d in seq_len(n)) Sd[d, , ] <- diag(nd) * 0.5

  postf <- list(B = Bf, Sigma = Sf)
  postd <- list(B = Bd, Sigma = Sd)

  # old (pre-fix) diagnostic: foreign marginal only
  stable_old <- mean(vapply(seq_len(min(n, 200)), function(d)
    max_eig_mod(matrix(postf$B[d, , ], ncol = nf), nf, p) < 1.05, logical(1)))
  # new (fixed) diagnostic: implied joint reduced form via .conj_br_joint
  stable_new <- mean(vapply(seq_len(min(n, 200)), function(d) {
    jr <- .conj_br_joint(postf$B[d, , ], postf$Sigma[d, , ],
                         postd$B[d, , ], postd$Sigma[d, , ], M, p, nf)
    max_eig_mod(jr$B, M, p) < 1.05
  }, logical(1)))

  expect_equal(stable_old, 1)      # the bug: foreign-only looked perfectly stable
  expect_lt(stable_new, 0.5)       # the fix: the joint system is flagged explosive
})

test_that("fit_conj_br runs end-to-end on well-behaved data (regression check)", {
  set.seed(21)
  Tn <- 60
  y <- cbind(f1 = cumsum(rnorm(Tn, sd = 0.1)) + rnorm(Tn, sd = 0.3),
            f2 = rnorm(Tn), d1 = rnorm(Tn), d2 = rnorm(Tn))
  spec_m <- data.frame(variable = c("f1", "f2", "d1", "d2"),
                       block = c("foreign", "foreign", "domestic", "domestic"),
                       delta = c(0, 0, 0, 0), transform = "level",
                       stringsAsFactors = FALSE)
  member <- list(name = "cb", lags = 1, block_exog = TRUE)
  cfg <- list(mcmc = list(forecast_draws = 100))
  prior <- list(lambda = 0.2, soc = FALSE, dio = FALSE)
  post <- fit_conj_br(y, member, spec_m, cfg, prior)
  expect_true(is.finite(post$diagnostics$stable_share))
  expect_gte(post$diagnostics$stable_share, 0)
  expect_lte(post$diagnostics$stable_share, 1)
})

test_that("fit_conj_br handles the nf == 1 / nd == 1 degenerate block sizes", {
  # a bare `arr[d, , ]` silently drops to a vector when a block has exactly
  # one variable; the joint-system diagnostic must reshape explicitly.
  set.seed(22)
  Tn <- 50
  y <- cbind(f = cumsum(rnorm(Tn, sd = 0.1)) + rnorm(Tn, sd = 0.3),
            d = rnorm(Tn))
  spec_m <- data.frame(variable = c("f", "d"),
                       block = c("foreign", "domestic"),
                       delta = c(0, 0), transform = "level",
                       stringsAsFactors = FALSE)
  member <- list(name = "cb", lags = 1, block_exog = TRUE)
  cfg <- list(mcmc = list(forecast_draws = 100))
  prior <- list(lambda = 0.2, soc = FALSE, dio = FALSE)
  post <- fit_conj_br(y, member, spec_m, cfg, prior)
  expect_true(is.finite(post$diagnostics$stable_share))
  expect_gte(post$diagnostics$stable_share, 0)
  expect_lte(post$diagnostics$stable_share, 1)
})

# ---- FIX 3: derive_seed unchanged (28-bit space is documented, not altered) --

test_that("derive_seed computation is unchanged", {
  expect_equal(derive_seed(12345, "member::small::2020-01-01"),
              (strtoi(substr(digest::digest("12345::member::small::2020-01-01",
                                            algo = "xxhash32"), 1, 7), base = 16L) +
                 12345) %% .Machine$integer.max)
})
