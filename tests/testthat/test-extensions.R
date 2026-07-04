# Tests for the 2026-07 extensions: Waggoner-Zha conditional forecasts (D21),
# combination weight-training robustness + crps scheme (D20), the unrestricted
# member (D19), the ar4 delta prior, PIT moment tests, event probabilities.

source_project <- function() {
  root <- testthat::test_path("..", "..")
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
    sys.source(f, envir = globalenv())
}
source_project()

tiny_cfg <- function() {
  cfg <- load_config(testthat::test_path("..", "..", "config", "config.yml"))
  cfg$mcmc$ndraw <- 120; cfg$mcmc$nburn <- 30
  cfg$mcmc$forecast_draws <- 120; cfg$mcmc$store_draws <- 60
  cfg$evaluation$max_origins <- 4
  cfg$glp$enabled <- FALSE
  cfg
}

# ---- Waggoner-Zha conditional forecasting (D21) --------------------------------

test_that("WZ conditioning hits the conditioned path exactly (gibbs)", {
  cfg <- tiny_cfg()
  spec <- build_transform_spec(cfg)
  raw <- generate_synthetic_data(cfg, spec)
  td <- transform_data(raw, spec)
  spec_s <- vars_for_set(spec, "small")
  y <- as.matrix(td[, spec_s$variable])
  member <- list(name = "t", kind = "var", engine = "gibbs", set = "small",
                 lags = 2, prior = list(lambda = 0.2, soc = FALSE, dio = FALSE))
  set.seed(1)
  post <- fit_var_member(y, member, spec_s, cfg)
  path <- c(3.0, 2.5, 2.0, NA)          # partial path: h=4 unconstrained
  cond <- list(variable = "cash_rate", path = path, method = "wz")
  set.seed(2)
  paths <- simulate_paths(post, y, 4, 50, condition = cond)
  cr <- paths[, , "cash_rate"]
  for (s in 1:3) expect_equal(unname(cr[, s]), rep(path[s], 50), tolerance = 1e-8)
  # unconstrained step still stochastic
  expect_gt(sd(cr[, 4]), 1e-3)
  # non-conditioned variables remain stochastic and finite
  expect_true(all(is.finite(paths)))
  expect_gt(sd(paths[, 1, "gdp_growth"]), 1e-3)
})

test_that("WZ conditional mean matches the analytic Gaussian update", {
  # y_t = c + A y_{t-1} + u, u ~ N(0, Sigma). Conditioning y1_{t+1} = r gives
  # E[y2_{t+1} | y1] = m2 + (S21/S11) (r - m1). Build a degenerate 'posterior'
  # with a single known draw and check the WZ machinery against the formula.
  M <- 2; p <- 1
  A <- matrix(c(0.5, 0.1, 0.2, 0.4), 2, 2, byrow = TRUE)   # rows = equations
  cvec <- c(0.3, -0.1)
  Sigma <- matrix(c(1.0, 0.6, 0.6, 0.8), 2, 2)
  B <- rbind(cvec, t(A))                                    # K x M convention
  ndraw <- 4000
  post <- structure(list(engine = "gibbs", M = M, p = p,
                         varnames = c("a", "b"), blocks = c("foreign", "domestic"),
                         B = array(rep(B, each = 1), c(1, 3, 2)),
                         Sigma = array(Sigma, c(1, 2, 2)), ndraw = 1,
                         diagnostics = list()),
                    class = c("post_gibbs", "var_posterior"))
  post$B[1, , ] <- B; post$Sigma[1, , ] <- Sigma
  y <- matrix(c(1, 2), 1, 2, dimnames = list(NULL, c("a", "b")))
  m <- cvec + drop(A %*% y[1, ])                            # unconditional mean
  r <- m[1] + 1.5                                           # condition 1.5 above mean
  set.seed(42)
  paths <- simulate_paths(post, y, 1, ndraw,
                          condition = list(variable = "a", path = r, method = "wz"))
  expect_equal(unname(paths[1, 1, "a"]), r, tolerance = 1e-8)
  mu2 <- m[2] + Sigma[2, 1] / Sigma[1, 1] * (r - m[1])
  v2  <- Sigma[2, 2] - Sigma[2, 1]^2 / Sigma[1, 1]
  expect_equal(mean(paths[, 1, "b"]), mu2, tolerance = 4 * sqrt(v2 / ndraw))
  expect_equal(var(paths[, 1, "b"]), v2, tolerance = 0.08)
})

test_that("conj_br joint-form reconstruction reproduces the two-block simulation", {
  # the joint (B, Sigma) must give the same one-step conditional moments as
  # simulating foreign-then-domestic
  set.seed(3)
  nf <- 2; nd <- 2; M <- 4; p <- 1
  Bf <- rbind(c(0.1, -0.2), matrix(runif(nf * nf, -0.2, 0.4), nf, nf))
  Sf <- crossprod(matrix(rnorm(nf * nf), nf)) / 2 + diag(nf) * 0.3
  Kd <- 1 + M * p + nf
  Bd <- matrix(runif(Kd * nd, -0.2, 0.3), Kd, nd)
  Sd <- crossprod(matrix(rnorm(nd * nd), nd)) / 4 + diag(nd) * 0.2
  jr <- .conj_br_joint(Bf, Sf, Bd, Sd, M, p, nf)
  ystate <- matrix(rnorm(M), 1, M)
  x <- c(1, as.vector(t(ystate)))
  # direct: foreign mean, then domestic mean given expected y_f
  xf <- c(1, ystate[1, seq_len(nf)])
  mf <- drop(crossprod(Bf, xf))
  md <- drop(crossprod(Bd, c(x, mf)))
  mj <- drop(crossprod(jr$B, x))
  expect_equal(mj[seq_len(nf)], mf, tolerance = 1e-10)
  expect_equal(mj[nf + seq_len(nd)], md, tolerance = 1e-10)
  # joint covariance: MC check against the two-stage simulation
  G <- Bd[(1 + M * p + 1):Kd, , drop = FALSE]
  n <- 60000
  ef <- matrix(rnorm(n * nf), n) %*% chol(Sf)
  ed <- matrix(rnorm(n * nd), n) %*% chol(Sd)
  ud <- ef %*% G + ed
  emp <- cov(cbind(ef, ud))
  expect_equal(max(abs(emp - jr$Sigma)), 0, tolerance = 0.05)
})

# ---- combination weight training (D20) -----------------------------------------

.mk_scores <- function(members, origins, hs, v = "v", measure = "q",
                       logdens = -1, crps = 1,
                       date0 = as.Date("2000-01-01")) {
  g <- expand.grid(member = members, origin = origins, h = hs,
                   stringsAsFactors = FALSE)
  g$variable <- v; g$measure <- measure
  g$date <- date0 + 91 * (g$origin + g$h)
  g$logdens <- logdens; g$crps <- crps
  g
}

test_that("crps scheme favours the lower-CRPS member and stays on the simplex", {
  sc <- .mk_scores(c("A", "B"), 1:20, 1:2)
  sc$crps <- ifelse(sc$member == "A", 0.5, 2.0)
  cfg <- list(combination = list(min_train_origins = 3, forgetting = 1,
                                 shrink_kappa = 0))
  w <- combo_weights("crps", sc, "v", 1:2, t = 25, c("A", "B"), cfg)
  expect_equal(sum(w), 1, tolerance = 1e-12)
  expect_equal(unname(w["A"]), 0.8, tolerance = 1e-8)   # (1/0.5)/(1/0.5+1/2)
})

test_that("exclude_covid_train removes COVID realizations from weight training", {
  # member B is catastrophic ONLY in the COVID window; with exclusion its
  # weight must equal member A's
  covq <- as.Date("2020-03-01")
  sc <- .mk_scores(c("A", "B"), 1:20, 1, date0 = as.Date("2015-01-01"))
  # place origins 19,20 realizations inside the exclusion window
  sc$date[sc$origin >= 19] <- covq + 10   # same quarter as 2020Q1
  sc$logdens <- ifelse(sc$member == "B" & sc$origin >= 19, -60, -1)
  base <- list(min_train_origins = 3, forgetting = 1, shrink_kappa = 0)
  cfg_in  <- list(combination = base,
                  covid = list(quarters = list("2020-03-01")))
  cfg_ex  <- modifyList(cfg_in, list(combination = c(base, list(exclude_covid_train = TRUE))))
  w_in <- combo_weights("logscore", sc, "v", 1, t = 25, c("A", "B"), cfg_in)
  w_ex <- combo_weights("logscore", sc, "v", 1, t = 25, c("A", "B"), cfg_ex)
  expect_lt(w_in["B"], 0.01)                    # contaminated: B destroyed
  expect_equal(unname(w_ex["B"]), 0.5, tolerance = 1e-8)  # excluded: equal again
})

test_that("train_measure level trains dlog-variable weights on cum scores", {
  # A wins on q, B wins on cum; level training must follow cum
  members <- c("A", "B")
  q  <- .mk_scores(members, 1:15, 1:4, measure = "q")
  q$logdens <- ifelse(q$member == "A", -1, -3)
  cm <- .mk_scores(members, 1:15, 2:4, measure = "cum")
  cm$logdens <- ifelse(cm$member == "A", -3, -1)
  sc <- rbind(q, cm)
  spec <- data.frame(variable = "v", transform = "dlog")
  base <- list(min_train_origins = 3, forgetting = 1, shrink_kappa = 0)
  cfg_q <- list(combination = base)
  cfg_l <- list(combination = c(base, list(train_measure = "level")))
  w_q <- combo_weights("logscore", sc, "v", 1:4, t = 25, members, cfg_q, spec)
  w_l <- combo_weights("logscore", sc, "v", 1:4, t = 25, members, cfg_l, spec)
  expect_gt(w_q["A"], 0.99)
  expect_gt(w_l["B"], 0.95)   # cum dominates (only h=1 rows come from q)
})

# ---- unrestricted member (D19) --------------------------------------------------

test_that("block_exog: false relaxes the restriction and is gate-exempt", {
  cfg <- tiny_cfg()
  spec <- build_transform_spec(cfg)
  raw <- generate_synthetic_data(cfg, spec)
  td <- transform_data(raw, spec)
  spec_s <- vars_for_set(spec, "small")
  y <- as.matrix(td[, spec_s$variable])
  member <- list(name = "u", kind = "var", engine = "gibbs", set = "small",
                 lags = 2, block_exog = FALSE,
                 prior = list(lambda = 0.2, soc = FALSE, dio = FALSE))
  set.seed(4)
  post <- fit_var_member(y, member, spec_s, cfg)
  expect_true(post$diagnostics$block_exog_exempt)
  # prior no longer pins the coefficients to ~0 (they are just shrunk)
  expect_gt(post$diagnostics$block_exog_max, 1e-3)
  # the gate accepts an exempt violation but rejects a non-exempt one
  dtab <- data.frame(block_exog_max = c(0.2, 1e-6),
                     block_exog_exempt = c(TRUE, FALSE),
                     converged_all = TRUE, sanity_all = TRUE,
                     no_lookahead = TRUE, reproducible = TRUE)
  expect_true(assert_diagnostics(dtab))
  dtab$block_exog_exempt <- FALSE
  expect_error(assert_diagnostics(dtab), "block exogeneity")
})

# ---- ar4 delta prior -------------------------------------------------------------

test_that("fit_ar centres the lag-1 prior on delta", {
  set.seed(5)
  # persistent level series, tiny sample: prior mean matters
  z <- as.numeric(filter(rnorm(28), 0.95, method = "recursive")) + 5
  y <- matrix(z, ncol = 1, dimnames = list(NULL, "u"))
  f0 <- fit_ar(y, list(), p = 4, delta = 0)
  f1 <- fit_ar(y, list(), p = 4, delta = 1)
  expect_gt(f1$fits[[1]]$bhat[2], f0$fits[[1]]$bhat[2])
})

# ---- PIT moment tests ------------------------------------------------------------

test_that("pit_moment_tests flags biased PITs and passes uniform ones", {
  set.seed(6)
  n <- 60
  mk <- function(member, pit) data.frame(
    member = member, variable = "v", measure = "q", h = 2,
    origin = seq_len(n), pit = pit)
  sc <- rbind(mk("good", runif(n)), mk("biased", rbeta(n, 6, 2)))
  pt <- pit_moment_tests(sc)
  expect_lt(pt$p_location[pt$member == "biased"], 0.01)
  expect_gt(pt$p_location[pt$member == "good"], 0.05)
})

# ---- event probabilities -----------------------------------------------------------

test_that("event probabilities are valid and respond to the draws", {
  cfg <- tiny_cfg()
  H <- cfg$horizons
  td <- data.frame(date = seq(as.Date("2000-01-01"), by = "quarter", length.out = 30),
                   gdp_growth = rnorm(30, 0.7, 0.5),
                   cpi_inflation = rnorm(30, 0.6, 0.2),
                   unemp_rate = rnorm(30, 4.5, 0.2))
  mk_dr <- function(mu) {
    a <- array(NA_real_, c(200, H, 3),
               dimnames = list(NULL, NULL, c("gdp_growth", "cpi_inflation", "unemp_rate")))
    a[, , "gdp_growth"] <- rnorm(200 * H, mu, 0.5)
    a[, , "cpi_inflation"] <- rnorm(200 * H, 0.625, 0.15)
    a[, , "unemp_rate"] <- rnorm(200 * H, 4.6, 0.3)
    a
  }
  wtab <- expand.grid(scheme = "equal", variable = c("gdp_growth", "cpi_inflation", "unemp_rate"),
                      bucket = c("near", "medium", "far"), member = c("m1", "m2"),
                      stringsAsFactors = FALSE)
  wtab$weight <- 0.5
  ff_hi <- list(member_draws = list(m1 = mk_dr(0.7), m2 = mk_dr(0.7)), weights = wtab)
  ff_lo <- list(member_draws = list(m1 = mk_dr(-2.0), m2 = mk_dr(-2.0)), weights = wtab)
  ep_hi <- event_probabilities(ff_hi, td, NULL, cfg, out_dir = tempdir())
  ep_lo <- event_probabilities(ff_lo, td, NULL, cfg, out_dir = tempdir())
  expect_true(all(ep_hi$prob >= 0 & ep_hi$prob <= 1))
  rec_hi <- ep_hi$prob[ep_hi$variable == "gdp_growth" & ep_hi$h == 8]
  rec_lo <- ep_lo$prob[ep_lo$variable == "gdp_growth" & ep_lo$h == 8]
  expect_lt(rec_hi, 0.3)   # healthy growth: low contraction probability
  expect_gt(rec_lo, 0.9)   # deep contraction draws: near-certain
})
