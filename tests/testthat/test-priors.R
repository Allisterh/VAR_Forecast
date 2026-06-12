source_project <- function() {
  root <- testthat::test_path("..", "..")
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
    sys.source(f, envir = globalenv())
}
source_project()

test_that("Minnesota prior has the documented structure", {
  M <- 4; p <- 2
  sigma <- c(1, 2, 0.5, 1)
  delta <- c(0, 1, 1, 0)
  blocks <- c("foreign", "foreign", "domestic", "domestic")
  pr <- minnesota_prior(M, p, sigma, delta, blocks, lambda = 0.2, cross = 0.5)
  expect_equal(dim(pr$b0), c(1 + M * p, M))
  # own first lag mean = delta
  for (i in seq_len(M)) expect_equal(pr$b0[1 + i, i], delta[i])
  # own-lag sd = lambda, cross-lag scaled by cross * sigma ratio
  expect_equal(pr$s0[1 + 1, 1], 0.2)
  expect_equal(pr$s0[1 + 2, 1], 0.2 * 0.5 * sigma[1] / sigma[2])
  # lag decay: lag-2 own sd = lambda / 2
  expect_equal(pr$s0[1 + M + 1, 1], 0.2 / 2)
})

test_that("block exogeneity restriction zeroes domestic lags in foreign equations", {
  M <- 4; p <- 2
  blocks <- c("foreign", "foreign", "domestic", "domestic")
  pr <- minnesota_prior(M, p, rep(1, M), rep(0, M), blocks,
                        lambda = 0.2, block_exog_sd = 1e-6)
  idx <- coef_index(M, p)
  dom_rows <- idx$row[idx$var %in% 3:4]
  expect_true(all(pr$s0[dom_rows, 1:2] == 1e-6))
  expect_true(all(pr$b0[dom_rows, 1:2] == 0))
  # domestic equations are NOT restricted
  expect_true(all(pr$s0[dom_rows, 3:4] > 1e-3))
})

test_that("inverse-Wishart sampler has the right mean", {
  set.seed(42)
  S <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  nu <- 12
  draws <- replicate(4000, riwish(nu, S))
  m <- apply(draws, c(1, 2), mean)
  expect_equal(m, S / (nu - 2 - 1), tolerance = 0.08)
})

test_that("SOC and DIO dummies have documented shapes and scaling", {
  M <- 3; p <- 2
  ybar <- c(1, 2, 3); delta <- c(1, 1, 0)
  d <- soc_dio_dummies(M, p, ybar, delta, soc = TRUE, soc_mu = 2,
                       dio = TRUE, dio_delta = 4)
  # SOC: one row per delta==1 variable; DIO: one row
  expect_equal(nrow(d$Y), 2 + 1)
  expect_equal(d$Y[1, 1], ybar[1] / 2)
  expect_equal(d$X[1, 1], 0)                      # no intercept in SOC rows
  expect_equal(d$X[3, 1], 1 / 4)                  # DIO intercept
  expect_equal(d$Y[3, ], ybar / 4)
})

test_that("conjugate marginal likelihood is finite and lambda selection works", {
  set.seed(7)
  y <- matrix(rnorm(200), 100, 2)
  y[, 1] <- filter(y[, 1], 0.5, method = "recursive")
  lam <- select_lambda(y, p = 2, sigma = ar_sigmas(y), delta = c(0, 0),
                       grid = c(0.05, 0.2, 0.6))
  expect_true(lam %in% c(0.05, 0.2, 0.6))
  pr <- conjugate_prior(2, 2, ar_sigmas(y), c(0, 0), lambda = 0.2)
  xy <- build_XY(y, 2)
  ml <- log_marginal_likelihood(xy$Y, xy$X, pr)
  expect_true(is.finite(ml))
})

test_that("steady-state mean update collapses to the analytic posterior", {
  # With A = 0 and Sigma = I, the model is y_t = Psi + e_t and the posterior
  # for Psi is the textbook normal-normal update.
  set.seed(11)
  M <- 2; p <- 1; Tn <- 60
  y <- matrix(rnorm(Tn * M, mean = c(1, -1)), Tn, M, byrow = FALSE)
  y[, 1] <- rnorm(Tn, 1); y[, 2] <- rnorm(Tn, -1)
  spec_m <- data.frame(variable = c("a", "b"), block = c("foreign", "domestic"),
                       delta = c(0, 0), ss_mean = c(0, 0), ss_sd = c(1, 1),
                       transform = "level")
  cfg <- list(mcmc = list(ndraw = 3000, nburn = 300, block_exog_prior_sd = 1e-4),
              master_seed = 1)
  member <- list(name = "t", engine = "ss", set = "x", lags = 1,
                 prior = list(lambda = 1e-6))   # lambda ~ 0 pins A at 0
  post <- fit_ss(y, member, spec_m, cfg,
                 prior = list(lambda = 1e-6, soc = FALSE, dio = FALSE))
  # analytic posterior: precision = 1/ss_sd^2 + (T-p)/sigma2_hat
  Te <- Tn - 1
  for (j in 1:2) {
    s2 <- var(y[(2):Tn, j])
    prec <- 1 + Te / s2
    mean_an <- (0 * 1 + sum(y[2:Tn, j]) / s2) / prec
    expect_equal(mean(post$Psi[, j]), mean_an, tolerance = 0.1)
    expect_equal(sd(post$Psi[, j]), sqrt(1 / prec), tolerance = 0.06)
  }
})
