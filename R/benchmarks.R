# benchmarks.R -- univariate benchmark members: rw, ar4, ucsv, ucmean ----------
#
# Benchmarks implement the same interface as VAR members: fit_benchmark()
# returns an object with a simulate_paths() method producing [ndraw, h, nvar]
# arrays over the *target* variables. They are pool members and the bar every
# VAR must clear.

#' Random walk on the modelled (transformed) series with Gaussian increments.
fit_rw <- function(y, cfg) {
  structure(list(engine = "rw", varnames = colnames(y),
                 yT = y[nrow(y), ], sdd = apply(diff(y), 2, sd),
                 diagnostics = list(converged = TRUE, ess_min = Inf,
                                    stable_share = 1, block_exog_max = 0)),
            class = c("post_rw", "var_posterior"))
}

simulate_paths.post_rw <- function(post, y, h, ndraw, condition = NULL) {
  nv <- length(post$yT)
  paths <- array(NA_real_, c(ndraw, h, nv), dimnames = list(NULL, NULL, post$varnames))
  for (j in seq_len(nv)) {
    inc <- matrix(rnorm(ndraw * h, 0, post$sdd[j]), ndraw, h)
    paths[, , j] <- post$yT[j] + t(apply(inc, 1, cumsum))
  }
  paths
}

#' Bayesian AR(p) per variable: conjugate normal-inverse-gamma with a loose
#' ridge prior; iterated density forecasts.
fit_ar <- function(y, cfg, p = 4) {
  fits <- lapply(seq_len(ncol(y)), function(j) {
    z <- y[, j]
    xy <- build_XY(matrix(z, ncol = 1), p)
    X <- xy$X; Y <- drop(xy$Y)
    K <- ncol(X)
    V0inv <- diag(c(1e-4, rep(1, p)))         # loose on intercept, unit ridge on lags
    P <- crossprod(X) + V0inv
    cP <- chol(P)
    bhat <- backsolve(cP, forwardsolve(t(cP), crossprod(X, Y)))
    resid <- Y - drop(X %*% bhat)
    s2 <- sum(resid^2) / (length(Y) - K)
    list(bhat = bhat, cP = cP, s2 = s2, df = length(Y) - K, p = p)
  })
  structure(list(engine = "ar4", varnames = colnames(y), fits = fits, p = p,
                 diagnostics = list(converged = TRUE, ess_min = Inf,
                                    stable_share = 1, block_exog_max = 0)),
            class = c("post_ar", "var_posterior"))
}

simulate_paths.post_ar <- function(post, y, h, ndraw, condition = NULL) {
  nv <- length(post$fits); p <- post$p
  paths <- array(NA_real_, c(ndraw, h, nv), dimnames = list(NULL, NULL, post$varnames))
  for (j in seq_len(nv)) {
    f <- post$fits[[j]]
    z <- y[, j]
    for (d in seq_len(ndraw)) {
      sig2 <- f$s2 * f$df / rchisq(1, f$df)
      beta <- f$bhat + sqrt(sig2) * backsolve(f$cP, rnorm(p + 1))
      st <- z[length(z) - seq_len(p) + 1]      # most recent first
      for (s in seq_len(h)) {
        ynew <- sum(c(1, st) * beta) + sqrt(sig2) * rnorm(1)
        paths[d, s, j] <- ynew
        st <- c(ynew, st[-p])
      }
    }
  }
  paths
}

# ---- UCSV (Stock-Watson unobserved components with twin SV) --------------------

#' FFBS for the local level model with time-varying variances:
#' y_t = tau_t + e_t, e_t ~ N(0, s2e_t); tau_t = tau_{t-1} + u_t, u_t ~ N(0, s2u_t).
ffbs_local_level <- function(y, s2e, s2u, tau0 = y[1], P0 = 10 * var(y)) {
  Tn <- length(y)
  af <- numeric(Tn); Pf <- numeric(Tn)
  a <- tau0; P <- P0
  for (t in seq_len(Tn)) {
    P <- P + s2u[t]
    Kk <- P / (P + s2e[t])
    a <- a + Kk * (y[t] - a)
    P <- (1 - Kk) * P
    af[t] <- a; Pf[t] <- P
  }
  tau <- numeric(Tn)
  tau[Tn] <- rnorm(1, af[Tn], sqrt(Pf[Tn]))
  for (t in (Tn - 1):1) {
    Pp <- Pf[t] + s2u[t + 1]
    J <- Pf[t] / Pp
    m <- af[t] + J * (tau[t + 1] - af[t])
    V <- Pf[t] * (1 - J)
    tau[t] <- rnorm(1, m, sqrt(max(V, 1e-12)))
  }
  tau
}

#' UCSV per variable via Gibbs: FFBS for the trend, stochvol updates for both
#' log-variance processes.
fit_ucsv <- function(y, cfg) {
  ndraw <- cfg$mcmc$ndraw; nburn <- cfg$mcmc$nburn
  pspec <- stochvol::specify_priors(
    mu = stochvol::sv_normal(0, 10),
    phi = stochvol::sv_beta(20, 1.5),
    sigma2 = stochvol::sv_gamma(0.5, 0.5))
  fits <- lapply(seq_len(ncol(y)), function(j) {
    z <- y[, j]; Tn <- length(z)
    tau <- stats::filter(z, rep(1 / 8, 8), sides = 1)
    tau[is.na(tau)] <- z[is.na(tau)]
    tau <- as.numeric(tau)
    he <- rep(log(var(z) / 2 + 1e-8), Tn); hu <- rep(log(var(z) / 20 + 1e-8), Tn)
    pe <- list(mu = he[1], phi = 0.9, sigma = 0.2, nu = Inf, rho = 0, beta = 0,
               latent0 = he[1])
    pu <- list(mu = hu[1], phi = 0.9, sigma = 0.2, nu = Inf, rho = 0, beta = 0,
               latent0 = hu[1])
    tauT <- numeric(ndraw); heT <- numeric(ndraw); huT <- numeric(ndraw)
    pe_d <- matrix(NA_real_, ndraw, 3); pu_d <- matrix(NA_real_, ndraw, 3)
    for (it in seq_len(nburn + ndraw)) {
      tau <- ffbs_local_level(z, exp(he), exp(hu))
      eres <- z - tau
      ures <- c(diff(tau)[1], diff(tau))   # u_1 approximated by u_2
      ue <- stochvol::svsample_fast_cpp(eres, draws = 1, burnin = 0,
        priorspec = pspec, startpara = pe, startlatent = he)
      he <- drop(ue$latent[1, ])
      pe$mu <- ue$para[1, "mu"]; pe$phi <- ue$para[1, "phi"]
      pe$sigma <- ue$para[1, "sigma"]; pe$latent0 <- he[1]
      uu <- stochvol::svsample_fast_cpp(ures, draws = 1, burnin = 0,
        priorspec = pspec, startpara = pu, startlatent = hu)
      hu <- drop(uu$latent[1, ])
      pu$mu <- uu$para[1, "mu"]; pu$phi <- uu$para[1, "phi"]
      pu$sigma <- uu$para[1, "sigma"]; pu$latent0 <- hu[1]
      if (it > nburn) {
        d <- it - nburn
        tauT[d] <- tau[Tn]; heT[d] <- he[Tn]; huT[d] <- hu[Tn]
        pe_d[d, ] <- c(pe$mu, pe$phi, pe$sigma)
        pu_d[d, ] <- c(pu$mu, pu$phi, pu$sigma)
      }
    }
    list(tauT = tauT, heT = heT, huT = huT, pe = pe_d, pu = pu_d)
  })
  ess <- min(vapply(fits, function(f) as.numeric(coda::effectiveSize(f$tauT)),
                    numeric(1)))
  structure(list(engine = "ucsv", varnames = colnames(y), fits = fits,
                 ndraw = ndraw,
                 diagnostics = list(converged = ess > 30, ess_min = ess,
                                    stable_share = 1, block_exog_max = 0)),
            class = c("post_ucsv", "var_posterior"))
}

simulate_paths.post_ucsv <- function(post, y, h, ndraw, condition = NULL) {
  nv <- length(post$fits)
  paths <- array(NA_real_, c(ndraw, h, nv), dimnames = list(NULL, NULL, post$varnames))
  for (j in seq_len(nv)) {
    f <- post$fits[[j]]
    for (d in seq_len(ndraw)) {
      k <- ((d - 1) %% post$ndraw) + 1
      tau <- f$tauT[k]; he <- f$heT[k]; hu <- f$huT[k]
      pe <- f$pe[k, ]; pu <- f$pu[k, ]
      for (s in seq_len(h)) {
        hu <- pu[1] + pu[2] * (hu - pu[1]) + pu[3] * rnorm(1)
        he <- pe[1] + pe[2] * (he - pe[1]) + pe[3] * rnorm(1)
        tau <- tau + exp(hu / 2) * rnorm(1)
        paths[d, s, j] <- tau + exp(he / 2) * rnorm(1)
      }
    }
  }
  paths
}

#' Unconditional mean with Gaussian predictive (expanding moments).
fit_ucmean <- function(y, cfg) {
  structure(list(engine = "ucmean", varnames = colnames(y),
                 mu = colMeans(y), sd = apply(y, 2, sd),
                 diagnostics = list(converged = TRUE, ess_min = Inf,
                                    stable_share = 1, block_exog_max = 0)),
            class = c("post_ucmean", "var_posterior"))
}

simulate_paths.post_ucmean <- function(post, y, h, ndraw, condition = NULL) {
  nv <- length(post$mu)
  paths <- array(rnorm(ndraw * h * nv), c(ndraw, h, nv),
                 dimnames = list(NULL, NULL, post$varnames))
  for (j in seq_len(nv)) paths[, , j] <- post$mu[j] + paths[, , j] * post$sd[j]
  paths
}

#' Dispatcher mirroring fit_var_member.
fit_benchmark <- function(y_targets, name, cfg) {
  switch(name,
         rw     = fit_rw(y_targets, cfg),
         ar4    = fit_ar(y_targets, cfg, p = 4),
         ucsv   = fit_ucsv(y_targets, cfg),
         ucmean = fit_ucmean(y_targets, cfg),
         stop("unknown benchmark: ", name))
}

#' Hook for externally supplied forecasts (e.g. published RBA forecasts):
#' if a CSV with columns (origin_date, variable, h, q05..q95 or point, sd)
#' exists at `path`, it is read and returned for comparison in the report.
read_external_forecasts <- function(path = "data/external_forecasts.csv") {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}
