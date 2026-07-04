# forecast.R -- iterated predictive simulation for every engine ----------------
#
# simulate_paths(post, y, h, ndraw, condition = NULL) -> array [ndraw, h, M]
# of future values in model units, integrating over parameter draws (recycled
# if ndraw > stored draws) and future shocks.
#
# condition: optional list(variable =, path =, method = "wz"|"substitute").
# Conditioning on a future path of one variable (e.g. a market-implied cash
# rate). path may be shorter than h and may contain NAs (those steps are
# unconstrained). Two methods (README.md D21):
#  * "wz" (default, Gaussian engines gibbs/ss/conj_br): Waggoner-Zha (1999)
#    conditional simulation -- for each posterior draw the joint future
#    shocks are drawn from their EXACT Gaussian conditional distribution
#    given the constraints, so non-conditioned variables respond through
#    both the lag dynamics and the error covariance and the conditional
#    predictive density is exact given (B, Sigma).
#  * "substitute": overwrite the conditioned variable each step before it
#    feeds the recursion. The only available method for the SV engine
#    (non-Gaussian shocks) and the univariate benchmarks; a documented
#    approximation (it ignores the contemporaneous shock correlation).

simulate_paths <- function(post, y, h, ndraw, condition = NULL,
                           shock_scale = NULL) {
  UseMethod("simulate_paths")
}

#' Normalise a condition: pad/trim path to h (NA = unconstrained step),
#' resolve the default method.
.norm_condition <- function(condition, h) {
  if (is.null(condition)) return(NULL)
  p <- as.numeric(condition$path)
  condition$path <- c(p, rep(NA_real_, max(0, h - length(p))))[seq_len(h)]
  if (is.null(condition$method)) condition$method <- "wz"
  condition
}

.apply_condition <- function(ystep, s, condition, varnames) {
  if (is.null(condition)) return(ystep)
  j <- match(condition$variable, varnames)
  if (!is.na(j) && s <= length(condition$path) &&
      is.finite(condition$path[s])) ystep[j] <- condition$path[s]
  ystep
}

# ---- Waggoner-Zha conditional machinery -----------------------------------------

#' MA weights Psi_0..Psi_{H-1} (list of M x M) of the VAR y_t = c + A(L)y + u:
#' Psi_0 = I, Psi_j = sum_{l=1..min(j,p)} A_l Psi_{j-l}. B is K x M with the
#' intercept in row 1 and lag blocks below (the build_XY convention).
.ma_weights <- function(B, M, p, H) {
  A <- lapply(seq_len(p), function(l)
    t(B[1 + ((l - 1) * M + 1):(l * M), , drop = FALSE]))
  Psi <- vector("list", H)
  Psi[[1]] <- diag(M)
  if (H > 1) for (k in 2:H) {
    S <- matrix(0, M, M)
    for (l in seq_len(min(k - 1, p))) S <- S + A[[l]] %*% Psi[[k - l]]
    Psi[[k]] <- S
  }
  Psi
}

#' Deterministic (zero-shock) path from state ystate under coefficients B.
.deterministic_path <- function(B, ystate, h, M, p) {
  out <- matrix(NA_real_, h, M)
  st <- ystate
  for (s in seq_len(h)) {
    x <- c(1, as.vector(t(st[seq_len(p), , drop = FALSE])))
    out[s, ] <- drop(crossprod(B, x))
    st <- rbind(out[s, ], st[-p, , drop = FALSE])
  }
  out
}

#' One Waggoner-Zha conditional path draw for a Gaussian VAR (B, Sigma).
#' cj: column index of the conditioned variable; path: length-h with NAs for
#' unconstrained steps; ss_: per-step shock scale (COVID decay path).
#' The future stacked shocks u ~ N(0, blockdiag(ss_s^2 Sigma)) are drawn from
#' their exact conditional distribution given R u = r via
#'   u* = u + C R' (R C R')^{-1} (r - R u),
#' then the path is iterated forward with u*.
.wz_conditional_path <- function(B, Sigma, cS, ystate, h, cj, path, M, p, ss_) {
  m <- .deterministic_path(B, ystate, h, M, p)
  cons <- which(is.finite(path))
  Psi <- .ma_weights(B, M, p, h)
  # R: (#constraints x M*h); constraint s row: y_{cj,s} - m_{cj,s} =
  #   sum_{j<=s} [Psi_{s-j}]_{cj,.} u_j
  R <- matrix(0, length(cons), M * h)
  for (i in seq_along(cons)) {
    s <- cons[i]
    for (j in seq_len(s))
      R[i, ((j - 1) * M + 1):(j * M)] <- Psi[[s - j + 1]][cj, ]
  }
  r <- path[cons] - m[cons, cj]
  # unconditional shock draw and blockdiag covariance
  u <- as.vector(vapply(seq_len(h), function(s)
    ss_[s] * drop(crossprod(cS, rnorm(M))), numeric(M)))
  C_Rt <- matrix(0, M * h, length(cons))          # C R', C = blockdiag(ss^2 Sigma)
  for (j in seq_len(h)) {
    ix <- ((j - 1) * M + 1):(j * M)
    C_Rt[ix, ] <- (ss_[j]^2) * (Sigma %*% t(R[, ix, drop = FALSE]))
  }
  RCRt <- R %*% C_Rt
  adj <- C_Rt %*% solve(RCRt, r - drop(R %*% u))
  ustar <- u + drop(adj)
  # iterate forward with the conditional shocks
  out <- matrix(NA_real_, h, M)
  st <- ystate
  for (s in seq_len(h)) {
    x <- c(1, as.vector(t(st[seq_len(p), , drop = FALSE])))
    out[s, ] <- drop(crossprod(B, x)) + ustar[((s - 1) * M + 1):(s * M)]
    st <- rbind(out[s, ], st[-p, , drop = FALSE])
  }
  out
}

#' TRUE when a condition should take the WZ route for this posterior.
.use_wz <- function(condition, varnames) {
  !is.null(condition) && identical(condition$method, "wz") &&
    condition$variable %in% varnames && any(is.finite(condition$path))
}

#' Iterate one VAR path: B (K x M, intercept first), cS = chol(Sigma),
#' ystate = matrix p x M with row 1 = most recent observation.
.iterate_path <- function(B, cS, ystate, h, M, p, condition, varnames,
                          shock_scale = NULL) {
  out <- matrix(NA_real_, h, M)
  ss_ <- if (is.null(shock_scale)) rep(1, h) else shock_scale
  for (s in seq_len(h)) {
    x <- c(1, as.vector(t(ystate[seq_len(p), , drop = FALSE])))
    mu_s <- drop(crossprod(B, x))
    ynew <- mu_s + ss_[s] * drop(crossprod(cS, rnorm(M)))
    ynew <- .apply_condition(ynew, s, condition, varnames)
    out[s, ] <- ynew
    ystate <- rbind(ynew, ystate[-p, , drop = FALSE])
  }
  out
}

# state helper: last p observations, most recent first
.ystate <- function(y, p) y[nrow(y) - seq_len(p) + 1, , drop = FALSE]

# ---- gibbs ----------------------------------------------------------------------

simulate_paths.post_gibbs <- function(post, y, h, ndraw, condition = NULL,
                                      shock_scale = NULL) {
  M <- post$M; p <- post$p
  condition <- .norm_condition(condition, h)
  wz <- .use_wz(condition, post$varnames)
  cj <- if (wz) match(condition$variable, post$varnames) else NA_integer_
  ss_ <- if (is.null(shock_scale)) rep(1, h) else shock_scale
  paths <- array(NA_real_, c(ndraw, h, M))
  st0 <- .ystate(y, p)
  for (d in seq_len(ndraw)) {
    k <- floor((d - 1) * post$ndraw / ndraw) + 1   # seed paths across the full posterior
    B <- post$B[k, , ]
    Sig <- post$Sigma[k, , ]
    cS <- chol(Sig)
    paths[d, , ] <- if (wz) {
      .wz_conditional_path(B, Sig, cS, st0, h, cj, condition$path, M, p, ss_)
    } else {
      .iterate_path(B, cS, st0, h, M, p, condition, post$varnames,
                    shock_scale = shock_scale)
    }
  }
  dimnames(paths)[[3]] <- post$varnames
  paths
}

# ---- ss -------------------------------------------------------------------------

simulate_paths.post_ss <- function(post, y, h, ndraw, condition = NULL,
                                   shock_scale = NULL) {
  M <- post$M; p <- post$p
  condition <- .norm_condition(condition, h)
  wz <- .use_wz(condition, post$varnames)
  cj <- if (wz) match(condition$variable, post$varnames) else NA_integer_
  ss_ <- if (is.null(shock_scale)) rep(1, h) else shock_scale
  paths <- array(NA_real_, c(ndraw, h, M))
  for (d in seq_len(ndraw)) {
    k <- floor((d - 1) * post$ndraw / ndraw) + 1   # seed paths across the full posterior
    A <- post$A[k, , ]; Psi <- post$Psi[k, ]
    Sig <- post$Sigma[k, , ]
    cS <- chol(Sig)
    z <- sweep(y, 2, Psi)
    st <- .ystate(z, p)
    B <- rbind(0, A)                      # zero intercept in demeaned form
    # condition in DEMEANED units so the conditioned variable feeds the
    # recursion each step (same contract as the other engines)
    cond_z <- condition
    if (!is.null(cond_z)) {
      j <- match(cond_z$variable, post$varnames)
      if (!is.na(j)) cond_z$path <- cond_z$path - Psi[j]
    }
    zp <- if (wz) {
      .wz_conditional_path(B, Sig, cS, st, h, cj, cond_z$path, M, p, ss_)
    } else {
      .iterate_path(B, cS, st, h, M, p, condition = cond_z, post$varnames,
                    shock_scale = shock_scale)
    }
    paths[d, , ] <- sweep(zp, 2, Psi, "+")
  }
  dimnames(paths)[[3]] <- post$varnames
  paths
}

# ---- conj_br ---------------------------------------------------------------------

#' Joint reduced form (B_joint K x M, Sigma_joint M x M) of a block-recursive
#' draw: substitute the foreign equations into the domestic conditional.
#' With y_f = Bf'xf + e_f and y_d = b0d + Blag'x + G'y_f + e_d (G the
#' contemporaneous-foreign coefficient rows of Bd):
#'   u_f = e_f,  u_d = G'e_f + e_d,
#'   Sigma_joint = [[Sf, Sf G], [G'Sf, G'Sf G + Sd]].
.conj_br_joint <- function(Bf, Sf, Bd, Sd, M, p, nf) {
  nd <- M - nf
  K <- 1 + M * p
  G <- Bd[(1 + M * p + 1):(1 + M * p + nf), , drop = FALSE]  # nf x nd
  Bj <- matrix(0, K, M)
  # foreign equations: intercept + foreign-lag coefficients in joint rows
  Bj[1, seq_len(nf)] <- Bf[1, ]
  for (l in seq_len(p)) for (j in seq_len(nf))
    Bj[1 + (l - 1) * M + j, seq_len(nf)] <- Bf[1 + (l - 1) * nf + j, ]
  # domestic equations: own lag block + foreign equations routed through G
  Bj[seq_len(K), nf + seq_len(nd)] <- Bd[seq_len(K), , drop = FALSE]
  Bj[, nf + seq_len(nd)] <- Bj[, nf + seq_len(nd), drop = FALSE] +
    Bj[, seq_len(nf), drop = FALSE] %*% G
  SfG <- Sf %*% G
  Sj <- rbind(cbind(Sf, SfG), cbind(t(SfG), t(G) %*% SfG + Sd))
  list(B = Bj, Sigma = (Sj + t(Sj)) / 2)
}

simulate_paths.post_conj_br <- function(post, y, h, ndraw, condition = NULL,
                                        shock_scale = NULL) {
  M <- post$M; p <- post$p; nf <- post$nf
  nd <- M - nf
  condition <- .norm_condition(condition, h)
  wz <- .use_wz(condition, post$varnames)
  cj <- if (wz) match(condition$variable, post$varnames) else NA_integer_
  ss_ <- if (is.null(shock_scale)) rep(1, h) else shock_scale
  paths <- array(NA_real_, c(ndraw, h, M))
  st0 <- .ystate(y, p)                       # all vars, most recent first
  for (d in seq_len(ndraw)) {
    k <- floor((d - 1) * post$ndraw / ndraw) + 1   # seed paths across the full posterior
    Bf <- post$foreign$B[k, , ];  cSf <- chol(post$foreign$Sigma[k, , ])
    Bd <- post$domestic$B[k, , ]; cSd <- chol(post$domestic$Sigma[k, , ])
    if (wz) {
      jr <- .conj_br_joint(Bf, post$foreign$Sigma[k, , ], Bd,
                           post$domestic$Sigma[k, , ], M, p, nf)
      paths[d, , ] <- .wz_conditional_path(jr$B, jr$Sigma, chol(jr$Sigma),
                                           st0, h, cj, condition$path, M, p, ss_)
      next
    }
    st <- st0
    for (s in seq_len(h)) {
      xf <- c(1, as.vector(t(st[seq_len(p), seq_len(nf), drop = FALSE])))
      yf <- drop(crossprod(Bf, xf)) + ss_[s] * drop(crossprod(cSf, rnorm(nf)))
      xd <- c(1, as.vector(t(st[seq_len(p), , drop = FALSE])), yf)
      yd <- drop(crossprod(Bd, xd)) + ss_[s] * drop(crossprod(cSd, rnorm(nd)))
      ynew <- c(yf, yd)
      ynew <- .apply_condition(ynew, s, condition, post$varnames)
      paths[d, s, ] <- ynew
      st <- rbind(ynew, st[-p, , drop = FALSE])
    }
  }
  dimnames(paths)[[3]] <- post$varnames
  paths
}

# ---- sv --------------------------------------------------------------------------

simulate_paths.post_sv <- function(post, y, h, ndraw, condition = NULL,
                                   shock_scale = NULL) {
  # shock_scale ignored: the SV engine handles COVID via t-errors + SV
  M <- post$M; p <- post$p
  condition <- .norm_condition(condition, h)
  if (.use_wz(condition, post$varnames)) {
    # WZ requires jointly Gaussian future shocks; the SV engine's are a
    # t/scale mixture with a stochastic variance path, so fall back to
    # substitution (documented approximation, README.md D21).
    log_warn("sv engine: conditional method 'wz' unavailable, using substitution")
    condition$method <- "substitute"
  }
  paths <- array(NA_real_, c(ndraw, h, M))
  st0 <- .ystate(y, p)
  for (d in seq_len(ndraw)) {
    k <- floor((d - 1) * post$ndraw / ndraw) + 1   # seed paths across the full posterior
    st <- st0
    # pull per-equation parameters for this draw
    betas <- lapply(post$eqs, function(e) e$beta[k, ])
    hs    <- vapply(post$eqs, function(e) e$hT[k], numeric(1))
    svp   <- lapply(post$eqs, function(e) e$svpara[k, ])
    for (s in seq_len(h)) {
      ynew <- numeric(M)
      lags_vec <- as.vector(t(st[seq_len(p), , drop = FALSE]))
      for (i in seq_len(M)) {
        eq <- post$eqs[[i]]$design
        # rebuild the design vector for this equation: 1, allowed lags, contemp
        allow <- eq$lag_meta
        xlag <- lags_vec[(allow$lag - 1) * M + allow$var]
        xx <- c(1, xlag, if (length(eq$contemp)) ynew[eq$contemp])
        hs[i] <- svp[[i]]["mu"] + svp[[i]]["phi"] * (hs[i] - svp[[i]]["mu"]) +
                 svp[[i]]["sigma"] * rnorm(1)
        nu_i <- svp[[i]]["nu"]
        # standardise the t-draw to UNIT variance: estimation calibrates exp(h)
        # to be the conditional variance, but Var[rt(nu)] = nu/(nu-2) > 1, so a
        # raw rt() would over-disperse every SV predictive density.
        eps <- if (!is.na(nu_i) && is.finite(nu_i) && nu_i > 2)
          rt(1, df = nu_i) * sqrt((nu_i - 2) / nu_i) else rnorm(1)
        ynew[i] <- sum(xx * betas[[i]]) + exp(hs[i] / 2) * eps
      }
      ynew <- .apply_condition(ynew, s, condition, post$varnames)
      paths[d, s, ] <- ynew
      st <- rbind(ynew, st[-p, , drop = FALSE])
    }
  }
  dimnames(paths)[[3]] <- post$varnames
  paths
}

# ---- fan-chart quantiles + sanity checks ------------------------------------------

#' Quantiles per horizon/variable from a path array.
fan_quantiles <- function(paths, probs) {
  qs <- apply(paths, c(2, 3), quantile, probs = probs)
  # qs: [prob, h, var] -> long data.frame
  out <- expand.grid(prob = probs, h = seq_len(dim(paths)[2]),
                     variable = dimnames(paths)[[3]], stringsAsFactors = FALSE)
  out$value <- as.vector(qs)
  out
}

#' Forecast sanity (section 9): finite, non-explosive in the predictive
#' tails, and non-explosive in the median. The "converged" test requires the
#' MEDIAN path to stay within a generous band around the historical data range
#' (range +/- 4 sd) at every horizon. This is the robust, Monte-Carlo-stable
#' way to express "the iterated forecast settles rather than diverging": a
#' stable VAR's median converges to its model-implied unconditional mean (which
#' legitimately differs from the raw sample mean -- e.g. forecasting from the
#' 2020Q2 trough, or with COVID observations downweighted), and a random walk's
#' median is flat; both stay in band, while a diverging/explosive median leaves
#' it. It deliberately does NOT use increment-to-increment comparisons (too
#' sensitive to median MC noise under COVID variance inflation) nor reversion to
#' the sample mean (the wrong target for a downweighted/extreme-origin model).
check_forecasts <- function(paths, y, label = "member", delta = NULL) {
  ok_finite <- all(is.finite(paths))
  rng <- apply(abs(y), 2, max)
  ok_bound <- TRUE
  mean_path <- apply(paths, c(2, 3), median)
  for (j in seq_len(ncol(y))) {
    bound <- 5 * max(rng[j], 1) + 50
    # explosiveness is about the central mass, not individual fat-tail draws
    ql <- apply(paths[, , j, drop = FALSE], 2, quantile, probs = c(0.005, 0.995))
    if (any(abs(ql) > bound)) ok_bound <- FALSE
  }
  sdv <- apply(y, 2, sd)
  lo <- apply(y, 2, min) - 4 * sdv
  hi <- apply(y, 2, max) + 4 * sdv
  rev_ok <- vapply(seq_len(ncol(y)), function(j)
    all(mean_path[, j] >= lo[j] & mean_path[, j] <= hi[j]), logical(1))
  if (!is.null(delta)) rev_ok <- rev_ok | (delta == 1)
  ok_converge <- all(rev_ok)
  ok <- ok_finite && ok_bound && ok_converge
  if (!ok) log_warn(
    "forecast sanity FAILED for {label}: finite={ok_finite} bounded={ok_bound} converge={ok_converge}")
  list(finite = ok_finite, bounded = ok_bound, converged = ok_converge, ok = ok)
}
