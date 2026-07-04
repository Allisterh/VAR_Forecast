# evaluate.R -- pseudo-real-time recursive OOS evaluation -----------------------
#
# No-look-ahead is enforced structurally: forecast_at_origin() receives only
# td[1:t, ] and everything downstream (GLP lambda selection, estimation,
# prediction) sees that truncated panel alone. test_no_lookahead() verifies it.

# scoringRules is OPTIONAL: scoring goes through safe_crps()/safe_logs() from
# aaa_capabilities.R (NA when it is absent), so it is not attached here.
suppressPackageStartupMessages({ library(furrr) })

#' Forecast origins (row indices of the transformed panel).
oos_origins <- function(td, cfg) {
  T_n <- nrow(td)
  t0 <- ceiling(cfg$evaluation$first_origin_frac * T_n)
  origins <- t0:(T_n - 1)
  if (length(origins) > cfg$evaluation$max_origins)
    origins <- tail(origins, cfg$evaluation$max_origins)
  origins
}

#' All pool members: suite entries + benchmarks, in a common list format.
all_members <- function(cfg) {
  suite <- lapply(cfg$suite, function(m) { m$kind <- "var"; m })
  bench <- lapply(cfg$benchmarks, function(b)
    list(name = b, kind = "benchmark", engine = b))
  members <- c(suite, bench)
  if (!has_stochvol()) {   # drop SV members + ucsv when stochvol is unavailable
    drop <- vapply(members, function(m) m$engine %in% c("sv", "ucsv"), logical(1))
    if (any(drop)) {
      nm <- vapply(members[drop], `[[`, "", "name")
      log_warn("stochvol unavailable: dropping SV members {paste(nm, collapse = ', ')}")
      members <- members[!drop]
    }
  }
  members
}

#' Fit + forecast a single member at one origin using ONLY data up to t.
#' Returns thinned predictive draws [store_draws, h, n_targets] plus sanity info.
#' condition: optional conditioning spec passed through to simulate_paths
#' (list(variable=, path=, method=); README.md D21) -- used for scenario
#' forecasts at the final origin, never inside the OOS evaluation.
forecast_at_origin <- function(member, td_t, spec, cfg, condition = NULL) {
  H <- cfg$horizons
  tgt <- spec$variable[spec$target]
  nfd <- cfg$mcmc$forecast_draws
  if (member$kind == "var") {
    set_name <- member$set
    spec_m <- vars_for_set(spec, set_name)
    y <- as.matrix(td_t[, spec_m$variable])
    # COVID treatment (LP scaling / dummy): scales estimated from THIS
    # member's data up to the origin only -- no look-ahead. The SV engine
    # ignores the weights (t-errors instead).
    cov <- covid_treatment(y, td_t$date, member$lags, ar_sigmas(y),
                           spec_m$delta, cfg, H)
    if (!is.null(cov$weights) && member$engine != "sv") {
      sig_w <- ar_sigmas(y, weights = cov$weights)
      cov2 <- covid_treatment(y, td_t$date, member$lags, sig_w,
                              spec_m$delta, cfg, H)
      if (!is.null(cov2$weights)) cov <- cov2
    }
    use_w <- if (member$engine == "sv") NULL else cov$weights
    s_fut <- if (member$engine == "sv") NULL else cov$s_future
    if (identical(member$prior$lambda, "auto") && isTRUE(cfg$glp$enabled)) {
      glp_lambda <- select_lambda(y, member$lags,
                                  ar_sigmas(y, weights = use_w), spec_m$delta,
                                  grid = unlist(cfg$glp$lambda_grid),
                                  weights = use_w)
    } else glp_lambda <- 0.2
    post <- fit_var_member(y, member, spec_m, cfg, glp_lambda = glp_lambda,
                           weights = use_w)
    paths <- simulate_paths(post, y, H, nfd, condition = condition,
                            shock_scale = s_fut)
    keep <- intersect(tgt, dimnames(paths)[[3]])
    paths <- paths[, , keep, drop = FALSE]
    covid_info <- list(scales = cov$scales, rho = cov$rho)
  } else {
    cfgb <- cfg
    cfgb$mcmc$ndraw <- cfg$mcmc$bench_ndraw
    cfgb$mcmc$nburn <- cfg$mcmc$bench_nburn
    y <- as.matrix(td_t[, tgt])
    cov <- covid_treatment(y, td_t$date, 4, ar_sigmas(y),
                           spec$delta[match(tgt, spec$variable)], cfg, H)
    post <- fit_benchmark(y, member$engine, cfgb, weights = cov$weights,
                          delta = spec$delta[match(tgt, spec$variable)])
    paths <- simulate_paths(post, y, H, nfd, condition = condition,
                            shock_scale = cov$s_future)
    glp_lambda <- NA_real_
    covid_info <- list(scales = cov$scales, rho = cov$rho)
  }
  vnames <- dimnames(paths)[[3]]
  sanity <- check_forecasts(paths, as.matrix(td_t[, vnames]),
                            label = member$name,
                            delta = spec$delta[match(vnames, spec$variable)])
  thin <- round(seq(1, dim(paths)[1], length.out = cfg$mcmc$store_draws))
  list(draws = paths[thin, , , drop = FALSE],
       diagnostics = post$diagnostics, sanity = sanity,
       glp_lambda = glp_lambda, covid = covid_info)
}

#' The single harness entry point that slices the panel at an origin: this is
#' the ONLY place the estimation window is cut, so the no-look-ahead test can
#' exercise the same code path the evaluation uses. Honors the configured
#' window type (expanding or rolling).
harness_forecast <- function(member, td, t, spec, cfg) {
  first <- if (identical(cfg$evaluation$window, "rolling")) {
    max(1L, t - cfg$evaluation$rolling_length + 1L)
  } else 1L
  set.seed(derive_seed(cfg$master_seed, paste0(member$name, "-", t)))
  out <- forecast_at_origin(member, td[first:t, , drop = FALSE], spec, cfg)
  out$origin <- t
  out
}

#' Run the recursive loop for one member over all origins, with disk caching
#' keyed by the config hash. Parallel over origins.
run_oos_member <- function(member, td, spec, cfg, cache_root = "cache") {
  hash <- config_hash(cfg)
  cdir <- file.path(cache_root, paste0("oos_", hash))
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  origins <- oos_origins(td, cfg)
  res <- furrr::future_map(origins, function(t) {
    ensure_project_loaded()
    f <- file.path(cdir, sprintf("%s_o%03d.rds", member$name, t))
    if (file.exists(f)) return(readRDS(f))
    out <- harness_forecast(member, td, t, spec, cfg)
    saveRDS(out, f)
    out
  }, .options = furrr::furrr_options(seed = TRUE))
  names(res) <- as.character(origins)
  log_info("OOS done: {member$name} ({length(origins)} origins)")
  res
}

# ---- scoring --------------------------------------------------------------------

#' Score one member's OOS results against realizations. Returns a long
#' data.frame: origin, date, variable, measure, h, point, real, logdens, crps,
#' pit. Measures, per target variable:
#'   q   - the modelled variable at quarter t+h (quarterly growth for dlog
#'         variables; the level for level variables).
#'   ye  - year-ended: the 4-quarter sum ending at t+h (dlog targets only) --
#'         the RBA's headline concept (e.g. year-ended GDP growth / trimmed-mean
#'         inflation, the 2-3% target).
#'   cum - cumulative level from the origin: the h-quarter sum t+1..t+h (dlog
#'         variables only) = 100*(log level_{t+h} - log level_t), i.e. where the
#'         level lands h quarters out. Redundant with q at h=1; for level
#'         variables it equals a constant shift of q and is omitted.
score_member <- function(member_name, oos_res, td, spec, cfg) {
  T_n <- nrow(td)
  tgt <- spec$variable[spec$target]
  rows <- list()
  for (res in oos_res) {
    t <- res$origin
    dr <- res$draws                       # [ndraw, H, var]
    H <- dim(dr)[2]
    for (v in dimnames(dr)[[3]]) {
      vt <- spec[spec$variable == v, ]
      real_q <- td[[v]]
      for (h in seq_len(H)) {
        if (t + h > T_n) next
        x <- dr[, h, v]
        y_real <- real_q[t + h]
        rows[[length(rows) + 1]] <- data.frame(
          member = member_name, origin = t, date = td$date[t + h],
          variable = v, measure = "q", h = h,
          point = mean(x), real = y_real,
          logdens = -safe_logs(y_real, x),
          crps = safe_crps(y_real, x),
          pit = mean(x <= y_real))
        # year-ended: sum of 4 quarterly outcomes ending at t+h
        if (isTRUE(vt$year_ended) && vt$transform == "dlog") {
          k <- (h - 3):h
          hist_part <- sum(real_q[t + k[k <= 0]])
          fc_idx <- k[k >= 1]
          xye <- rowSums(dr[, fc_idx, v, drop = FALSE]) + hist_part
          ye_real <- sum(real_q[t + k])
          rows[[length(rows) + 1]] <- data.frame(
            member = member_name, origin = t, date = td$date[t + h],
            variable = v, measure = "ye", h = h,
            point = mean(xye), real = ye_real,
            logdens = -safe_logs(ye_real, xye),
            crps = safe_crps(ye_real, xye),
            pit = mean(xye <= ye_real))
        }
        # cumulative level from the origin: sum of quarterly draws t+1..t+h
        if (vt$transform == "dlog" && h >= 2) {
          xcum <- rowSums(dr[, seq_len(h), v, drop = FALSE])
          cum_real <- sum(real_q[(t + 1):(t + h)])
          rows[[length(rows) + 1]] <- data.frame(
            member = member_name, origin = t, date = td$date[t + h],
            variable = v, measure = "cum", h = h,
            point = mean(xcum), real = cum_real,
            logdens = -safe_logs(cum_real, xcum),
            crps = safe_crps(cum_real, xcum),
            pit = mean(xcum <= cum_real))
        }
      }
    }
  }
  do.call(rbind, rows)
}

#' Stored draws in a flat structure for the combination layer:
#' list keyed "origin|variable|h" -> named list(member -> draw vector).
collect_draws <- function(oos_all, cfg) {
  out <- new.env(parent = emptyenv())
  for (m in names(oos_all)) {
    for (res in oos_all[[m]]) {
      t <- res$origin
      dr <- res$draws
      for (v in dimnames(dr)[[3]]) for (h in seq_len(dim(dr)[2])) {
        key <- paste(t, v, h, sep = "|")
        cur <- if (exists(key, out)) get(key, out) else list()
        cur[[m]] <- dr[, h, v]
        assign(key, cur, out)
      }
    }
  }
  out
}

#' Summary score table: mean logdens / CRPS / RMSE by member, variable,
#' measure, horizon. Aggregated separately so an NA in one score (e.g. a
#' combo row without stored draws for CRPS) does not drop the row from the
#' others. exclude_dates: optional realization dates to exclude (e.g. COVID
#' quarters, whose extreme realizations dominate mean log scores).
summarise_scores <- function(scores, exclude_dates = NULL) {
  if (!is.null(exclude_dates))
    # match at QUARTER granularity: realizations are stamped at quarter start
    # (real) or mid-quarter (synthetic), so exact-date matching silently
    # excludes nothing (the same .qidx gotcha covid.R guards against).
    scores <- scores[!(.qidx(scores$date) %in% .qidx(exclude_dates)), ]
  agg_one <- function(v) {
    out <- aggregate(scores[[v]],
                     by = scores[, c("member", "variable", "measure", "h")],
                     FUN = function(x) mean(x, na.rm = TRUE))
    names(out)[5] <- v
    out
  }
  ld <- agg_one("logdens"); cr <- agg_one("crps")
  rmse <- aggregate((scores$point - scores$real)^2,
                    by = scores[, c("member", "variable", "measure", "h")],
                    FUN = function(x) sqrt(mean(x, na.rm = TRUE)))
  names(rmse)[5] <- "rmse"
  Reduce(merge, list(ld, cr, rmse))
}

# ---- Diebold-Mariano ---------------------------------------------------------------

#' DM test with Newey-West variance (lag h-1) and the Harvey small-sample
#' correction. loss1/loss2 aligned by origin; H0: equal predictive accuracy.
dm_test <- function(loss1, loss2, h = 1) {
  d <- loss1 - loss2
  d <- d[is.finite(d)]
  n <- length(d)
  if (n < 8 || sd(d) < 1e-12) return(c(stat = NA_real_, p = NA_real_))
  dbar <- mean(d)
  L <- max(0, h - 1)
  g0 <- mean((d - dbar)^2)
  v <- g0
  if (L > 0) for (l in seq_len(min(L, n - 1))) {
    gl <- mean((d[(l + 1):n] - dbar) * (d[1:(n - l)] - dbar))
    v <- v + 2 * (1 - l / (L + 1)) * gl
  }
  v <- max(v, 1e-12)
  stat <- dbar / sqrt(v / n)
  k <- sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)   # Harvey et al. correction
  stat <- stat * k
  p <- 2 * pt(-abs(stat), df = n - 1)
  c(stat = stat, p = p)
}

#' DM comparisons of every member against a reference member, per variable,
#' measure and horizon, for both squared-error and CRPS losses.
dm_vs_reference <- function(scores, reference) {
  combos <- unique(scores[, c("variable", "measure", "h")])
  members <- setdiff(unique(scores$member), reference)
  # the Newey-West correction below indexes autocovariances positionally, which
  # assumes contiguous (consecutive) forecast origins; warn if that fails.
  ro <- sort(unique(scores$origin[scores$member == reference]))
  if (length(ro) > 1 && any(diff(ro) != 1))
    log_warn("dm_vs_reference: non-contiguous origins -> Newey-West lags are approximate")
  rows <- list()
  for (i in seq_len(nrow(combos))) {
    cb <- combos[i, ]
    base <- scores[scores$member == reference & scores$variable == cb$variable &
                   scores$measure == cb$measure & scores$h == cb$h, ]
    base <- base[order(base$origin), ]
    for (m in members) {
      alt <- scores[scores$member == m & scores$variable == cb$variable &
                    scores$measure == cb$measure & scores$h == cb$h, ]
      alt <- alt[order(alt$origin), ]
      common <- intersect(base$origin, alt$origin)
      if (length(common) < 8) next
      b <- base[match(common, base$origin), ]; a <- alt[match(common, alt$origin), ]
      # year-ended losses come from overlapping 4-quarter sums: serial
      # correlation extends to order h+3, not h-1
      h_nw <- if (cb$measure == "ye") cb$h + 3 else cb$h
      d_se   <- dm_test((a$point - a$real)^2, (b$point - b$real)^2, h = h_nw)
      d_crps <- dm_test(a$crps, b$crps, h = h_nw)
      rows[[length(rows) + 1]] <- data.frame(
        member = m, reference = reference, variable = cb$variable,
        measure = cb$measure, h = cb$h,
        dm_se = d_se["stat"], p_se = d_se["p"],
        dm_crps = d_crps["stat"], p_crps = d_crps["p"])
    }
  }
  do.call(rbind, rows)
}

# ---- section 9 self-checks ----------------------------------------------------------

#' No-look-ahead test: corrupting all data AFTER the origin must not change
#' the forecast. The corrupted FULL panel goes through harness_forecast --
#' the same entry point the evaluation uses -- so the test exercises the
#' actual slice point rather than pre-sliced data (which could never fail).
#' One representative member per distinct engine (covers each engine's RNG /
#' estimation path, which a single-member test would miss).
.one_per_engine <- function(cfg) {
  ms <- all_members(cfg)
  ms[!duplicated(vapply(ms, function(m) m$engine, ""))]
}

test_no_lookahead <- function(td, spec, cfg, members = NULL) {
  if (is.null(members)) members <- .one_per_engine(cfg)
  else if (!is.null(members$name)) members <- list(members)  # a single member
  origins <- oos_origins(td, cfg)
  # test the first origin AND, if available, an origin whose forecast horizon
  # spans a COVID quarter (exercises the no-leak of the COVID-scale path).
  covq <- if (!is.null(cfg$covid$quarters)) .qidx(unlist(cfg$covid$quarters)) else integer(0)
  spans <- function(t) length(covq) > 0 &&
    any(.qidx(td$date[(t + 1):min(t + cfg$horizons, nrow(td))]) %in% covq)
  cov_origin <- Filter(function(t) t < nrow(td) && spans(t), origins)
  test_origins <- unique(c(origins[1], if (length(cov_origin)) cov_origin[1]))
  ok <- TRUE
  for (member in members) for (t in test_origins) {
    if (t >= nrow(td)) next
    td_bad <- td
    td_bad[(t + 1):nrow(td), -1] <- td_bad[(t + 1):nrow(td), -1] * 1e6 + 999
    f1 <- harness_forecast(member, td, t, spec, cfg)
    f2 <- harness_forecast(member, td_bad, t, spec, cfg)
    if (!identical(f1$draws, f2$draws)) {
      ok <- FALSE
      log_error("NO-LOOK-AHEAD TEST FAILED (member {member$name}, origin {t})")
    }
  }
  if (ok) log_info("no-look-ahead test passed ({length(members)} engines x {length(test_origins)} origins)")
  ok
}

#' Reproducibility: same seed twice -> identical draws, for every engine.
test_reproducibility <- function(td, spec, cfg) {
  t <- oos_origins(td, cfg)[1]
  ok <- TRUE
  for (member in .one_per_engine(cfg)) {
    set.seed(derive_seed(cfg$master_seed, "repro"))
    f1 <- forecast_at_origin(member, td[seq_len(t), , drop = FALSE], spec, cfg)
    set.seed(derive_seed(cfg$master_seed, "repro"))
    f2 <- forecast_at_origin(member, td[seq_len(t), , drop = FALSE], spec, cfg)
    if (!identical(f1$draws, f2$draws)) {
      ok <- FALSE
      log_error("REPRODUCIBILITY TEST FAILED (member {member$name})")
    }
  }
  if (ok) log_info("reproducibility test passed ({length(.one_per_engine(cfg))} engines)")
  ok
}

#' Calibration check: PIT approximately uniform (chi-square on quintile bins).
check_calibration <- function(scores, alpha = 0.01) {
  out <- list()
  for (m in unique(scores$member)) {
    p <- scores$pit[scores$member == m & scores$measure == "q" & scores$h == 1]
    if (length(p) < 20) next
    cnt <- table(cut(p, seq(0, 1, 0.2), include.lowest = TRUE))
    chi <- suppressWarnings(chisq.test(cnt))
    out[[m]] <- chi$p.value
  }
  out
}

#' PIT moment tests at every horizon (Knuppel 2015 style). Multi-step PITs are
#' serially correlated by construction (overlapping forecast windows), so the
#' h=1 chi-square test above is wrong for h>1; instead test the first two
#' moments of the PIT with a Newey-West (lag h-1) variance:
#'   location:   E[pit] = 1/2   (bias: forecasts systematically too high/low)
#'   dispersion: E[(pit-1/2)^2] = 1/12  (intervals too wide / too narrow)
#' Returns a long data.frame (member, variable, h, n, mean_pit, p_location,
#' var_pit, p_dispersion); written to output/tables/pit_tests.csv.
pit_moment_tests <- function(scores, min_n = 15) {
  nw_p <- function(z, L) {
    z <- z[is.finite(z)]
    n <- length(z)
    if (n < min_n || sd(z) < 1e-12) return(NA_real_)
    zb <- mean(z)
    v <- mean((z - zb)^2)
    if (L > 0) for (l in seq_len(min(L, n - 1))) {
      gl <- mean((z[(l + 1):n] - zb) * (z[1:(n - l)] - zb))
      v <- v + 2 * (1 - l / (L + 1)) * gl
    }
    v <- max(v, 1e-12)
    2 * pt(-abs(zb / sqrt(v / n)), df = n - 1)
  }
  d <- scores[scores$measure == "q" & is.finite(scores$pit), ]
  combos <- unique(d[, c("member", "variable", "h")])
  rows <- lapply(seq_len(nrow(combos)), function(i) {
    cb <- combos[i, ]
    s <- d[d$member == cb$member & d$variable == cb$variable & d$h == cb$h, ]
    s <- s[order(s$origin), ]
    L <- cb$h - 1
    data.frame(member = cb$member, variable = cb$variable, h = cb$h,
               n = nrow(s), mean_pit = mean(s$pit),
               p_location = nw_p(s$pit - 0.5, L),
               var_pit = mean((s$pit - 0.5)^2),
               p_dispersion = nw_p((s$pit - 0.5)^2 - 1 / 12, L))
  })
  do.call(rbind, rows)
}
