# probabilities.R -- policy-relevant event probabilities from the final pool ---
#
# Computes P(event) by horizon from the final-origin predictive draws:
#   * P(year-ended GDP growth < threshold)          -- "year-ended contraction"
#   * P(year-ended trimmed-mean inflation in band)  -- the 2-3% target band
#   * P(unemployment rises by >= delta pp from the jump-off)
#
# Pooling is done at the PROBABILITY level: for a linear pool the event
# probability is exactly the weighted average of member event probabilities,
# and member probabilities are computed from each member's JOINT path draws
# (year-ended sums need within-path dependence across horizons, which the
# per-horizon mixture-resampled pooled draws do not preserve).
#
# NOTE this file is deliberately NOT in config_hash()'s estimation-file list
# (R/utils.R): it is pure reporting and must not invalidate the OOS cache.

#' Weight vector for (scheme, variable, h) from the final-forecast weight
#' table (ff$weights), resolving h to its horizon bucket.
.event_weights <- function(wtab, scheme, v, h, members, cfg) {
  buckets <- cfg$combination$horizon_buckets
  bn <- names(buckets)[vapply(buckets, function(hs) h %in% unlist(hs),
                              logical(1))][1]
  w <- wtab[wtab$scheme == scheme & wtab$variable == v & wtab$bucket == bn, ]
  out <- setNames(rep(1 / length(members), length(members)), members)
  if (nrow(w)) {
    m <- match(members, w$member)
    ok <- !is.na(m)
    out[ok] <- w$weight[m[ok]]
    out <- out / sum(out)
  }
  out
}

#' Year-ended draws for variable v at horizon h from one member's joint path
#' draws [ndraw, H, var]: 4-quarter sum ending at t+h, splicing realized
#' history for h < 4 (same convention as score_member / .ye_from_q).
.ye_draws <- function(dr, td, v, h) {
  t <- nrow(td)
  k <- (h - 3):h
  hist_part <- sum(td[[v]][t + k[k <= 0]])
  fc_idx <- k[k >= 1]
  rowSums(dr[, fc_idx, v, drop = FALSE]) + hist_part
}

#' Event-probability table from the final forecasts (ff = final_forecasts()).
#' Events and thresholds come from cfg$report$events; scheme selects whose
#' weights to pool with (default the equal-weight pool, the production
#' recommendation). Returns a long data.frame and writes it to out_dir.
event_probabilities <- function(ff, td, spec, cfg, scheme = "equal",
                                out_dir = "output/tables") {
  ev <- cfg$report$events
  if (is.null(ev)) return(NULL)
  H <- cfg$horizons
  members <- names(ff$member_draws)
  rows <- list()
  add <- function(v, event, h, prob) {
    rows[[length(rows) + 1]] <<- data.frame(
      scheme = scheme, variable = v, event = event, h = h, prob = prob)
  }
  pool_prob <- function(v, h, member_prob) {
    w <- .event_weights(ff$weights, scheme, v, h, members, cfg)
    sum(w * member_prob[members])
  }

  # 1. year-ended GDP contraction
  if (!is.null(ev$recession_ye_gdp_below) && "gdp_growth" %in% members_vars(ff)) {
    thr <- as.numeric(ev$recession_ye_gdp_below)
    for (h in seq_len(H)) {
      pm <- vapply(members, function(m)
        mean(.ye_draws(ff$member_draws[[m]], td, "gdp_growth", h) < thr),
        numeric(1))
      add("gdp_growth", sprintf("P(year-ended GDP growth < %g%%)", thr),
          h, pool_prob("gdp_growth", h, pm))
    }
  }
  # 2. year-ended trimmed-mean inflation in the target band
  if (!is.null(ev$inflation_band) && "cpi_inflation" %in% members_vars(ff)) {
    band <- as.numeric(unlist(ev$inflation_band))
    for (h in seq_len(H)) {
      pm <- vapply(members, function(m) {
        ye <- .ye_draws(ff$member_draws[[m]], td, "cpi_inflation", h)
        mean(ye >= band[1] & ye <= band[2])
      }, numeric(1))
      add("cpi_inflation",
          sprintf("P(year-ended trimmed-mean inflation in [%g, %g]%%)",
                  band[1], band[2]),
          h, pool_prob("cpi_inflation", h, pm))
    }
  }
  # 3. unemployment rises by >= delta pp from the jump-off level
  if (!is.null(ev$unemp_rise) && "unemp_rate" %in% members_vars(ff)) {
    dlt <- as.numeric(ev$unemp_rise)
    u0 <- td$unemp_rate[nrow(td)]
    for (h in seq_len(H)) {
      pm <- vapply(members, function(m)
        mean(ff$member_draws[[m]][, h, "unemp_rate"] - u0 >= dlt), numeric(1))
      add("unemp_rate", sprintf("P(unemployment rises >= %gpp)", dlt),
          h, pool_prob("unemp_rate", h, pm))
    }
  }
  tab <- do.call(rbind, rows)
  if (is.null(tab)) return(NULL)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(tab, file.path(out_dir, "event_probabilities.csv"),
            row.names = FALSE)
  tab
}

#' Variables available in the member draw arrays.
members_vars <- function(ff) dimnames(ff$member_draws[[1]])[[3]]

#' Figure: event probabilities by horizon.
plot_event_probs <- function(ep, out_dir = "output/figures") {
  if (is.null(ep) || !nrow(ep)) return(NULL)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  g <- ggplot2::ggplot(ep, ggplot2::aes(h, prob)) +
    ggplot2::geom_line(color = "steelblue4") +
    ggplot2::geom_point(size = 0.9, color = "steelblue4") +
    ggplot2::facet_wrap(~event, scales = "free_y") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = "Event probabilities by horizon (final pooled forecast)",
                  x = "horizon (quarters ahead)", y = "probability") +
    ggplot2::theme_minimal(base_size = 11)
  f <- file.path(out_dir, "event_probabilities.png")
  ggplot2::ggsave(f, g, width = 9, height = 4, dpi = 130)
  f
}
