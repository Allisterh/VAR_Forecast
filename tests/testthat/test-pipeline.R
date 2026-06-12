source_project <- function() {
  root <- testthat::test_path("..", "..")
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
    sys.source(f, envir = globalenv())
}
source_project()

# tiny config for fast structural tests
tiny_cfg <- function() {
  cfg <- load_config(testthat::test_path("..", "..", "config", "config.yml"))
  cfg$mcmc$ndraw <- 80; cfg$mcmc$nburn <- 20
  cfg$mcmc$forecast_draws <- 80; cfg$mcmc$store_draws <- 40
  cfg$mcmc$bench_ndraw <- 80; cfg$mcmc$bench_nburn <- 20
  cfg$evaluation$max_origins <- 4
  cfg$glp$enabled <- FALSE
  cfg
}

tiny_data <- function(cfg) {
  spec <- build_transform_spec(cfg)
  raw <- generate_synthetic_data(cfg, spec)
  list(spec = spec, td = transform_data(raw, spec))
}

test_that("synthetic DGP is block-exogenous by construction", {
  cfg <- tiny_cfg()
  d <- tiny_data(cfg)
  dgp <- attr(d$td, "dgp")
  dom <- which(d$spec$block == "domestic")
  expect_true(all(dgp$A1[seq_len(dgp$nf), dom] == 0))
  expect_true(all(dgp$A2[seq_len(dgp$nf), dom] == 0))
  # and domestic shocks never move foreign contemporaneously
  expect_true(all(dgp$G[seq_len(dgp$nf), dom] == 0))
})

test_that("gibbs engine recovers block exogeneity on synthetic ground truth", {
  cfg <- tiny_cfg(); cfg$mcmc$ndraw <- 200; cfg$mcmc$nburn <- 50
  d <- tiny_data(cfg)
  spec_s <- vars_for_set(d$spec, "small")
  y <- as.matrix(d$td[, spec_s$variable])
  member <- list(name = "t", kind = "var", engine = "gibbs", set = "small",
                 lags = 2, prior = list(lambda = 0.2, soc = FALSE, dio = FALSE))
  set.seed(1)
  post <- fit_var_member(y, member, spec_s, cfg)
  expect_lt(post$diagnostics$block_exog_max, 1e-2)
  expect_true(post$diagnostics$converged)
})

test_that("no-look-ahead: corrupted future data does not change forecasts", {
  cfg <- tiny_cfg()
  d <- tiny_data(cfg)
  expect_true(test_no_lookahead(d$td, d$spec, cfg))
})

test_that("reproducibility: same seed gives identical draws", {
  cfg <- tiny_cfg()
  d <- tiny_data(cfg)
  expect_true(test_reproducibility(d$td, d$spec, cfg))
})

test_that("combination weights are recursive (only past realizations used)", {
  # synthetic score table: member B only becomes good AFTER origin 20; weights
  # computed at t=20 must not see that
  sc <- expand.grid(member = c("A", "B"), origin = 1:30, h = 1:2,
                    stringsAsFactors = FALSE)
  sc$variable <- "v"; sc$measure <- "q"
  sc$logdens <- ifelse(sc$member == "A", -1,
                       ifelse(sc$origin + sc$h <= 20, -5, 10))
  cfg <- list(combination = list(min_train_origins = 3, forgetting = 1,
                                 shrink_kappa = 0))
  w <- combo_weights("logscore", sc, "v", 1:2, t = 20, c("A", "B"), cfg)
  expect_gt(w["A"], 0.99)   # B's future scores must be invisible at t=20
})

test_that("full mini pipeline: OOS -> scores -> combination runs and scores sanely", {
  cfg <- tiny_cfg()
  cfg$suite <- cfg$suite[1]            # one VAR
  cfg$benchmarks <- c("rw", "ucmean")
  cfg$combination$min_train_origins <- 2
  d <- tiny_data(cfg)
  future::plan(future::sequential)
  mem <- all_members(cfg)
  oos <- lapply(mem, function(m) run_oos_member(m, d$td, d$spec, cfg,
                                                cache_root = tempfile()))
  names(oos) <- vapply(mem, `[[`, "", "name")
  scores <- do.call(rbind, lapply(names(oos), function(m)
    score_member(m, oos[[m]], d$td, d$spec, cfg)))
  expect_true(all(c("logdens", "crps", "pit") %in% names(scores)))
  expect_true(all(is.finite(scores$crps)))
  expect_true(all(scores$pit >= 0 & scores$pit <= 1))
  draws_env <- collect_draws(oos, cfg)
  cmb <- combine_all(scores, draws_env, d$td, d$spec, cfg)
  expect_true(nrow(cmb$scores) > 0)
  # pooled log density must be a proper mixture: between min and max member ld
  one <- cmb$scores[1, ]
  mems <- scores[scores$origin == one$origin & scores$variable == one$variable &
                 scores$h == one$h & scores$measure == one$measure, ]
  expect_gte(one$logdens, min(mems$logdens) - 1e-6)
  expect_lte(one$logdens, max(mems$logdens) + log(nrow(mems)) + 1e-6)
})

test_that("transform spec ordering puts the foreign block first", {
  cfg <- tiny_cfg()
  spec <- build_transform_spec(cfg)
  fb <- which(spec$block == "foreign")
  expect_equal(fb, seq_along(fb))
})
