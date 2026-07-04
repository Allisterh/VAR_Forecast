# SOE-BVAR Suite — Model Scorecard

**Data source:** real (RBA + ABS + FRED)  
**Panel:** 1997Q4-2026Q1 | **Targets:** gdp_growth, cpi_inflation, unemp_rate, cash_rate | **Horizons:** 1-12 quarters  
**Evaluation:** expanding-window pseudo-real-time, 60 forecast origins; densities scored by CRPS and log predictive density, points by RMSE, all by horizon.  
**Two views, level first.** Performance is reported two ways. The headline is the **level** error (§2): the forecast error of each series' *level* at quarter t+h — for GDP and inflation (modelled as growth) the cumulative level from the forecast origin, where real GDP and the price level land, and for unemployment and the cash rate the rate level itself. Alongside it (§3) is the **quarterly-growth** view: each target's single-quarter outcome at t+h. The year-ended *growth* scores remain in `output/tables/scores_by_horizon.csv`.  
**Diagnostics (§9):** all green (block exogeneity, MCMC convergence, forecast sanity, no-look-ahead, reproducibility).

## 0. How to read this scorecard

**What this document is.** A league table for 12 forecasting models (plus their pooled combinations) that were all asked to do the same job: at each of 60 past quarters ("forecast origins"), using only data available at that date, forecast the Australian economy 1-12 quarters ahead. Every forecast was then scored against what actually happened. Models that look good here look good because they *predicted*, not because they *fitted*.

**The forecasts are densities, not points.** Each model produces a full probability distribution per variable per horizon. So we score three things: whether the centre was right (**RMSE** — root-mean-squared error of the point forecast), and whether the whole distribution was right (**CRPS** and **log score**). Intuition: CRPS is the average distance between the forecast distribution and the realized outcome — like an absolute error that also rewards honest uncertainty; *lower is better*. The log score is the log of the probability density the model assigned to what happened; *higher is better*, and it punishes a model brutally for calling an outcome "nearly impossible" that then occurs. CRPS and the log score usually agree; they diverge on outlier episodes (COVID), which is why both are shown.

**Horizons and buckets.** h = 1 means one quarter ahead, h = 12 three years ahead. Weights and summaries group horizons into buckets: **near** (h 1-4), **medium** (5-8), **far** (9-12).

**The two views.** §2 scores where the *level* of each series lands (the cumulative path — "where is the price level in two years?"), which is what policy usually cares about and is the headline. §3 scores each single quarter's outcome in isolation. A model can be good at one and poor at the other: the level view compounds any persistent bias, the quarterly view rewards getting the wiggles right.

**Members vs combinations.** Rows named `combo_*` are not models but *pools* — weighted mixtures of every model's density (§1c). The suite's premise is that no single model wins everywhere, so the production forecast is a pool; the members exist to make the pool good.

## 1. The models

Every VAR member is **block-exogenous** (the domestic block never feeds back into the foreign block) and produces **iterated** density forecasts. Members are designed to fail differently; see README.md for the full rationale.

### 1a. VAR members

| Model | Family | System | Lags | Shrinkage λ | Volatility | COVID |
|:--|:--|:--|:--|:--|:--|:--|
| `small_minn` | Independent Normal-inverse-Wishart (Gibbs) + SOC/DIO | small (8 var) | 4 | auto (GLP) | constant | LP scaling |
| `small_ss` | Steady-state (Villani) | small (8 var) | 4 | auto (GLP) | constant | LP scaling |
| `small_sv` | Stochastic volatility, equation-by-equation | small (8 var) | 4 | auto (GLP) | stochastic (SV) | t-errors (SV-t) |
| `small_loose_p5` | Independent Normal-inverse-Wishart (Gibbs) | small (8 var) | 5 | 0.4 | constant | LP scaling |
| `medium_minn` | Stochastic volatility, equation-by-equation | medium (13 var) | 2 | auto (GLP) | stochastic (SV) | t-errors (SV-t) |
| `medium_conj` | Block-recursive conjugate NIW + SOC/DIO | medium (13 var) | 4 | 0.1 | constant | LP scaling |
| `small_tight` | Block-recursive conjugate NIW + SOC/DIO | small (8 var) | 4 | 0.05 | constant | LP scaling |
| `small_unres` | Independent Normal-inverse-Wishart (Gibbs) + SOC/DIO | small (8 var) | 4 | auto (GLP) | constant | LP scaling |

### 1b. Benchmark members (the bar every VAR must clear)

| Model | Description | COVID |
|:--|:--|:--|
| `rw` | Random walk | LP scaling |
| `ar4` | Bayesian AR(4) | LP scaling |
| `ucsv` | Unobserved components + stochastic volatility | t-errors (robust) |
| `ucmean` | Unconditional mean | LP scaling |

### 1c. Combination schemes (density pools)

Weights estimated **per target variable and per horizon bucket** (near 1-4, medium 5-8, far 9-12), shrunk toward equal weights, strictly recursive (no look-ahead).

| Scheme | How weights are set |
|:--|:--|
| `combo_equal` | Equal weights (the benchmark pool — hard to beat) |
| `combo_logscore` | Recursive log-score weights, with forgetting |
| `combo_crps` | Inverse discounted-mean-CRPS weights (outlier-robust) |
| `combo_pool` | Optimal prediction pool (Hall-Mitchell / Geweke-Amisano) |
| `combo_bma` | Bayesian model averaging — reported as a diagnostic only |

Performance weights train on the **same loss the headline tables report** (the level view) and **exclude the COVID realization window** (2020Q1-2021Q2) from training — see README.md D20 for why.

## 2. Forecast performance — levels (the headline)

These tables score the **level** of each series at t+h: the cumulative real GDP and price level for the growth-modelled variables, the rate level for unemployment and the cash rate. Level errors accumulate the whole forecast path, so they grow with horizon and are not comparable across horizons — read down each column, not across. Lower CRPS / RMSE is better; higher log score is better. **Bold** = best in that column; models ordered best-first (mean over the shown horizons).

**Best single model by variable and horizon bucket (CRPS)** — lowest mean CRPS in the bucket; value in parentheses:

| Variable | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| gdp_growth | small_sv (0.718) | small_sv (1.169) | small_sv (1.482) |
| cpi_inflation | combo_logscore (0.326) | combo_equal (0.997) | combo_equal (1.801) |
| unemp_rate | small_ss (0.257) | ar4 (0.470) | ar4 (0.574) |
| cash_rate | medium_minn (0.268) | rw (0.746) | rw (1.096) |


**Density accuracy by variable and horizon (CRPS, lower better)**

**Real GDP growth (qtr %) (`gdp_growth`) — cumulative level from the origin**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_sv | **0.447** | **0.921** | **1.325** | 1.544 |
| combo_crps | 0.468 | 0.946 | 1.374 | 1.543 |
| combo_bma | 0.465 | 0.942 | 1.409 | 1.540 |
| small_tight | 0.492 | 0.969 | 1.400 | **1.497** |
| small_unres | 0.490 | 0.976 | 1.403 | 1.513 |
| combo_equal | 0.469 | 0.950 | 1.383 | 1.589 |
| small_ss | 0.483 | 0.976 | 1.429 | 1.508 |
| combo_logscore | 0.479 | 0.968 | 1.403 | 1.561 |
| combo_pool | 0.486 | 0.967 | 1.423 | 1.565 |
| medium_conj | 0.489 | 0.989 | 1.430 | 1.546 |
| medium_minn | 0.453 | 0.930 | 1.390 | 1.684 |
| small_minn | 0.488 | 0.997 | 1.450 | 1.524 |
| small_loose_p5 | 0.465 | 0.977 | 1.443 | 1.608 |
| ucsv | 0.465 | 0.958 | 1.416 | 1.708 |
| ar4 | 0.485 | 0.976 | 1.489 | 1.775 |
| ucmean | 0.491 | 1.010 | 1.607 | 1.997 |
| rw | 0.690 | 2.402 | 4.918 | 7.944 |

**Trimmed-mean CPI inflation (qtr %) (`cpi_inflation`) — cumulative level from the origin**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_equal | 0.133 | 0.556 | **1.287** | **2.112** |
| combo_crps | 0.134 | 0.555 | 1.312 | 2.157 |
| small_ss | 0.139 | 0.588 | 1.369 | 2.154 |
| small_minn | 0.137 | 0.570 | 1.344 | 2.225 |
| small_unres | 0.137 | 0.568 | 1.346 | 2.231 |
| combo_pool | 0.136 | 0.568 | 1.406 | 2.290 |
| small_sv | 0.141 | 0.603 | 1.390 | 2.271 |
| combo_logscore | **0.126** | **0.552** | 1.458 | 2.357 |
| rw | 0.131 | 0.561 | 1.403 | 2.420 |
| medium_conj | 0.150 | 0.613 | 1.412 | 2.372 |
| medium_minn | 0.137 | 0.582 | 1.410 | 2.442 |
| combo_bma | 0.147 | 0.629 | 1.562 | 2.281 |
| ucsv | 0.142 | 0.640 | 1.498 | 2.431 |
| small_tight | 0.155 | 0.651 | 1.473 | 2.432 |
| small_loose_p5 | 0.153 | 0.651 | 1.529 | 2.448 |
| ar4 | 0.142 | 0.679 | 1.554 | 2.465 |
| ucmean | 0.187 | 0.812 | 1.723 | 2.630 |

**Unemployment rate (%) (`unemp_rate`) — rate level at t+h**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| ar4 | 0.130 | **0.371** | **0.521** | 0.605 |
| small_ss | **0.121** | 0.376 | 0.544 | **0.596** |
| combo_equal | 0.126 | 0.374 | 0.550 | 0.609 |
| combo_crps | 0.127 | 0.373 | 0.552 | 0.612 |
| combo_logscore | 0.126 | 0.376 | 0.562 | 0.641 |
| combo_pool | 0.125 | 0.377 | 0.566 | 0.642 |
| combo_bma | 0.122 | 0.385 | 0.572 | 0.646 |
| small_unres | 0.125 | 0.402 | 0.581 | 0.622 |
| ucsv | 0.195 | 0.390 | 0.548 | 0.624 |
| rw | 0.133 | 0.396 | 0.574 | 0.667 |
| small_minn | 0.124 | 0.402 | 0.596 | 0.656 |
| small_tight | 0.128 | 0.409 | 0.597 | 0.644 |
| medium_minn | 0.123 | 0.385 | 0.610 | 0.693 |
| small_sv | 0.128 | 0.400 | 0.635 | 0.701 |
| medium_conj | 0.131 | 0.414 | 0.636 | 0.714 |
| small_loose_p5 | 0.128 | 0.441 | 0.662 | 0.721 |
| ucmean | 0.550 | 0.561 | 0.613 | 0.664 |

**Cash rate (%) (`cash_rate`) — rate level at t+h**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| rw | 0.160 | 0.495 | **0.887** | **1.208** |
| small_tight | 0.142 | 0.501 | 0.971 | 1.341 |
| ucsv | 0.153 | 0.528 | 0.972 | 1.307 |
| small_minn | 0.130 | 0.481 | 0.967 | 1.386 |
| combo_logscore | 0.106 | 0.474 | 1.016 | 1.419 |
| combo_pool | 0.104 | 0.475 | 1.029 | 1.418 |
| combo_bma | 0.105 | 0.487 | 1.021 | 1.417 |
| combo_crps | 0.129 | 0.502 | 1.021 | 1.444 |
| medium_conj | 0.157 | 0.499 | 0.984 | 1.459 |
| combo_equal | 0.142 | 0.523 | 1.025 | 1.430 |
| small_ss | 0.132 | 0.532 | 1.100 | 1.538 |
| small_unres | 0.134 | 0.549 | 1.118 | 1.585 |
| medium_minn | **0.091** | **0.473** | 1.158 | 1.791 |
| small_sv | 0.111 | 0.536 | 1.205 | 1.679 |
| ar4 | 0.160 | 0.674 | 1.287 | 1.724 |
| small_loose_p5 | 0.158 | 0.717 | 1.386 | 1.989 |
| ucmean | 1.554 | 1.751 | 1.924 | 2.058 |


**Point accuracy by variable and horizon (RMSE, lower better)**

**`gdp_growth`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_unres | 1.270 | 2.022 | **2.699** | **2.839** |
| small_ss | 1.262 | 2.026 | 2.774 | 2.858 |
| combo_bma | 1.255 | 2.019 | 2.728 | 2.923 |
| small_sv | 1.239 | **1.981** | 2.732 | 2.996 |
| small_tight | 1.269 | 2.051 | 2.760 | 2.868 |
| small_minn | 1.273 | 2.060 | 2.793 | 2.892 |
| small_loose_p5 | **1.239** | 2.028 | 2.780 | 3.081 |
| combo_logscore | 1.264 | 2.065 | 2.806 | 3.058 |
| medium_conj | 1.273 | 2.105 | 2.847 | 3.007 |
| combo_pool | 1.270 | 2.073 | 2.829 | 3.063 |
| ucsv | 1.275 | 2.050 | 2.818 | 3.131 |
| ar4 | 1.284 | 2.068 | 2.794 | 3.138 |
| combo_crps | 1.269 | 2.114 | 2.885 | 3.191 |
| medium_minn | 1.263 | 2.062 | 2.897 | 3.286 |
| ucmean | 1.274 | 2.080 | 2.873 | 3.314 |
| combo_equal | 1.276 | 2.181 | 3.080 | 3.603 |
| rw | 1.874 | 6.069 | 11.625 | 17.624 |

**`cpi_inflation`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| ar4 | 0.274 | 1.183 | **2.369** | **3.386** |
| ucmean | 0.361 | 1.338 | 2.481 | 3.446 |
| combo_equal | 0.273 | 1.182 | 2.532 | 3.845 |
| combo_crps | 0.273 | 1.188 | 2.588 | 3.948 |
| small_ss | 0.276 | 1.235 | 2.660 | 3.926 |
| small_unres | 0.273 | 1.235 | 2.663 | 3.972 |
| small_minn | 0.276 | 1.236 | 2.669 | 4.022 |
| combo_pool | 0.272 | 1.218 | 2.742 | 4.138 |
| combo_logscore | 0.247 | 1.168 | 2.844 | 4.244 |
| small_tight | 0.328 | 1.355 | 2.767 | 4.118 |
| rw | **0.228** | **1.084** | 2.737 | 4.540 |
| small_sv | 0.297 | 1.304 | 2.801 | 4.299 |
| combo_bma | 0.310 | 1.344 | 3.011 | 4.043 |
| medium_conj | 0.325 | 1.357 | 2.835 | 4.316 |
| ucsv | 0.297 | 1.321 | 2.876 | 4.386 |
| medium_minn | 0.283 | 1.249 | 2.823 | 4.541 |
| small_loose_p5 | 0.305 | 1.349 | 2.945 | 4.495 |

**`unemp_rate`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| ar4 | 0.301 | 0.716 | **0.943** | **1.007** |
| small_ss | **0.277** | **0.715** | 0.990 | 1.024 |
| combo_pool | 0.288 | 0.721 | 0.999 | 1.060 |
| combo_logscore | 0.290 | 0.726 | 1.006 | 1.061 |
| combo_bma | 0.281 | 0.749 | 1.029 | 1.058 |
| combo_crps | 0.288 | 0.730 | 1.021 | 1.083 |
| combo_equal | 0.292 | 0.725 | 1.025 | 1.087 |
| small_unres | 0.285 | 0.754 | 1.055 | 1.088 |
| small_tight | 0.294 | 0.765 | 1.093 | 1.134 |
| small_minn | 0.285 | 0.764 | 1.098 | 1.158 |
| ucsv | 0.483 | 0.783 | 1.015 | 1.071 |
| rw | 0.314 | 0.778 | 1.096 | 1.186 |
| small_sv | 0.291 | 0.760 | 1.132 | 1.210 |
| medium_minn | 0.280 | 0.771 | 1.156 | 1.291 |
| small_loose_p5 | 0.291 | 0.823 | 1.191 | 1.308 |
| medium_conj | 0.305 | 0.840 | 1.200 | 1.293 |
| ucmean | 1.005 | 1.032 | 1.080 | 1.131 |

**`cash_rate`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| rw | 0.297 | 0.998 | **1.678** | **2.115** |
| ucsv | 0.299 | 1.008 | 1.696 | 2.129 |
| combo_logscore | 0.227 | **0.890** | 1.739 | 2.277 |
| combo_pool | 0.231 | 0.900 | 1.741 | 2.277 |
| combo_bma | 0.228 | 0.915 | 1.791 | 2.282 |
| small_minn | 0.249 | 0.923 | 1.726 | 2.339 |
| small_ss | 0.248 | 0.944 | 1.762 | 2.321 |
| combo_crps | 0.269 | 0.941 | 1.769 | 2.405 |
| small_tight | 0.266 | 0.992 | 1.780 | 2.348 |
| combo_equal | 0.330 | 0.959 | 1.748 | 2.369 |
| small_unres | 0.252 | 0.995 | 1.836 | 2.434 |
| medium_conj | 0.385 | 1.028 | 1.856 | 2.628 |
| small_sv | 0.233 | 0.945 | 1.982 | 2.795 |
| ar4 | 0.285 | 1.170 | 2.025 | 2.515 |
| medium_minn | **0.216** | 0.894 | 2.022 | 3.136 |
| small_loose_p5 | 0.318 | 1.223 | 2.218 | 3.101 |
| ucmean | 2.448 | 2.613 | 2.747 | 2.896 |


**Density calibration — mean log predictive density (higher better)**

Averaged across the 4 targets. **−∞** marks at least one origin where the realization fell outside a member's predictive support; the combinations, which always assign positive density, never do.

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_crps | -1.067 | **-1.310** | -1.959 | -2.317 |
| combo_equal | -1.084 | -1.314 | **-1.953** | **-2.304** |
| combo_logscore | -1.013 | -1.329 | -2.044 | -2.423 |
| combo_pool | **-1.012** | -1.335 | -2.062 | -2.428 |
| combo_bma | -1.369 | -1.675 | -2.511 | -2.841 |
| small_sv | -2.727 | -2.604 | -3.457 | -2.781 |
| small_loose_p5 | -3.869 | -4.250 | -6.359 | -4.263 |
| ar4 | −∞ | -2.933 | -3.506 | -4.157 |
| medium_conj | −∞ | -3.164 | -6.267 | -6.395 |
| medium_minn | −∞ | -2.535 | -2.526 | -2.955 |
| rw | −∞ | -3.201 | -3.290 | -3.259 |
| small_minn | −∞ | -3.692 | -6.646 | -6.072 |
| small_ss | −∞ | -3.836 | -6.571 | -5.222 |
| small_tight | −∞ | -5.971 | -9.082 | -8.269 |
| small_unres | −∞ | -4.310 | -6.782 | -5.934 |
| ucmean | −∞ | -4.144 | -4.952 | -7.362 |
| ucsv | −∞ | -4.265 | -5.466 | -4.817 |


## 3. Forecast performance — quarterly growth

The same scoring on each target's **single-quarter outcome** at t+h: the quarterly growth rate for real GDP and CPI inflation, and — since they are modelled in levels — the rate level for unemployment and the cash rate (identical to §2 for those two). Quarterly growth is a noisier, more local target than the accumulated level: it does not compound the path, so errors do not grow mechanically with horizon and the rankings can differ from the level view. Same conventions as §2.

**Best single model by variable and horizon bucket (CRPS)** — lowest mean CRPS in the bucket; value in parentheses:

| Variable | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| gdp_growth | small_sv (0.483) | small_unres (0.514) | small_unres (0.533) |
| cpi_inflation | combo_crps (0.158) | combo_equal (0.207) | ucmean (0.229) |
| unemp_rate | small_ss (0.257) | ar4 (0.470) | ar4 (0.574) |
| cash_rate | medium_minn (0.268) | rw (0.746) | rw (1.096) |


**Density accuracy by variable and horizon (CRPS, lower better)**

**Real GDP growth (qtr %) (`gdp_growth`) — quarterly growth rate (%)**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_crps | 0.468 | 0.502 | 0.532 | 0.550 |
| small_unres | 0.490 | **0.501** | **0.519** | 0.543 |
| combo_equal | 0.469 | 0.505 | 0.535 | 0.550 |
| small_tight | 0.492 | 0.509 | 0.520 | 0.543 |
| small_ss | 0.483 | 0.510 | 0.526 | 0.547 |
| combo_bma | 0.465 | 0.512 | 0.535 | 0.555 |
| small_loose_p5 | 0.465 | 0.512 | 0.536 | 0.556 |
| medium_conj | 0.489 | 0.518 | 0.523 | **0.542** |
| small_sv | **0.447** | 0.504 | 0.554 | 0.567 |
| combo_pool | 0.486 | 0.512 | 0.526 | 0.548 |
| combo_logscore | 0.479 | 0.514 | 0.532 | 0.549 |
| small_minn | 0.488 | 0.512 | 0.528 | 0.546 |
| ar4 | 0.485 | 0.503 | 0.528 | 0.559 |
| medium_minn | 0.453 | 0.512 | 0.556 | 0.579 |
| ucmean | 0.491 | 0.510 | 0.541 | 0.565 |
| ucsv | 0.465 | 0.513 | 0.571 | 0.585 |
| rw | 0.690 | 0.849 | 0.984 | 1.030 |

**Trimmed-mean CPI inflation (qtr %) (`cpi_inflation`) — quarterly growth rate (%)**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_equal | 0.133 | **0.181** | 0.218 | 0.241 |
| combo_crps | 0.134 | 0.182 | 0.225 | 0.244 |
| small_ss | 0.139 | 0.191 | 0.231 | 0.237 |
| ar4 | 0.142 | 0.197 | 0.221 | 0.239 |
| small_minn | 0.137 | 0.188 | 0.233 | 0.248 |
| small_unres | 0.137 | 0.188 | 0.235 | 0.247 |
| combo_pool | 0.136 | 0.193 | 0.243 | 0.253 |
| small_sv | 0.141 | 0.195 | 0.237 | 0.257 |
| combo_logscore | **0.126** | 0.195 | 0.253 | 0.262 |
| ucmean | 0.187 | 0.198 | **0.216** | **0.235** |
| small_tight | 0.155 | 0.196 | 0.234 | 0.254 |
| medium_conj | 0.150 | 0.194 | 0.238 | 0.265 |
| medium_minn | 0.137 | 0.194 | 0.247 | 0.279 |
| combo_bma | 0.147 | 0.200 | 0.264 | 0.253 |
| ucsv | 0.142 | 0.202 | 0.252 | 0.268 |
| rw | 0.131 | 0.198 | 0.268 | 0.306 |
| small_loose_p5 | 0.153 | 0.209 | 0.267 | 0.279 |

**Unemployment rate (%) (`unemp_rate`) — rate level at t+h (as §2)**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| ar4 | 0.130 | **0.371** | **0.521** | 0.605 |
| small_ss | **0.121** | 0.376 | 0.544 | **0.596** |
| combo_equal | 0.126 | 0.374 | 0.550 | 0.609 |
| combo_crps | 0.127 | 0.373 | 0.552 | 0.612 |
| combo_logscore | 0.126 | 0.376 | 0.562 | 0.641 |
| combo_pool | 0.125 | 0.377 | 0.566 | 0.642 |
| combo_bma | 0.122 | 0.385 | 0.572 | 0.646 |
| small_unres | 0.125 | 0.402 | 0.581 | 0.622 |
| ucsv | 0.195 | 0.390 | 0.548 | 0.624 |
| rw | 0.133 | 0.396 | 0.574 | 0.667 |
| small_minn | 0.124 | 0.402 | 0.596 | 0.656 |
| small_tight | 0.128 | 0.409 | 0.597 | 0.644 |
| medium_minn | 0.123 | 0.385 | 0.610 | 0.693 |
| small_sv | 0.128 | 0.400 | 0.635 | 0.701 |
| medium_conj | 0.131 | 0.414 | 0.636 | 0.714 |
| small_loose_p5 | 0.128 | 0.441 | 0.662 | 0.721 |
| ucmean | 0.550 | 0.561 | 0.613 | 0.664 |

**Cash rate (%) (`cash_rate`) — rate level at t+h (as §2)**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| rw | 0.160 | 0.495 | **0.887** | **1.208** |
| small_tight | 0.142 | 0.501 | 0.971 | 1.341 |
| ucsv | 0.153 | 0.528 | 0.972 | 1.307 |
| small_minn | 0.130 | 0.481 | 0.967 | 1.386 |
| combo_logscore | 0.106 | 0.474 | 1.016 | 1.419 |
| combo_pool | 0.104 | 0.475 | 1.029 | 1.418 |
| combo_bma | 0.105 | 0.487 | 1.021 | 1.417 |
| combo_crps | 0.129 | 0.502 | 1.021 | 1.444 |
| medium_conj | 0.157 | 0.499 | 0.984 | 1.459 |
| combo_equal | 0.142 | 0.523 | 1.025 | 1.430 |
| small_ss | 0.132 | 0.532 | 1.100 | 1.538 |
| small_unres | 0.134 | 0.549 | 1.118 | 1.585 |
| medium_minn | **0.091** | **0.473** | 1.158 | 1.791 |
| small_sv | 0.111 | 0.536 | 1.205 | 1.679 |
| ar4 | 0.160 | 0.674 | 1.287 | 1.724 |
| small_loose_p5 | 0.158 | 0.717 | 1.386 | 1.989 |
| ucmean | 1.554 | 1.751 | 1.924 | 2.058 |


**Point accuracy by variable and horizon (RMSE, lower better)**

**`gdp_growth`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_sv | 1.239 | **1.283** | 1.355 | 1.388 |
| combo_bma | 1.255 | 1.296 | 1.340 | 1.387 |
| small_unres | 1.270 | 1.289 | 1.340 | 1.382 |
| small_tight | 1.269 | 1.299 | **1.336** | 1.380 |
| small_loose_p5 | **1.239** | 1.309 | 1.346 | 1.395 |
| small_ss | 1.262 | 1.295 | 1.346 | 1.388 |
| combo_logscore | 1.264 | 1.297 | 1.343 | 1.388 |
| combo_pool | 1.270 | 1.300 | 1.342 | 1.388 |
| combo_crps | 1.269 | 1.295 | 1.350 | 1.386 |
| small_minn | 1.273 | 1.298 | 1.343 | 1.388 |
| medium_conj | 1.273 | 1.312 | 1.343 | **1.378** |
| ucsv | 1.275 | 1.292 | 1.356 | 1.382 |
| ar4 | 1.284 | 1.291 | 1.338 | 1.403 |
| combo_equal | 1.276 | 1.296 | 1.358 | 1.390 |
| ucmean | 1.274 | 1.301 | 1.353 | 1.400 |
| medium_minn | 1.263 | 1.294 | 1.363 | 1.412 |
| rw | 1.874 | 1.811 | 1.961 | 1.941 |

**`cpi_inflation`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| ar4 | 0.274 | **0.366** | 0.388 | **0.403** |
| ucmean | 0.361 | 0.370 | **0.385** | 0.403 |
| combo_equal | 0.273 | 0.375 | 0.443 | 0.450 |
| small_ss | 0.276 | 0.393 | 0.448 | 0.435 |
| small_unres | 0.273 | 0.391 | 0.450 | 0.441 |
| combo_crps | 0.273 | 0.378 | 0.451 | 0.456 |
| small_minn | 0.276 | 0.392 | 0.451 | 0.453 |
| combo_pool | 0.272 | 0.394 | 0.470 | 0.463 |
| combo_logscore | 0.247 | 0.399 | 0.490 | 0.471 |
| small_tight | 0.328 | 0.401 | 0.447 | 0.456 |
| small_sv | 0.297 | 0.404 | 0.479 | 0.484 |
| rw | **0.228** | 0.390 | 0.512 | 0.538 |
| ucsv | 0.297 | 0.405 | 0.486 | 0.488 |
| combo_bma | 0.310 | 0.407 | 0.506 | 0.454 |
| medium_conj | 0.325 | 0.407 | 0.471 | 0.482 |
| medium_minn | 0.283 | 0.398 | 0.501 | 0.534 |
| small_loose_p5 | 0.305 | 0.425 | 0.516 | 0.503 |

**`unemp_rate`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| ar4 | 0.301 | 0.716 | **0.943** | **1.007** |
| small_ss | **0.277** | **0.715** | 0.990 | 1.024 |
| combo_pool | 0.288 | 0.721 | 0.999 | 1.060 |
| combo_logscore | 0.290 | 0.726 | 1.006 | 1.061 |
| combo_bma | 0.281 | 0.749 | 1.029 | 1.058 |
| combo_crps | 0.288 | 0.730 | 1.021 | 1.083 |
| combo_equal | 0.292 | 0.725 | 1.025 | 1.087 |
| small_unres | 0.285 | 0.754 | 1.055 | 1.088 |
| small_tight | 0.294 | 0.765 | 1.093 | 1.134 |
| small_minn | 0.285 | 0.764 | 1.098 | 1.158 |
| ucsv | 0.483 | 0.783 | 1.015 | 1.071 |
| rw | 0.314 | 0.778 | 1.096 | 1.186 |
| small_sv | 0.291 | 0.760 | 1.132 | 1.210 |
| medium_minn | 0.280 | 0.771 | 1.156 | 1.291 |
| small_loose_p5 | 0.291 | 0.823 | 1.191 | 1.308 |
| medium_conj | 0.305 | 0.840 | 1.200 | 1.293 |
| ucmean | 1.005 | 1.032 | 1.080 | 1.131 |

**`cash_rate`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| rw | 0.297 | 0.998 | **1.678** | **2.115** |
| ucsv | 0.299 | 1.008 | 1.696 | 2.129 |
| combo_logscore | 0.227 | **0.890** | 1.739 | 2.277 |
| combo_pool | 0.231 | 0.900 | 1.741 | 2.277 |
| combo_bma | 0.228 | 0.915 | 1.791 | 2.282 |
| small_minn | 0.249 | 0.923 | 1.726 | 2.339 |
| small_ss | 0.248 | 0.944 | 1.762 | 2.321 |
| combo_crps | 0.269 | 0.941 | 1.769 | 2.405 |
| small_tight | 0.266 | 0.992 | 1.780 | 2.348 |
| combo_equal | 0.330 | 0.959 | 1.748 | 2.369 |
| small_unres | 0.252 | 0.995 | 1.836 | 2.434 |
| medium_conj | 0.385 | 1.028 | 1.856 | 2.628 |
| small_sv | 0.233 | 0.945 | 1.982 | 2.795 |
| ar4 | 0.285 | 1.170 | 2.025 | 2.515 |
| medium_minn | **0.216** | 0.894 | 2.022 | 3.136 |
| small_loose_p5 | 0.318 | 1.223 | 2.218 | 3.101 |
| ucmean | 2.448 | 2.613 | 2.747 | 2.896 |


**Density calibration — mean log predictive density (higher better)**

Averaged across the 4 targets. **−∞** marks at least one origin where the realization fell outside a member's predictive support; the combinations, which always assign positive density, never do.

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_equal | -1.084 | **-1.089** | **-1.297** | **-1.422** |
| combo_crps | -1.067 | -1.091 | -1.315 | -1.456 |
| combo_logscore | -1.013 | -1.132 | -1.415 | -1.568 |
| combo_pool | **-1.012** | -1.139 | -1.420 | -1.563 |
| combo_bma | -1.369 | -1.473 | -1.719 | -1.939 |
| small_loose_p5 | -3.869 | -5.827 | -5.413 | -6.026 |
| ar4 | −∞ | -4.667 | -5.895 | −∞ |
| medium_conj | −∞ | −∞ | −∞ | -7.280 |
| medium_minn | −∞ | −∞ | -1.591 | -3.100 |
| rw | −∞ | -2.027 | -2.109 | -2.003 |
| small_minn | −∞ | −∞ | −∞ | −∞ |
| small_ss | −∞ | -6.033 | -7.017 | -5.720 |
| small_sv | -2.727 | −∞ | -4.410 | -2.313 |
| small_tight | −∞ | -6.520 | −∞ | -6.832 |
| small_unres | −∞ | −∞ | −∞ | −∞ |
| ucmean | −∞ | −∞ | −∞ | −∞ |
| ucsv | −∞ | −∞ | −∞ | -7.287 |


## 4. Do the combinations beat the best single model?

Mean CRPS over all 4 targets, by horizon bucket. The honest test of a pool is whether it beats both equal weights and the best individual member. (Buckets average errors of differing scale across horizons, so use this within a bucket to rank models, not to compare buckets.)


**4a. Levels**

| Model | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| combo_equal | 0.416 | **0.893** | **1.306** |
| combo_crps | 0.411 | 0.895 | 1.312 |
| small_minn | 0.423 | 0.916 | 1.327 |
| combo_pool | 0.410 | 0.927 | 1.347 |
| combo_logscore | 0.405 | 0.926 | 1.354 |
| small_ss | 0.424 | 0.933 | 1.337 |
| combo_bma | 0.416 | 0.955 | 1.347 |
| small_unres | 0.430 | 0.935 | 1.357 |
| small_tight | 0.439 | 0.937 | 1.350 |
| ucsv | 0.442 | 0.931 | 1.374 |
| medium_conj | 0.438 | 0.937 | 1.378 |
| small_sv | 0.417 | 0.942 | 1.408 |
| medium_minn | **0.400** | 0.932 | 1.467 |
| ar4 | 0.460 | 1.018 | 1.489 |
| small_loose_p5 | 0.467 | 1.052 | 1.539 |
| ucmean | 0.875 | 1.311 | 1.704 |
| rw | 0.619 | 1.553 | 2.638 |

**4b. Quarterly growth**

| Model | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| combo_crps | 0.303 | 0.517 | 0.668 |
| small_minn | 0.309 | **0.516** | 0.666 |
| combo_equal | 0.308 | 0.517 | 0.667 |
| small_tight | 0.318 | 0.518 | **0.658** |
| combo_pool | 0.300 | 0.526 | 0.672 |
| combo_logscore | 0.299 | 0.527 | 0.673 |
| ucsv | 0.328 | 0.525 | 0.660 |
| combo_bma | 0.303 | 0.535 | 0.678 |
| small_ss | 0.311 | 0.531 | 0.687 |
| medium_conj | 0.321 | 0.527 | 0.693 |
| small_unres | 0.316 | 0.543 | 0.704 |
| medium_minn | **0.296** | 0.549 | 0.765 |
| small_sv | 0.310 | 0.568 | 0.756 |
| ar4 | 0.334 | 0.569 | 0.733 |
| rw | 0.388 | 0.603 | 0.759 |
| small_loose_p5 | 0.349 | 0.630 | 0.826 |
| ucmean | 0.727 | 0.801 | 0.861 |

## 5. Statistical significance (Diebold-Mariano)

How often each combination **significantly beats** the random-walk and AR(4) benchmarks on CRPS (Harvey-corrected, 10% level), counted over the 4 targets x 3 horizons {4, 8, 12} tested. A negative DM statistic means the combination is more accurate; significance is one-sided here. The tests have low power on these high-variance losses (dominated by a few episodes such as 2020), so the combinations beat the benchmarks on *average* (§4) more often than they do *significantly*.


**5a. Levels**

| Combination | beats ar4 | beats rw |
|:--|:--|:--|
| combo_bma | 1 / 12 | 3 / 12 |
| combo_crps | 1 / 12 | 3 / 12 |
| combo_equal | 1 / 12 | 3 / 12 |
| combo_logscore | 1 / 12 | 3 / 12 |
| combo_pool | 1 / 12 | 3 / 12 |

Strongest results (significant at 5% or better):

- `combo_crps` beats `ar4` on **cash_rate** at h=4 (DM -3.72***)
- `combo_equal` beats `ar4` on **cash_rate** at h=4 (DM -3.57***)
- `combo_logscore` beats `ar4` on **cash_rate** at h=4 (DM -3.23***)
- `combo_pool` beats `ar4` on **cash_rate** at h=4 (DM -3.19***)
- `combo_bma` beats `ar4` on **cash_rate** at h=4 (DM -2.51**)
- `combo_bma` beats `rw` on **gdp_growth** at h=8 (DM -2.27**)

**5b. Quarterly growth**

| Combination | beats ar4 | beats rw |
|:--|:--|:--|
| combo_bma | 1 / 12 | 2 / 12 |
| combo_crps | 1 / 12 | 2 / 12 |
| combo_equal | 1 / 12 | 3 / 12 |
| combo_logscore | 1 / 12 | 2 / 12 |
| combo_pool | 1 / 12 | 2 / 12 |

Strongest results (significant at 5% or better):

- `combo_crps` beats `ar4` on **cash_rate** at h=4 (DM -3.72***)
- `combo_equal` beats `ar4` on **cash_rate** at h=4 (DM -3.57***)
- `combo_logscore` beats `ar4` on **cash_rate** at h=4 (DM -3.23***)
- `combo_pool` beats `ar4` on **cash_rate** at h=4 (DM -3.19***)
- `combo_bma` beats `ar4` on **cash_rate** at h=4 (DM -2.51**)
- `combo_equal` beats `rw` on **gdp_growth** at h=4 (DM -2.41**)

## 6. Model profiles

One entry per model: its specification, what makes it distinct, the role it plays in the suite, its strengths and failure modes, and where it actually ranks in this evaluation (the eval line is computed, not asserted). Full rationale is in README.md.

### 6a. VAR members

**`small_minn`**  
*Spec:* Independent Normal-inverse-Wishart (Gibbs); 8-variable small SOE core; 4 lags; GLP marginal-likelihood shrinkage; constant volatility; LP scaling (COVID); sum-of-coefficients + dummy-initial-observation priors.  
*Distinctive:* The workhorse Minnesota BVAR, estimated by Gibbs with an *independent* Normal-inverse-Wishart prior — the engine required to impose block exogeneity, since asymmetric (equation-specific) prior variances cannot be represented by a Kronecker/conjugate prior. Shrinkage is data-driven (GLP marginal-likelihood), and sum-of-coefficients + dummy-initial-observation priors discipline the I(1) levels (rates, real TWI).  
*Role:* The central, representative small-SOE BVAR — the reference point the other small members are deliberate variations around (tighter, looser, steady-state, SV).  
*Strengths & failure modes:* Solid all-rounder at short-to-medium horizons. Constant volatility means it leans on the LP scaling for 2020; without it the 2020 outliers would distort the coefficients.  
*In this evaluation — levels:* best individual model for cpi_inflation at medium (5-8).  
*In this evaluation — quarterly growth:* strongest at cpi_inflation (near (1-4)), ranked 2 of 12 individual models.  
*See:* README.md D3, D4, D5, D17

**`small_ss`**  
*Spec:* Steady-state (Villani); 8-variable small SOE core; 4 lags; GLP marginal-likelihood shrinkage; constant volatility; LP scaling (COVID).  
*Distinctive:* Reparametrised around its *unconditional means* (Villani steady state), with informative priors placed directly on long-run levels — inflation 2.5% (target midpoint), NAIRU 4.5%, neutral cash rate 3.5%, US potential growth — and on the foreign block's steady states, which the domestic forecast inherits.  
*Role:* The **long-horizon anchor**. An iterated VAR reverts to its unconditional mean at h = 8-12; this member makes that mean economically grounded rather than the raw sample average.  
*Strengths & failure modes:* Strongest at medium/far horizons for the mean-reverting targets. Vulnerable if a steady-state anchor is stale (e.g. a shifted neutral rate drags the long end); constant volatility under-disperses around 2020 absent the LP correction.  
*In this evaluation — levels:* best individual model for cpi_inflation at far (9-12); unemp_rate at near (1-4).  
*In this evaluation — quarterly growth:* best individual model for unemp_rate at near (1-4).  
*See:* README.md D3, D4, D5, D17

**`small_sv`**  
*Spec:* Stochastic volatility, equation-by-equation; 8-variable small SOE core; 4 lags; GLP marginal-likelihood shrinkage; stochastic volatility; t-errors (SV-t) (COVID).  
*Distinctive:* Stochastic volatility estimated equation-by-equation (the Carriero-Clark-Marcellino triangular factorisation), with t-distributed errors (the CCMM SV-t COVID treatment). Block exogeneity is *exact* — foreign equations simply drop the domestic regressors.  
*Role:* The **density-calibration specialist**: time-varying volatility tracks changing macro uncertainty and the fat tails absorb outliers instead of letting them widen the whole history.  
*Strengths & failure modes:* Its SV captures the time-varying conditional variance, so it stays well-calibrated when uncertainty shifts — competitive on the GDP level at medium horizons. The volatility state at the jump-off can over/under-shoot if the last few quarters were unusual; the small system limits cross-variable information.  
*In this evaluation — levels:* best individual model for gdp_growth at near (1-4); gdp_growth at medium (5-8); gdp_growth at far (9-12).  
*In this evaluation — quarterly growth:* best individual model for gdp_growth at near (1-4).  
*See:* README.md D3, D6, D17

**`small_loose_p5`**  
*Spec:* Independent Normal-inverse-Wishart (Gibbs); 8-variable small SOE core; 5 lags; fixed shrinkage λ=0.4; constant volatility; LP scaling (COVID).  
*Distinctive:* A deliberately *under-shrunk*, longer-lag variant — fixed λ = 0.4 (vs the ~0.1-0.15 the GLP procedure selects) and 5 lags — so the data speak more and richer dynamics can show through, at the cost of estimation noise.  
*Role:* The **loose / long-lag** diversity axis: it fails differently from the tightly-shrunk members and can capture dynamics they shrink away.  
*Strengths & failure modes:* Occasionally best at near-horizon density when the extra flexibility pays; noisier and prone to wider intervals at long horizons (the cost of light shrinkage in a short sample).  
*In this evaluation — levels:* strongest at gdp_growth (near (1-4)), ranked 4 of 12 individual models.  
*In this evaluation — quarterly growth:* strongest at gdp_growth (near (1-4)), ranked 3 of 12 individual models.  
*See:* README.md D4, D8

**`medium_minn`**  
*Spec:* Stochastic volatility, equation-by-equation; 13-variable medium system; 2 lags; GLP marginal-likelihood shrinkage; stochastic volatility; t-errors (SV-t) (COVID).  
*Distinctive:* The larger system — 13 variables (adds terms of trade, wages, employment, consumption, the 10y yield) with stochastic volatility and t-errors, but only 2 lags to keep the parameter count feasible; equation-by-equation estimation keeps the recursive loop tractable.  
*Role:* The **medium-system** axis (Banbura-Giannone-Reichlin): medium systems often forecast best given enough shrinkage, and the extra variables bring cross-sectional information the small core lacks.  
*Strengths & failure modes:* Strong short-horizon density (it tends to win the near bucket). The short lag length limits long-horizon dynamics, and more parameters mean more estimation uncertainty at the far end.  
*In this evaluation — levels:* best individual model for cash_rate at near (1-4).  
*In this evaluation — quarterly growth:* best individual model for cash_rate at near (1-4).  
*See:* README.md D1, D6, D17

**`medium_conj`**  
*Spec:* Block-recursive conjugate NIW; 13-variable medium system; 4 lags; fixed shrinkage λ=0.1; constant volatility; LP scaling (COVID); sum-of-coefficients + dummy-initial-observation priors.  
*Distinctive:* The medium system estimated by the fast *block-recursive conjugate* scheme (foreign VAR + domestic block conditioned on contemporaneous foreign values, the RBNZ Bloor-Matheson approach) — closed-form, so cheap even at 13 variables x 4 lags. Block exogeneity is exact by the recursive structure, not the prior.  
*Role:* The cheap medium workhorse; it complements `medium_minn` (conjugate constant-volatility vs SV) on the same large system.  
*Strengths & failure modes:* Tends to lead the far-horizon GDP level error. Constant volatility leans on LP for 2020; the conjugate Kronecker prior cannot represent asymmetric shrinkage, which is why block exogeneity comes from the recursive structure.  
*In this evaluation — levels:* strongest at cash_rate (medium (5-8)), ranked 4 of 12 individual models.  
*In this evaluation — quarterly growth:* strongest at gdp_growth (medium (5-8)), ranked 3 of 12 individual models.  
*See:* README.md D3, D5, D8

**`small_tight`**  
*Spec:* Block-recursive conjugate NIW; 8-variable small SOE core; 4 lags; fixed shrinkage λ=0.05; constant volatility; LP scaling (COVID); sum-of-coefficients + dummy-initial-observation priors.  
*Distinctive:* The heavily-shrunk small model — fixed λ = 0.05, far tighter than the GLP selection — pulling hard toward the persistence/random-walk prior, so it is parsimonious and low-variance.  
*Role:* The **tight** diversity axis and the long-horizon robustness member: heavy shrinkage buys stability where lightly-parametrised models wander.  
*Strengths & failure modes:* Best or near-best at far-horizon GDP (the tight prior stops it over-reacting). Can be too rigid at short horizons, missing genuine dynamics the looser members catch.  
*In this evaluation — levels:* strongest at gdp_growth (far (9-12)), ranked 2 of 12 individual models.  
*In this evaluation — quarterly growth:* strongest at gdp_growth (medium (5-8)), ranked 2 of 12 individual models.  
*See:* README.md D4, D5, D8

**`small_unres`**  
*Spec:* Independent Normal-inverse-Wishart (Gibbs); 8-variable small SOE core; 4 lags; GLP marginal-likelihood shrinkage; constant volatility; LP scaling (COVID); sum-of-coefficients + dummy-initial-observation priors.  
*Distinctive:* Identical to `small_minn` in every respect — same engine, lags, data-driven shrinkage, SOC/DIO priors — **except that block exogeneity is NOT imposed**: the prior lets Australian variables feed back into world activity, commodity prices and the fed funds rate.  
*Role:* The **control experiment** for the suite's central identifying restriction: the score gap between `small_unres` and `small_minn` measures, inside this evaluation, what imposing the small-open-economy structure buys (RBA RDP 2013-06 argues it matters; this member lets our own data confirm or deny).  
*Strengths & failure modes:* Economically misspecified by construction (it implies Australia moves world demand). Watch whether its forecasts deteriorate at medium/far horizons, where the spurious feedback loops compound; it is exempt from the block-exogeneity gate.  
*In this evaluation — levels:* strongest at cpi_inflation (near (1-4)), ranked 2 of 12 individual models.  
*In this evaluation — quarterly growth:* best individual model for gdp_growth at medium (5-8); gdp_growth at far (9-12).  
*See:* README.md D19

### 6b. Benchmark members

**`rw`**  
*Spec:* Random walk; LP scaling (COVID).  
*Distinctive:* The no-change forecast: the last observed value persists, with Gaussian increments scaled to the historical change.  
*Role:* The universal hard-to-beat short-horizon bar for persistent/level variables, and a pool member.  
*Strengths & failure modes:* Competitive at h = 1 for level variables; fails badly at long horizons — its level path runs away, which the level-error tables (§2) expose brutally.  
*In this evaluation — levels:* best individual model for cpi_inflation at near (1-4); cash_rate at medium (5-8); cash_rate at far (9-12).  
*In this evaluation — quarterly growth:* best individual model for cpi_inflation at near (1-4); cash_rate at medium (5-8); cash_rate at far (9-12).  
*See:* README.md D8

**`ar4`**  
*Spec:* Bayesian AR(4); LP scaling (COVID).  
*Distinctive:* A Bayesian AR(4) per variable with Minnesota-style lag shrinkage and a stationarity-truncated posterior.  
*Role:* The univariate-persistence bar — it isolates how much of the forecast is just own-history dynamics.  
*Strengths & failure modes:* Surprisingly strong for inflation at near/medium horizons, where univariate dynamics dominate; it cannot use cross-variable information, so it lags when that matters.  
*In this evaluation — levels:* best individual model for unemp_rate at medium (5-8); unemp_rate at far (9-12).  
*In this evaluation — quarterly growth:* best individual model for unemp_rate at medium (5-8); unemp_rate at far (9-12).  
*See:* README.md D8

**`ucsv`**  
*Spec:* Unobserved components + stochastic volatility; t-errors (robust) (COVID).  
*Distinctive:* Stock-Watson unobserved-components stochastic volatility per variable: a random-walk trend plus transitory noise, both with time-varying variances and outlier-robust t-errors.  
*Role:* The canonical inflation benchmark and a genuine density anchor for the other targets.  
*Strengths & failure modes:* Strong for inflation (its native use case); weaker for variables with richer multivariate dynamics. The trend/noise split is weakly identified, so it is gated on Monte-Carlo precision, not raw ESS.  
*In this evaluation — levels:* strongest at cash_rate (far (9-12)), ranked 2 of 12 individual models.  
*In this evaluation — quarterly growth:* strongest at cash_rate (far (9-12)), ranked 2 of 12 individual models.  
*See:* README.md D9, D17

**`ucmean`**  
*Spec:* Unconditional mean; LP scaling (COVID).  
*Distinctive:* A Gaussian density centred on the expanding-sample mean with the sample variance — the simplest possible density forecast.  
*Role:* The floor: the 'did the model beat just predicting the long-run average' bar.  
*Strengths & failure modes:* Unexpectedly competitive at long horizons for mean-reverting growth (everything reverts to the mean eventually); useless at short horizons where dynamics matter.  
*In this evaluation — levels:* strongest at unemp_rate (far (9-12)), ranked 8 of 12 individual models.  
*In this evaluation — quarterly growth:* best individual model for cpi_inflation at medium (5-8); cpi_inflation at far (9-12).  
*See:* README.md D8

### 6c. Combination schemes

**`combo_equal`**  
*Weights:* Equal weights on every member, per variable x horizon bucket.  
*In this evaluation:* The forecast-combination-puzzle benchmark and the recommended robust default. Hard to beat because it never over-fits weights.

**`combo_logscore`**  
*Weights:* Weights proportional to each member's recent log predictive score, with a forgetting factor, shrunk toward equal.  
*In this evaluation:* Adapts to which members are forecasting well lately; the shrinkage and forgetting guard against over-concentrating on a member that was lucky. Log scores are tail-sensitive, so this scheme reacts hardest to outlier episodes (weights train on the level loss with COVID realizations excluded; D20).

**`combo_crps`**  
*Weights:* Weights proportional to the inverse of each member's discounted mean CRPS, shrunk toward equal.  
*In this evaluation:* The outlier-robust performance weighting (D20): CRPS grows linearly rather than logarithmically in the miss distance, so one extreme quarter cannot dominate the weights the way it can for the log-score scheme.

**`combo_pool`**  
*Weights:* Optimal prediction pool (Hall-Mitchell / Geweke-Amisano): weights on the simplex that maximise the historical *pooled* log score, shrunk toward equal.  
*In this evaluation:* Unlike BMA it does not degenerate to a single model ('all models are false but useful').

**`combo_bma`**  
*Weights:* Bayesian model averaging by predictive likelihood — no shrinkage.  
*In this evaluation:* **Diagnostic only.** It concentrates weight on the single best-fitting member, so it answers 'which model does the data favour' rather than serving as a robust combination; reported, not recommended.


## 7. How to read this

- **Two complementary metrics.** The **level** error (§2) is the headline: for GDP and inflation it is the cumulative level from the forecast origin (where real GDP and the price level land h quarters out), for unemployment and the cash rate the rate level at t+h. It accumulates the whole forecast path, so it grows with horizon and exposes a model that gets the persistent/drift component wrong (e.g. the random walk). The **quarterly-growth** view (§3) scores each target's single-quarter outcome — a noisier, more local target that does not compound the path, so it does not grow mechanically with horizon and can rank models differently. For unemployment and the cash rate (modelled in levels) the two views coincide. The year-ended *growth* scores live in `output/tables/scores_by_horizon.csv`.

- **Point gains over the best member are modest by design** — equal weights are hard to beat (the forecast-combination puzzle). The pool's payoff is *calibration and robustness*: it insures against any single member failing, rather than always winning on accuracy.

- **CRPS and log score can disagree** on outlier-heavy windows (log score is far more sensitive to tail events); both are reported. A COVID-excluded variant is in `output/tables/scores_by_horizon_excovid.csv`. See the full Quarto report for fan charts, PIT calibration, and weight-evolution plots.

