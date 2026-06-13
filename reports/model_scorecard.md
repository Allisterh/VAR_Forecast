# SOE-BVAR Suite — Model Scorecard

**Data source:** real (RBA + ABS + FRED)  
**Panel:** 1997Q4-2026Q1 | **Targets:** gdp_growth, cpi_inflation, unemp_rate, cash_rate | **Horizons:** 1-12 quarters  
**Evaluation:** expanding-window pseudo-real-time, 36 forecast origins; densities scored by CRPS and log predictive density, points by RMSE, all by horizon.  
**Diagnostics (§9):** all green (block exogeneity, MCMC convergence, forecast sanity, no-look-ahead, reproducibility).

## 1. The models

Every VAR member is **block-exogenous** (the domestic block never feeds back into the foreign block) and produces **iterated** density forecasts. Members are designed to fail differently; see DECISIONS.md for the full rationale.

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
| `combo_pool` | Optimal prediction pool (Hall-Mitchell / Geweke-Amisano) |
| `combo_bma` | Bayesian model averaging — reported as a diagnostic only |

## 2. Forecast performance

Lower CRPS / RMSE is better; higher log score is better. **Bold** = best in that column. Models ordered best-first (by mean over the shown horizons).

### 2a. Who forecasts best, by variable and horizon (CRPS)

Best single model (lowest mean CRPS in the bucket; value in parentheses):

| Variable | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| gdp_growth | small_sv (0.681) | small_tight (0.773) | small_tight (0.845) |
| cpi_inflation | rw (0.188) | ar4 (0.273) | ucmean (0.302) |
| unemp_rate | medium_minn (0.314) | small_ss (0.630) | small_ss (0.720) |
| cash_rate | small_sv (0.302) | combo_equal (0.887) | ar4 (1.442) |

### 2b. Density accuracy by variable and horizon (CRPS, lower better)

**Real GDP growth (qtr %) (`gdp_growth`)**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| combo_logscore | 0.127 | 0.207 | 0.237 | 0.129 | 0.528 |
| ucsv | 0.114 | 0.200 | 0.234 | 0.123 | 0.575 |
| small_tight | 0.125 | 0.195 | 0.268 | **0.109** | 0.572 |
| medium_conj | 0.118 | 0.185 | 0.273 | 0.116 | 0.578 |
| combo_pool | 0.114 | 0.184 | 0.236 | 0.137 | 0.612 |
| small_ss | **0.101** | 0.195 | 0.246 | 0.122 | 0.621 |
| combo_bma | 0.126 | 0.210 | 0.245 | 0.129 | 0.580 |
| ar4 | 0.123 | 0.147 | 0.199 | 0.110 | 0.719 |
| medium_minn | 0.143 | 0.227 | 0.248 | 0.165 | 0.516 |
| ucmean | 0.129 | **0.133** | **0.161** | 0.116 | 0.770 |
| small_minn | 0.124 | 0.207 | 0.262 | 0.118 | 0.611 |
| combo_equal | 0.122 | 0.201 | 0.244 | 0.153 | 0.602 |
| small_loose_p5 | 0.144 | 0.213 | 0.240 | 0.121 | 0.607 |
| small_sv | 0.131 | 0.245 | 0.271 | 0.148 | 0.545 |
| rw | 0.188 | 0.354 | 0.428 | 0.447 | **0.511** |

**Trimmed-mean CPI inflation (qtr %) (`cpi_inflation`)**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| medium_conj | 0.038 | 0.048 | 0.042 | 0.099 | 0.049 |
| ucsv | 0.048 | **0.045** | 0.049 | **0.093** | 0.053 |
| small_tight | **0.036** | 0.065 | 0.042 | 0.109 | **0.041** |
| small_sv | 0.045 | 0.058 | 0.043 | 0.093 | 0.062 |
| small_ss | 0.040 | 0.045 | 0.040 | 0.132 | 0.044 |
| small_minn | 0.037 | 0.051 | 0.046 | 0.123 | 0.050 |
| medium_minn | 0.047 | 0.058 | **0.039** | 0.108 | 0.063 |
| combo_logscore | 0.040 | 0.065 | 0.050 | 0.112 | 0.057 |
| combo_equal | 0.040 | 0.066 | 0.045 | 0.123 | 0.056 |
| combo_pool | 0.041 | 0.064 | 0.044 | 0.128 | 0.060 |
| combo_bma | 0.045 | 0.063 | 0.048 | 0.138 | 0.058 |
| small_loose_p5 | 0.058 | 0.050 | 0.048 | 0.138 | 0.063 |
| rw | 0.043 | 0.076 | 0.092 | 0.158 | 0.144 |
| ar4 | 0.062 | 0.150 | 0.099 | 0.257 | 0.085 |
| ucmean | 0.082 | 0.153 | 0.087 | 0.256 | 0.083 |

**Unemployment rate (%) (`unemp_rate`)**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| ar4 | **0.076** | **0.106** | **0.091** | 0.318 | **0.211** |
| small_minn | 0.114 | 0.137 | 0.136 | 0.334 | 0.287 |
| small_loose_p5 | 0.148 | 0.175 | 0.145 | **0.297** | 0.248 |
| small_tight | 0.130 | 0.155 | 0.137 | 0.365 | 0.266 |
| combo_logscore | 0.134 | 0.155 | 0.163 | 0.457 | 0.365 |
| small_ss | 0.142 | 0.174 | 0.164 | 0.456 | 0.339 |
| combo_pool | 0.114 | 0.165 | 0.168 | 0.479 | 0.352 |
| ucsv | 0.111 | 0.147 | 0.147 | 0.500 | 0.374 |
| combo_equal | 0.130 | 0.180 | 0.161 | 0.463 | 0.363 |
| combo_bma | 0.132 | 0.179 | 0.173 | 0.464 | 0.386 |
| rw | 0.141 | 0.181 | 0.162 | 0.514 | 0.382 |
| medium_conj | 0.146 | 0.207 | 0.205 | 0.496 | 0.376 |
| ucmean | 0.224 | 0.227 | 0.213 | 0.452 | 0.383 |
| medium_minn | 0.180 | 0.242 | 0.288 | 0.757 | 0.683 |
| small_sv | 0.185 | 0.261 | 0.313 | 0.778 | 0.681 |

**Cash rate (%) (`cash_rate`)**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| small_sv | 0.042 | 0.084 | 0.150 | 0.261 | 0.547 |
| ucsv | 0.084 | 0.122 | 0.191 | **0.249** | 0.487 |
| rw | 0.097 | 0.143 | 0.171 | 0.260 | 0.471 |
| small_tight | 0.074 | 0.119 | 0.165 | 0.287 | 0.518 |
| medium_conj | 0.074 | 0.123 | 0.216 | 0.313 | **0.444** |
| medium_minn | **0.030** | **0.060** | **0.134** | 0.291 | 0.734 |
| small_minn | 0.074 | 0.133 | 0.232 | 0.395 | 0.734 |
| combo_bma | 0.102 | 0.141 | 0.260 | 0.383 | 0.733 |
| combo_equal | 0.090 | 0.151 | 0.232 | 0.406 | 0.804 |
| combo_logscore | 0.087 | 0.157 | 0.253 | 0.364 | 0.838 |
| combo_pool | 0.097 | 0.148 | 0.237 | 0.389 | 0.864 |
| small_ss | 0.076 | 0.149 | 0.281 | 0.505 | 1.233 |
| ar4 | 0.100 | 0.173 | 0.300 | 0.497 | 1.241 |
| small_loose_p5 | 0.081 | 0.182 | 0.393 | 0.739 | 1.441 |
| ucmean | 2.024 | 2.145 | 2.138 | 2.120 | 2.983 |


### 2c. Point accuracy by variable and horizon (RMSE, lower better)

**`gdp_growth`**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| small_minn | — | — | — | — | — |
| small_ss | — | — | — | — | — |
| small_sv | — | — | — | — | — |
| small_loose_p5 | — | — | — | — | — |
| medium_minn | — | — | — | — | — |
| medium_conj | — | — | — | — | — |
| small_tight | — | — | — | — | — |
| rw | — | — | — | — | — |
| ar4 | — | — | — | — | — |
| ucsv | — | — | — | — | — |
| ucmean | — | — | — | — | — |
| combo_equal | — | — | — | — | — |
| combo_logscore | — | — | — | — | — |
| combo_pool | — | — | — | — | — |
| combo_bma | — | — | — | — | — |

**`cpi_inflation`**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| small_minn | — | — | — | — | — |
| small_ss | — | — | — | — | — |
| small_sv | — | — | — | — | — |
| small_loose_p5 | — | — | — | — | — |
| medium_minn | — | — | — | — | — |
| medium_conj | — | — | — | — | — |
| small_tight | — | — | — | — | — |
| rw | — | — | — | — | — |
| ar4 | — | — | — | — | — |
| ucsv | — | — | — | — | — |
| ucmean | — | — | — | — | — |
| combo_equal | — | — | — | — | — |
| combo_logscore | — | — | — | — | — |
| combo_pool | — | — | — | — | — |
| combo_bma | — | — | — | — | — |

**`unemp_rate`**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| small_minn | — | — | — | — | — |
| small_ss | — | — | — | — | — |
| small_sv | — | — | — | — | — |
| small_loose_p5 | — | — | — | — | — |
| medium_minn | — | — | — | — | — |
| medium_conj | — | — | — | — | — |
| small_tight | — | — | — | — | — |
| rw | — | — | — | — | — |
| ar4 | — | — | — | — | — |
| ucsv | — | — | — | — | — |
| ucmean | — | — | — | — | — |
| combo_equal | — | — | — | — | — |
| combo_logscore | — | — | — | — | — |
| combo_pool | — | — | — | — | — |
| combo_bma | — | — | — | — | — |

**`cash_rate`**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| small_minn | — | — | — | — | — |
| small_ss | — | — | — | — | — |
| small_sv | — | — | — | — | — |
| small_loose_p5 | — | — | — | — | — |
| medium_minn | — | — | — | — | — |
| medium_conj | — | — | — | — | — |
| small_tight | — | — | — | — | — |
| rw | — | — | — | — | — |
| ar4 | — | — | — | — | — |
| ucsv | — | — | — | — | — |
| ucmean | — | — | — | — | — |
| combo_equal | — | — | — | — | — |
| combo_logscore | — | — | — | — | — |
| combo_pool | — | — | — | — | — |
| combo_bma | — | — | — | — | — |


### 2d. Density calibration — mean log predictive density (higher better)

Averaged across the 4 targets. The mean log score is brutally sensitive to tail events: **−∞** means at least one origin where the realization fell outside that member's predictive support (an individual model can catastrophically miss a tail). The **combinations never do** — the linear pool always assigns positive density — which is the clearest single piece of evidence that pooling buys calibration and robustness.

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| combo_equal | -1.852 | **-1.504** | **-1.524** | **-1.757** | **-1.872** |
| combo_pool | **-1.837** | -1.548 | -1.670 | -1.954 | -2.028 |
| combo_logscore | -1.840 | -1.537 | -1.649 | -2.013 | -2.042 |
| combo_bma | -2.056 | -1.870 | -1.986 | -2.379 | -2.192 |
| medium_minn | -3.127 | -4.714 | -4.551 | -7.785 | -3.681 |
| ar4 | −∞ | −∞ | -7.791 | -10.428 | −∞ |
| medium_conj | −∞ | −∞ | −∞ | −∞ | -13.455 |
| rw | −∞ | -5.098 | -2.978 | -2.956 | -2.692 |
| small_loose_p5 | -6.527 | −∞ | -9.128 | -9.524 | -8.624 |
| small_minn | -6.573 | −∞ | −∞ | −∞ | −∞ |
| small_ss | −∞ | −∞ | −∞ | −∞ | -12.302 |
| small_sv | -4.600 | −∞ | -6.722 | -6.114 | -3.030 |
| small_tight | −∞ | −∞ | -10.973 | −∞ | -12.616 |
| ucmean | −∞ | −∞ | −∞ | −∞ | −∞ |
| ucsv | -5.404 | −∞ | -6.487 | -10.739 | -9.902 |


## 3. Do the combinations beat the best single model?

Mean CRPS over all 4 targets, by horizon bucket. The honest test of a pool is whether it beats both equal weights and the best individual member.

| Model | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| combo_equal | 0.386 | **0.647** | **0.858** |
| small_ss | 0.393 | 0.662 | 0.860 |
| ar4 | 0.411 | 0.666 | 0.868 |
| small_minn | 0.399 | 0.675 | 0.892 |
| combo_pool | 0.394 | 0.713 | 0.931 |
| small_tight | 0.418 | 0.704 | 0.934 |
| small_sv | 0.385 | 0.717 | 0.986 |
| medium_conj | 0.420 | 0.703 | 0.966 |
| ucsv | 0.435 | 0.722 | 0.944 |
| small_loose_p5 | 0.415 | 0.721 | 0.965 |
| combo_logscore | 0.425 | 0.761 | 0.942 |
| medium_minn | **0.383** | 0.714 | 1.038 |
| rw | 0.506 | 0.799 | 1.040 |
| combo_bma | 0.441 | 0.900 | 1.012 |
| ucmean | 0.830 | 0.895 | 0.950 |

## 4. Statistical significance (Diebold-Mariano)

How often each combination **significantly beats** the random-walk and AR(4) benchmarks on CRPS (Harvey-corrected, 10% level), counted over the 4 targets x 3 horizons {1, 4, 8} tested. A negative DM statistic means the combination is more accurate; significance is one-sided here.

| Combination | beats ar4 | beats rw |
|:--|:--|:--|
| combo_bma | 0 / 12 | 0 / 12 |
| combo_equal | 0 / 12 | 0 / 12 |
| combo_logscore | 0 / 12 | 0 / 12 |
| combo_pool | 0 / 12 | 0 / 12 |

## 5. How to read this

- **Point gains over the best member are modest by design** — equal weights are hard to beat (the forecast-combination puzzle). The pool's payoff is *calibration and robustness*: it insures against any single member failing, rather than always winning on accuracy.

- **CRPS and log score can disagree** on outlier-heavy windows (log score is far more sensitive to tail events); both are reported. A COVID-excluded variant is in `output/tables/scores_by_horizon_excovid.csv`.

- Year-ended (policy-relevant) transforms are scored too; this scorecard shows quarterly. See the full Quarto report for fan charts, PIT calibration, and weight-evolution plots.

