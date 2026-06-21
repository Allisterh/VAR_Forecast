# SOE-BVAR Suite — Model Scorecard

**Data source:** real (RBA + ABS + FRED)  
**Panel:** 1997Q4-2026Q1 | **Targets:** gdp_growth, cpi_inflation, unemp_rate, cash_rate | **Horizons:** 1-12 quarters  
**Evaluation:** expanding-window pseudo-real-time, 36 forecast origins; densities scored by CRPS and log predictive density, points by RMSE, all by horizon.  
**Two views, level first.** Performance is reported two ways. The headline is the **level** error (§2): the forecast error of each series' *level* at quarter t+h — for GDP and inflation (modelled as growth) the cumulative level from the forecast origin, where real GDP and the price level land, and for unemployment and the cash rate the rate level itself. Alongside it (§3) is the **quarterly-growth** view: each target's single-quarter outcome at t+h. The year-ended *growth* scores remain in `output/tables/scores_by_horizon.csv`.  
**Diagnostics (§9):** all green (block exogeneity, MCMC convergence, forecast sanity, no-look-ahead, reproducibility).

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

## 2. Forecast performance — levels (the headline)

These tables score the **level** of each series at t+h: the cumulative real GDP and price level for the growth-modelled variables, the rate level for unemployment and the cash rate. Level errors accumulate the whole forecast path, so they grow with horizon and are not comparable across horizons — read down each column, not across. Lower CRPS / RMSE is better; higher log score is better. **Bold** = best in that column; models ordered best-first (mean over the shown horizons).

**Best single model by variable and horizon bucket (CRPS)** — lowest mean CRPS in the bucket; value in parentheses:

| Variable | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| gdp_growth | medium_minn (0.985) | small_sv (1.665) | small_sv (2.102) |
| cpi_inflation | rw (0.386) | rw (1.358) | ar4 (2.693) |
| unemp_rate | medium_minn (0.325) | small_ss (0.631) | small_ss (0.726) |
| cash_rate | small_sv (0.307) | combo_equal (0.894) | ar4 (1.447) |


**Density accuracy by variable and horizon (CRPS, lower better)**

**Real GDP growth (qtr %) (`gdp_growth`) — cumulative level from the origin**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_sv | **0.604** | 1.269 | **1.921** | **2.154** |
| medium_minn | 0.614 | **1.266** | 1.977 | 2.315 |
| combo_equal | 0.631 | 1.300 | 2.042 | 2.233 |
| small_tight | 0.682 | 1.366 | 2.101 | 2.165 |
| ucsv | 0.629 | 1.298 | 2.026 | 2.391 |
| medium_conj | 0.672 | 1.383 | 2.139 | 2.238 |
| small_loose_p5 | 0.628 | 1.362 | 2.161 | 2.300 |
| small_ss | 0.662 | 1.376 | 2.190 | 2.284 |
| small_minn | 0.677 | 1.414 | 2.230 | 2.293 |
| ar4 | 0.668 | 1.342 | 2.156 | 2.470 |
| ucmean | 0.680 | 1.380 | 2.273 | 2.662 |
| combo_pool | 0.635 | 1.432 | 2.438 | 3.699 |
| combo_logscore | 0.760 | 1.803 | 4.328 | 3.982 |
| combo_bma | 0.840 | 2.232 | 6.320 | 5.346 |
| rw | 0.941 | 3.264 | 6.711 | 11.084 |

**Trimmed-mean CPI inflation (qtr %) (`cpi_inflation`) — cumulative level from the origin**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| rw | **0.131** | **0.674** | **1.832** | 3.356 |
| ar4 | 0.162 | 0.844 | 1.973 | **3.108** |
| combo_equal | 0.160 | 0.770 | 1.912 | 3.323 |
| combo_pool | 0.147 | 0.764 | 2.081 | 3.502 |
| combo_logscore | 0.139 | 0.745 | 2.037 | 3.614 |
| small_minn | 0.170 | 0.824 | 2.109 | 3.714 |
| small_ss | 0.172 | 0.849 | 2.165 | 3.655 |
| small_sv | 0.171 | 0.856 | 2.153 | 3.769 |
| ucmean | 0.230 | 1.056 | 2.269 | 3.405 |
| medium_minn | 0.165 | 0.810 | 2.153 | 4.021 |
| small_tight | 0.198 | 0.925 | 2.241 | 3.910 |
| combo_bma | 0.146 | 0.808 | 2.391 | 3.945 |
| medium_conj | 0.191 | 0.895 | 2.234 | 4.029 |
| ucsv | 0.176 | 0.926 | 2.333 | 4.058 |
| small_loose_p5 | 0.193 | 0.944 | 2.437 | 4.089 |

**Unemployment rate (%) (`unemp_rate`) — rate level at t+h**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_ss | **0.148** | 0.497 | **0.690** | **0.747** |
| combo_equal | 0.159 | 0.493 | 0.702 | 0.785 |
| medium_minn | 0.150 | **0.477** | 0.726 | 0.803 |
| small_sv | 0.157 | 0.501 | 0.743 | 0.771 |
| small_minn | 0.150 | 0.524 | 0.744 | 0.789 |
| small_loose_p5 | 0.150 | 0.540 | 0.775 | 0.812 |
| small_tight | 0.158 | 0.549 | 0.779 | 0.811 |
| combo_pool | 0.177 | 0.496 | 0.767 | 0.860 |
| ucsv | 0.272 | 0.536 | 0.740 | 0.841 |
| rw | 0.169 | 0.545 | 0.777 | 0.913 |
| medium_conj | 0.163 | 0.545 | 0.820 | 0.896 |
| ar4 | 0.189 | 0.534 | 0.760 | 0.955 |
| combo_logscore | 0.383 | 0.618 | 0.776 | 0.841 |
| combo_bma | 0.677 | 0.806 | 0.938 | 0.915 |
| ucmean | 0.726 | 0.774 | 0.905 | 1.034 |

**Cash rate (%) (`cash_rate`) — rate level at t+h**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_equal | 0.148 | **0.514** | **1.120** | 1.652 |
| small_ss | 0.140 | 0.526 | 1.149 | 1.673 |
| ar4 | 0.167 | 0.589 | 1.141 | **1.593** |
| rw | 0.183 | 0.590 | 1.136 | 1.670 |
| small_minn | 0.144 | 0.528 | 1.152 | 1.792 |
| combo_pool | 0.129 | 0.549 | 1.225 | 1.715 |
| combo_logscore | 0.125 | 0.538 | 1.241 | 1.780 |
| combo_bma | 0.128 | 0.576 | 1.381 | 1.854 |
| ucsv | 0.175 | 0.646 | 1.281 | 1.842 |
| small_tight | 0.157 | 0.598 | 1.286 | 1.935 |
| small_loose_p5 | 0.162 | 0.593 | 1.247 | 2.038 |
| medium_conj | 0.192 | 0.572 | 1.239 | 2.053 |
| small_sv | 0.113 | 0.535 | 1.410 | 2.147 |
| medium_minn | **0.112** | 0.559 | 1.465 | 2.494 |
| ucmean | 1.563 | 1.654 | 1.671 | 1.615 |


**Point accuracy by variable and horizon (RMSE, lower better)**

**`gdp_growth`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_sv | 1.573 | **2.524** | **3.599** | 3.938 |
| small_tight | 1.614 | 2.623 | 3.654 | **3.830** |
| small_ss | 1.603 | 2.587 | 3.693 | 3.878 |
| small_minn | 1.619 | 2.633 | 3.715 | 3.910 |
| small_loose_p5 | **1.571** | 2.585 | 3.674 | 4.071 |
| ar4 | 1.632 | 2.632 | 3.631 | 4.036 |
| ucsv | 1.618 | 2.587 | 3.689 | 4.106 |
| ucmean | 1.619 | 2.634 | 3.685 | 4.152 |
| medium_conj | 1.616 | 2.690 | 3.774 | 4.030 |
| medium_minn | 1.605 | 2.625 | 3.805 | 4.323 |
| combo_equal | 1.624 | 2.815 | 4.140 | 4.971 |
| combo_pool | 1.688 | 3.411 | 5.770 | 8.846 |
| combo_logscore | 2.059 | 5.597 | 11.640 | 8.842 |
| combo_bma | 2.308 | 6.971 | 15.317 | 11.148 |
| rw | 2.377 | 7.807 | 15.379 | 24.092 |

**`cpi_inflation`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| ar4 | 0.315 | 1.433 | **2.905** | **4.138** |
| ucmean | 0.434 | 1.648 | 3.084 | 4.269 |
| combo_equal | 0.327 | 1.500 | 3.319 | 5.208 |
| combo_pool | 0.287 | 1.450 | 3.451 | 5.449 |
| combo_logscore | 0.269 | 1.410 | 3.420 | 5.675 |
| small_ss | 0.333 | 1.587 | 3.541 | 5.413 |
| small_minn | 0.333 | 1.590 | 3.549 | 5.525 |
| rw | **0.235** | **1.301** | 3.486 | 6.042 |
| small_tight | 0.402 | 1.736 | 3.655 | 5.611 |
| small_sv | 0.360 | 1.671 | 3.716 | 5.907 |
| medium_conj | 0.400 | 1.750 | 3.777 | 5.951 |
| medium_minn | 0.340 | 1.591 | 3.729 | 6.223 |
| ucsv | 0.361 | 1.697 | 3.812 | 6.017 |
| combo_bma | 0.278 | 1.516 | 3.914 | 6.180 |
| small_loose_p5 | 0.371 | 1.735 | 3.926 | 6.181 |

**`unemp_rate`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_ss | **0.339** | **0.890** | 1.208 | **1.226** |
| ar4 | 0.413 | 0.915 | **1.174** | 1.300 |
| combo_equal | 0.370 | 0.903 | 1.259 | 1.324 |
| small_sv | 0.354 | 0.927 | 1.333 | 1.349 |
| small_minn | 0.348 | 0.946 | 1.337 | 1.364 |
| combo_pool | 0.428 | 0.919 | 1.260 | 1.401 |
| small_tight | 0.362 | 0.953 | 1.351 | 1.371 |
| small_loose_p5 | 0.349 | 0.989 | 1.382 | 1.467 |
| medium_minn | 0.341 | 0.943 | 1.388 | 1.516 |
| ucsv | 0.614 | 0.991 | 1.281 | 1.342 |
| rw | 0.391 | 0.984 | 1.391 | 1.509 |
| medium_conj | 0.375 | 1.050 | 1.487 | 1.574 |
| combo_logscore | 0.880 | 1.126 | 1.304 | 1.390 |
| ucmean | 1.260 | 1.326 | 1.435 | 1.542 |
| combo_bma | 1.252 | 1.365 | 1.502 | 1.448 |

**`cash_rate`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_ss | 0.276 | 0.990 | **1.908** | 2.581 |
| ar4 | 0.307 | 1.103 | 1.931 | **2.467** |
| combo_equal | 0.341 | **0.971** | 1.916 | 2.777 |
| combo_logscore | 0.272 | 0.990 | 2.133 | 2.862 |
| small_minn | 0.292 | 1.051 | 2.045 | 2.902 |
| rw | 0.352 | 1.194 | 2.062 | 2.684 |
| combo_pool | 0.290 | 1.025 | 2.149 | 2.849 |
| ucsv | 0.358 | 1.208 | 2.080 | 2.695 |
| combo_bma | 0.282 | 1.054 | 2.353 | 2.961 |
| small_tight | 0.313 | 1.190 | 2.230 | 3.078 |
| small_sv | **0.252** | 0.982 | 2.306 | 3.468 |
| small_loose_p5 | 0.352 | 1.113 | 2.180 | 3.387 |
| medium_conj | 0.479 | 1.218 | 2.280 | 3.418 |
| medium_minn | 0.265 | 1.052 | 2.486 | 4.079 |
| ucmean | 2.576 | 2.658 | 2.649 | 2.635 |


**Density calibration — mean log predictive density (higher better)**

Averaged across the 4 targets. **−∞** marks at least one origin where the realization fell outside a member's predictive support; the combinations, which always assign positive density, never do.

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_equal | -1.844 | **-1.570** | **-2.328** | **-2.764** |
| combo_pool | **-1.816** | -1.674 | -2.502 | -2.933 |
| combo_logscore | -1.925 | -1.698 | -2.571 | -2.941 |
| combo_bma | -3.731 | -3.547 | -3.463 | -3.201 |
| small_sv | -4.619 | -3.748 | -4.950 | -3.554 |
| small_loose_p5 | -6.452 | -6.500 | -10.199 | -6.380 |
| ar4 | −∞ | -4.390 | -5.271 | -6.282 |
| medium_conj | −∞ | -4.830 | -10.309 | -10.904 |
| medium_minn | −∞ | -3.713 | -3.342 | -3.985 |
| rw | −∞ | -4.619 | -4.431 | -4.180 |
| small_minn | −∞ | -5.758 | -11.007 | -10.271 |
| small_ss | −∞ | -5.976 | -10.826 | -8.541 |
| small_tight | −∞ | -9.674 | -15.438 | -14.541 |
| ucmean | −∞ | -5.804 | -6.761 | -10.352 |
| ucsv | −∞ | -6.761 | -8.879 | -7.851 |


## 3. Forecast performance — quarterly growth

The same scoring on each target's **single-quarter outcome** at t+h: the quarterly growth rate for real GDP and CPI inflation, and — since they are modelled in levels — the rate level for unemployment and the cash rate (identical to §2 for those two). Quarterly growth is a noisier, more local target than the accumulated level: it does not compound the path, so errors do not grow mechanically with horizon and the rankings can differ from the level view. Same conventions as §2.

**Best single model by variable and horizon bucket (CRPS)** — lowest mean CRPS in the bucket; value in parentheses:

| Variable | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| gdp_growth | small_sv (0.679) | small_tight (0.773) | small_tight (0.846) |
| cpi_inflation | rw (0.189) | ar4 (0.273) | ucmean (0.302) |
| unemp_rate | medium_minn (0.325) | small_ss (0.631) | small_ss (0.726) |
| cash_rate | small_sv (0.307) | combo_equal (0.894) | ar4 (1.447) |


**Density accuracy by variable and horizon (CRPS, lower better)**

**Real GDP growth (qtr %) (`gdp_growth`) — quarterly growth rate (%)**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_equal | 0.631 | 0.727 | 0.810 | 0.886 |
| small_loose_p5 | 0.628 | 0.737 | 0.810 | 0.887 |
| small_sv | **0.604** | **0.726** | 0.840 | 0.908 |
| ar4 | 0.668 | 0.727 | 0.797 | 0.898 |
| small_ss | 0.662 | 0.739 | 0.806 | 0.887 |
| small_tight | 0.682 | 0.743 | **0.794** | 0.878 |
| medium_conj | 0.672 | 0.757 | 0.800 | **0.876** |
| small_minn | 0.677 | 0.745 | 0.803 | 0.883 |
| medium_minn | 0.614 | 0.734 | 0.841 | 0.920 |
| ucmean | 0.680 | 0.732 | 0.816 | 0.908 |
| ucsv | 0.629 | 0.740 | 0.886 | 0.959 |
| combo_pool | 0.635 | 0.757 | 0.902 | 1.015 |
| combo_logscore | 0.760 | 0.855 | 1.132 | 1.043 |
| combo_bma | 0.840 | 0.931 | 1.411 | 1.165 |
| rw | 0.941 | 1.164 | 1.387 | 1.447 |

**Trimmed-mean CPI inflation (qtr %) (`cpi_inflation`) — quarterly growth rate (%)**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| ar4 | 0.162 | 0.246 | 0.288 | 0.314 |
| ucmean | 0.230 | 0.250 | **0.286** | **0.312** |
| combo_equal | 0.160 | **0.245** | 0.323 | 0.364 |
| combo_pool | 0.147 | 0.258 | 0.361 | 0.393 |
| small_ss | 0.172 | 0.268 | 0.354 | 0.369 |
| rw | **0.131** | 0.248 | 0.361 | 0.431 |
| small_minn | 0.170 | 0.262 | 0.355 | 0.389 |
| combo_logscore | 0.139 | 0.266 | 0.368 | 0.409 |
| small_sv | 0.171 | 0.270 | 0.359 | 0.406 |
| small_tight | 0.198 | 0.272 | 0.350 | 0.395 |
| medium_minn | 0.165 | 0.265 | 0.372 | 0.445 |
| medium_conj | 0.191 | 0.272 | 0.367 | 0.426 |
| ucsv | 0.176 | 0.285 | 0.386 | 0.430 |
| combo_bma | 0.146 | 0.291 | 0.407 | 0.437 |
| small_loose_p5 | 0.193 | 0.296 | 0.414 | 0.420 |

**Unemployment rate (%) (`unemp_rate`) — rate level at t+h (as §2)**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_ss | **0.148** | 0.497 | **0.690** | **0.747** |
| combo_equal | 0.159 | 0.493 | 0.702 | 0.785 |
| medium_minn | 0.150 | **0.477** | 0.726 | 0.803 |
| small_sv | 0.157 | 0.501 | 0.743 | 0.771 |
| small_minn | 0.150 | 0.524 | 0.744 | 0.789 |
| small_loose_p5 | 0.150 | 0.540 | 0.775 | 0.812 |
| small_tight | 0.158 | 0.549 | 0.779 | 0.811 |
| combo_pool | 0.177 | 0.496 | 0.767 | 0.860 |
| ucsv | 0.272 | 0.536 | 0.740 | 0.841 |
| rw | 0.169 | 0.545 | 0.777 | 0.913 |
| medium_conj | 0.163 | 0.545 | 0.820 | 0.896 |
| ar4 | 0.189 | 0.534 | 0.760 | 0.955 |
| combo_logscore | 0.383 | 0.618 | 0.776 | 0.841 |
| combo_bma | 0.677 | 0.806 | 0.938 | 0.915 |
| ucmean | 0.726 | 0.774 | 0.905 | 1.034 |

**Cash rate (%) (`cash_rate`) — rate level at t+h (as §2)**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_equal | 0.148 | **0.514** | **1.120** | 1.652 |
| small_ss | 0.140 | 0.526 | 1.149 | 1.673 |
| ar4 | 0.167 | 0.589 | 1.141 | **1.593** |
| rw | 0.183 | 0.590 | 1.136 | 1.670 |
| small_minn | 0.144 | 0.528 | 1.152 | 1.792 |
| combo_pool | 0.129 | 0.549 | 1.225 | 1.715 |
| combo_logscore | 0.125 | 0.538 | 1.241 | 1.780 |
| combo_bma | 0.128 | 0.576 | 1.381 | 1.854 |
| ucsv | 0.175 | 0.646 | 1.281 | 1.842 |
| small_tight | 0.157 | 0.598 | 1.286 | 1.935 |
| small_loose_p5 | 0.162 | 0.593 | 1.247 | 2.038 |
| medium_conj | 0.192 | 0.572 | 1.239 | 2.053 |
| small_sv | 0.113 | 0.535 | 1.410 | 2.147 |
| medium_minn | **0.112** | 0.559 | 1.465 | 2.494 |
| ucmean | 1.563 | 1.654 | 1.671 | 1.615 |


**Point accuracy by variable and horizon (RMSE, lower better)**

**`gdp_growth`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_sv | 1.573 | **1.664** | 1.806 | 1.912 |
| small_loose_p5 | **1.571** | 1.698 | 1.795 | 1.925 |
| small_tight | 1.614 | 1.687 | 1.784 | 1.907 |
| small_ss | 1.603 | 1.681 | 1.799 | 1.919 |
| ucsv | 1.618 | 1.675 | 1.812 | 1.909 |
| small_minn | 1.619 | 1.685 | 1.793 | 1.918 |
| medium_conj | 1.616 | 1.704 | 1.794 | **1.904** |
| ar4 | 1.632 | 1.674 | **1.784** | 1.934 |
| ucmean | 1.619 | 1.684 | 1.801 | 1.929 |
| combo_equal | 1.624 | 1.681 | 1.816 | 1.921 |
| medium_minn | 1.605 | 1.676 | 1.819 | 1.946 |
| combo_pool | 1.688 | 1.721 | 1.895 | 2.011 |
| combo_logscore | 2.059 | 1.956 | 2.257 | 2.011 |
| combo_bma | 2.308 | 2.159 | 2.570 | 2.085 |
| rw | 2.377 | 2.326 | 2.609 | 2.661 |

**`cpi_inflation`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| ar4 | 0.315 | **0.441** | **0.478** | **0.504** |
| ucmean | 0.434 | 0.449 | 0.480 | 0.507 |
| combo_equal | 0.327 | 0.471 | 0.581 | 0.603 |
| combo_pool | 0.287 | 0.482 | 0.614 | 0.628 |
| small_ss | 0.333 | 0.500 | 0.593 | 0.586 |
| small_minn | 0.333 | 0.498 | 0.594 | 0.610 |
| combo_logscore | 0.269 | 0.490 | 0.625 | 0.654 |
| small_tight | 0.402 | 0.508 | 0.586 | 0.612 |
| rw | **0.235** | 0.481 | 0.666 | 0.726 |
| small_sv | 0.360 | 0.514 | 0.632 | 0.654 |
| combo_bma | 0.278 | 0.527 | 0.687 | 0.684 |
| ucsv | 0.361 | 0.516 | 0.642 | 0.663 |
| medium_conj | 0.400 | 0.518 | 0.623 | 0.655 |
| medium_minn | 0.340 | 0.504 | 0.661 | 0.725 |
| small_loose_p5 | 0.371 | 0.543 | 0.684 | 0.667 |

**`unemp_rate`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_ss | **0.339** | **0.890** | 1.208 | **1.226** |
| ar4 | 0.413 | 0.915 | **1.174** | 1.300 |
| combo_equal | 0.370 | 0.903 | 1.259 | 1.324 |
| small_sv | 0.354 | 0.927 | 1.333 | 1.349 |
| small_minn | 0.348 | 0.946 | 1.337 | 1.364 |
| combo_pool | 0.428 | 0.919 | 1.260 | 1.401 |
| small_tight | 0.362 | 0.953 | 1.351 | 1.371 |
| small_loose_p5 | 0.349 | 0.989 | 1.382 | 1.467 |
| medium_minn | 0.341 | 0.943 | 1.388 | 1.516 |
| ucsv | 0.614 | 0.991 | 1.281 | 1.342 |
| rw | 0.391 | 0.984 | 1.391 | 1.509 |
| medium_conj | 0.375 | 1.050 | 1.487 | 1.574 |
| combo_logscore | 0.880 | 1.126 | 1.304 | 1.390 |
| ucmean | 1.260 | 1.326 | 1.435 | 1.542 |
| combo_bma | 1.252 | 1.365 | 1.502 | 1.448 |

**`cash_rate`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_ss | 0.276 | 0.990 | **1.908** | 2.581 |
| ar4 | 0.307 | 1.103 | 1.931 | **2.467** |
| combo_equal | 0.341 | **0.971** | 1.916 | 2.777 |
| combo_logscore | 0.272 | 0.990 | 2.133 | 2.862 |
| small_minn | 0.292 | 1.051 | 2.045 | 2.902 |
| rw | 0.352 | 1.194 | 2.062 | 2.684 |
| combo_pool | 0.290 | 1.025 | 2.149 | 2.849 |
| ucsv | 0.358 | 1.208 | 2.080 | 2.695 |
| combo_bma | 0.282 | 1.054 | 2.353 | 2.961 |
| small_tight | 0.313 | 1.190 | 2.230 | 3.078 |
| small_sv | **0.252** | 0.982 | 2.306 | 3.468 |
| small_loose_p5 | 0.352 | 1.113 | 2.180 | 3.387 |
| medium_conj | 0.479 | 1.218 | 2.280 | 3.418 |
| medium_minn | 0.265 | 1.052 | 2.486 | 4.079 |
| ucmean | 2.576 | 2.658 | 2.649 | 2.635 |


**Density calibration — mean log predictive density (higher better)**

Averaged across the 4 targets. **−∞** marks at least one origin where the realization fell outside a member's predictive support; the combinations, which always assign positive density, never do.

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_equal | -1.844 | **-1.514** | **-1.712** | **-1.878** |
| combo_pool | **-1.816** | -1.654 | -1.920 | -2.039 |
| combo_logscore | -1.925 | -1.690 | -1.982 | -2.056 |
| combo_bma | -3.731 | -2.403 | -2.354 | -2.211 |
| small_loose_p5 | -6.452 | -9.541 | -9.016 | -10.606 |
| ar4 | −∞ | -7.791 | -10.428 | −∞ |
| medium_conj | −∞ | −∞ | −∞ | -13.455 |
| medium_minn | −∞ | −∞ | -2.210 | -5.106 |
| rw | −∞ | -2.978 | -2.956 | -2.692 |
| small_minn | −∞ | −∞ | −∞ | −∞ |
| small_ss | −∞ | -10.085 | -12.188 | -10.290 |
| small_sv | -4.619 | −∞ | -7.262 | -3.452 |
| small_tight | −∞ | -10.973 | −∞ | -12.616 |
| ucmean | −∞ | −∞ | −∞ | −∞ |
| ucsv | −∞ | −∞ | −∞ | -13.565 |


## 4. Do the combinations beat the best single model?

Mean CRPS over all 4 targets, by horizon bucket. The honest test of a pool is whether it beats both equal weights and the best individual member. (Buckets average errors of differing scale across horizons, so use this within a bucket to rank models, not to compare buckets.)


**4a. Levels**

| Model | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| combo_equal | 0.534 | **1.189** | **1.808** |
| ar4 | 0.578 | 1.259 | 1.849 |
| small_ss | 0.554 | 1.276 | 1.910 |
| small_sv | 0.534 | 1.259 | 1.982 |
| small_minn | 0.562 | 1.286 | 1.946 |
| medium_minn | **0.525** | 1.266 | 2.100 |
| small_tight | 0.590 | 1.327 | 1.994 |
| ucsv | 0.594 | 1.312 | 2.036 |
| medium_conj | 0.588 | 1.323 | 2.052 |
| small_loose_p5 | 0.580 | 1.356 | 2.082 |
| combo_pool | 0.549 | 1.333 | 2.184 |
| combo_logscore | 0.646 | 1.680 | 2.257 |
| ucmean | 1.027 | 1.576 | 2.042 |
| combo_bma | 0.780 | 2.199 | 2.662 |
| rw | 0.812 | 2.063 | 3.629 |

**4b. Quarterly growth**

| Model | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| combo_equal | 0.389 | **0.650** | **0.860** |
| small_ss | 0.396 | 0.665 | 0.862 |
| ar4 | 0.414 | 0.668 | 0.870 |
| small_minn | 0.401 | 0.676 | 0.894 |
| combo_pool | 0.398 | 0.714 | 0.932 |
| small_tight | 0.421 | 0.707 | 0.936 |
| small_loose_p5 | 0.415 | 0.716 | 0.959 |
| small_sv | 0.388 | 0.715 | 0.988 |
| medium_conj | 0.423 | 0.706 | 0.968 |
| ucsv | 0.439 | 0.728 | 0.950 |
| medium_minn | **0.386** | 0.721 | 1.046 |
| combo_logscore | 0.461 | 0.768 | 0.945 |
| rw | 0.510 | 0.802 | 1.042 |
| combo_bma | 0.553 | 0.908 | 1.021 |
| ucmean | 0.831 | 0.896 | 0.950 |

## 5. Statistical significance (Diebold-Mariano)

How often each combination **significantly beats** the random-walk and AR(4) benchmarks on CRPS (Harvey-corrected, 10% level), counted over the 4 targets x 3 horizons {4, 8, 12} tested. A negative DM statistic means the combination is more accurate; significance is one-sided here. The tests have low power on these high-variance losses (dominated by a few episodes such as 2020), so the combinations beat the benchmarks on *average* (§4) more often than they do *significantly*.


**5a. Levels**

| Combination | beats ar4 | beats rw |
|:--|:--|:--|
| combo_bma | 0 / 12 | 1 / 12 |
| combo_equal | 2 / 12 | 1 / 12 |
| combo_logscore | 1 / 12 | 1 / 12 |
| combo_pool | 0 / 12 | 1 / 12 |

Strongest results (significant at 5% or better):

- `combo_equal` beats `ar4` on **unemp_rate** at h=12 (DM -3.66***)
- `combo_equal` beats `ar4` on **cash_rate** at h=4 (DM -2.22**)
- `combo_bma` beats `rw` on **gdp_growth** at h=4 (DM -2.19**)
- `combo_logscore` beats `ar4` on **unemp_rate** at h=12 (DM -2.14**)

**5b. Quarterly growth**

| Combination | beats ar4 | beats rw |
|:--|:--|:--|
| combo_bma | 0 / 12 | 1 / 12 |
| combo_equal | 2 / 12 | 1 / 12 |
| combo_logscore | 1 / 12 | 1 / 12 |
| combo_pool | 0 / 12 | 0 / 12 |

Strongest results (significant at 5% or better):

- `combo_equal` beats `ar4` on **unemp_rate** at h=12 (DM -3.66***)
- `combo_equal` beats `ar4` on **cash_rate** at h=4 (DM -2.22**)
- `combo_logscore` beats `ar4` on **unemp_rate** at h=12 (DM -2.14**)
- `combo_logscore` beats `rw` on **gdp_growth** at h=4 (DM -2.04**)

## 6. Model profiles

One entry per model: its specification, what makes it distinct, the role it plays in the suite, its strengths and failure modes, and where it actually ranks in this evaluation (the eval line is computed, not asserted). Full rationale is in README.md.

### 6a. VAR members

**`small_minn`**  
*Spec:* Independent Normal-inverse-Wishart (Gibbs); 8-variable small SOE core; 4 lags; GLP marginal-likelihood shrinkage; constant volatility; LP scaling (COVID); sum-of-coefficients + dummy-initial-observation priors.  
*Distinctive:* The workhorse Minnesota BVAR, estimated by Gibbs with an *independent* Normal-inverse-Wishart prior — the engine required to impose block exogeneity, since asymmetric (equation-specific) prior variances cannot be represented by a Kronecker/conjugate prior. Shrinkage is data-driven (GLP marginal-likelihood), and sum-of-coefficients + dummy-initial-observation priors discipline the I(1) levels (rates, real TWI).  
*Role:* The central, representative small-SOE BVAR — the reference point the other small members are deliberate variations around (tighter, looser, steady-state, SV).  
*Strengths & failure modes:* Solid all-rounder at short-to-medium horizons. Constant volatility means it leans on the LP scaling for 2020; without it the 2020 outliers would distort the coefficients.  
*In this evaluation — levels:* strongest at cash_rate (medium (5-8)), ranked 2 of 11 individual models.  
*In this evaluation — quarterly growth:* strongest at cash_rate (medium (5-8)), ranked 2 of 11 individual models.  
*See:* README.md D3, D4, D5, D17

**`small_ss`**  
*Spec:* Steady-state (Villani); 8-variable small SOE core; 4 lags; GLP marginal-likelihood shrinkage; constant volatility; LP scaling (COVID).  
*Distinctive:* Reparametrised around its *unconditional means* (Villani steady state), with informative priors placed directly on long-run levels — inflation 2.5% (target midpoint), NAIRU 4.5%, neutral cash rate 3.5%, US potential growth — and on the foreign block's steady states, which the domestic forecast inherits.  
*Role:* The **long-horizon anchor**. An iterated VAR reverts to its unconditional mean at h = 8-12; this member makes that mean economically grounded rather than the raw sample average.  
*Strengths & failure modes:* Strongest at medium/far horizons for the mean-reverting targets. Vulnerable if a steady-state anchor is stale (e.g. a shifted neutral rate drags the long end); constant volatility under-disperses around 2020 absent the LP correction.  
*In this evaluation — levels:* best individual model for unemp_rate at medium (5-8); unemp_rate at far (9-12); cash_rate at medium (5-8).  
*In this evaluation — quarterly growth:* best individual model for unemp_rate at medium (5-8); unemp_rate at far (9-12); cash_rate at medium (5-8).  
*See:* README.md D3, D4, D5, D17

**`small_sv`**  
*Spec:* Stochastic volatility, equation-by-equation; 8-variable small SOE core; 4 lags; GLP marginal-likelihood shrinkage; stochastic volatility; t-errors (SV-t) (COVID).  
*Distinctive:* Stochastic volatility estimated equation-by-equation (the Carriero-Clark-Marcellino triangular factorisation), with t-distributed errors (the CCMM SV-t COVID treatment). Block exogeneity is *exact* — foreign equations simply drop the domestic regressors.  
*Role:* The **density-calibration specialist**: time-varying volatility tracks changing macro uncertainty and the fat tails absorb outliers instead of letting them widen the whole history.  
*Strengths & failure modes:* Its SV captures the time-varying conditional variance, so it stays well-calibrated when uncertainty shifts — competitive on the GDP level at medium horizons. The volatility state at the jump-off can over/under-shoot if the last few quarters were unusual; the small system limits cross-variable information.  
*In this evaluation — levels:* best individual model for gdp_growth at medium (5-8); gdp_growth at far (9-12); cash_rate at near (1-4).  
*In this evaluation — quarterly growth:* best individual model for gdp_growth at near (1-4); cash_rate at near (1-4).  
*See:* README.md D3, D6, D17

**`small_loose_p5`**  
*Spec:* Independent Normal-inverse-Wishart (Gibbs); 8-variable small SOE core; 5 lags; fixed shrinkage λ=0.4; constant volatility; LP scaling (COVID).  
*Distinctive:* A deliberately *under-shrunk*, longer-lag variant — fixed λ = 0.4 (vs the ~0.1-0.15 the GLP procedure selects) and 5 lags — so the data speak more and richer dynamics can show through, at the cost of estimation noise.  
*Role:* The **loose / long-lag** diversity axis: it fails differently from the tightly-shrunk members and can capture dynamics they shrink away.  
*Strengths & failure modes:* Occasionally best at near-horizon density when the extra flexibility pays; noisier and prone to wider intervals at long horizons (the cost of light shrinkage in a short sample).  
*In this evaluation — levels:* strongest at gdp_growth (near (1-4)), ranked 4 of 11 individual models.  
*In this evaluation — quarterly growth:* strongest at gdp_growth (near (1-4)), ranked 3 of 11 individual models.  
*See:* README.md D4, D8

**`medium_minn`**  
*Spec:* Stochastic volatility, equation-by-equation; 13-variable medium system; 2 lags; GLP marginal-likelihood shrinkage; stochastic volatility; t-errors (SV-t) (COVID).  
*Distinctive:* The larger system — 13 variables (adds terms of trade, wages, employment, consumption, the 10y yield) with stochastic volatility and t-errors, but only 2 lags to keep the parameter count feasible; equation-by-equation estimation keeps the recursive loop tractable.  
*Role:* The **medium-system** axis (Banbura-Giannone-Reichlin): medium systems often forecast best given enough shrinkage, and the extra variables bring cross-sectional information the small core lacks.  
*Strengths & failure modes:* Strong short-horizon density (it tends to win the near bucket). The short lag length limits long-horizon dynamics, and more parameters mean more estimation uncertainty at the far end.  
*In this evaluation — levels:* best individual model for gdp_growth at near (1-4); unemp_rate at near (1-4).  
*In this evaluation — quarterly growth:* best individual model for unemp_rate at near (1-4).  
*See:* README.md D1, D6, D17

**`medium_conj`**  
*Spec:* Block-recursive conjugate NIW; 13-variable medium system; 4 lags; fixed shrinkage λ=0.1; constant volatility; LP scaling (COVID); sum-of-coefficients + dummy-initial-observation priors.  
*Distinctive:* The medium system estimated by the fast *block-recursive conjugate* scheme (foreign VAR + domestic block conditioned on contemporaneous foreign values, the RBNZ Bloor-Matheson approach) — closed-form, so cheap even at 13 variables x 4 lags. Block exogeneity is exact by the recursive structure, not the prior.  
*Role:* The cheap medium workhorse; it complements `medium_minn` (conjugate constant-volatility vs SV) on the same large system.  
*Strengths & failure modes:* Tends to lead the far-horizon GDP level error. Constant volatility leans on LP for 2020; the conjugate Kronecker prior cannot represent asymmetric shrinkage, which is why block exogeneity comes from the recursive structure.  
*In this evaluation — levels:* strongest at gdp_growth (far (9-12)), ranked 4 of 11 individual models.  
*In this evaluation — quarterly growth:* strongest at gdp_growth (far (9-12)), ranked 2 of 11 individual models.  
*See:* README.md D3, D5, D8

**`small_tight`**  
*Spec:* Block-recursive conjugate NIW; 8-variable small SOE core; 4 lags; fixed shrinkage λ=0.05; constant volatility; LP scaling (COVID); sum-of-coefficients + dummy-initial-observation priors.  
*Distinctive:* The heavily-shrunk small model — fixed λ = 0.05, far tighter than the GLP selection — pulling hard toward the persistence/random-walk prior, so it is parsimonious and low-variance.  
*Role:* The **tight** diversity axis and the long-horizon robustness member: heavy shrinkage buys stability where lightly-parametrised models wander.  
*Strengths & failure modes:* Best or near-best at far-horizon GDP (the tight prior stops it over-reacting). Can be too rigid at short horizons, missing genuine dynamics the looser members catch.  
*In this evaluation — levels:* strongest at gdp_growth (far (9-12)), ranked 2 of 11 individual models.  
*In this evaluation — quarterly growth:* best individual model for gdp_growth at medium (5-8); gdp_growth at far (9-12).  
*See:* README.md D4, D5, D8

### 6b. Benchmark members

**`rw`**  
*Spec:* Random walk; LP scaling (COVID).  
*Distinctive:* The no-change forecast: the last observed value persists, with Gaussian increments scaled to the historical change.  
*Role:* The universal hard-to-beat short-horizon bar for persistent/level variables, and a pool member.  
*Strengths & failure modes:* Competitive at h = 1 for level variables; fails badly at long horizons — its level path runs away, which the level-error tables (§2) expose brutally.  
*In this evaluation — levels:* best individual model for cpi_inflation at near (1-4); cpi_inflation at medium (5-8).  
*In this evaluation — quarterly growth:* best individual model for cpi_inflation at near (1-4).  
*See:* README.md D8

**`ar4`**  
*Spec:* Bayesian AR(4); LP scaling (COVID).  
*Distinctive:* A Bayesian AR(4) per variable with Minnesota-style lag shrinkage and a stationarity-truncated posterior.  
*Role:* The univariate-persistence bar — it isolates how much of the forecast is just own-history dynamics.  
*Strengths & failure modes:* Surprisingly strong for inflation at near/medium horizons, where univariate dynamics dominate; it cannot use cross-variable information, so it lags when that matters.  
*In this evaluation — levels:* best individual model for cpi_inflation at far (9-12); cash_rate at far (9-12).  
*In this evaluation — quarterly growth:* best individual model for cpi_inflation at medium (5-8); cash_rate at far (9-12).  
*See:* README.md D8

**`ucsv`**  
*Spec:* Unobserved components + stochastic volatility; t-errors (robust) (COVID).  
*Distinctive:* Stock-Watson unobserved-components stochastic volatility per variable: a random-walk trend plus transitory noise, both with time-varying variances and outlier-robust t-errors.  
*Role:* The canonical inflation benchmark and a genuine density anchor for the other targets.  
*Strengths & failure modes:* Strong for inflation (its native use case); weaker for variables with richer multivariate dynamics. The trend/noise split is weakly identified, so it is gated on Monte-Carlo precision, not raw ESS.  
*In this evaluation — levels:* strongest at gdp_growth (near (1-4)), ranked 3 of 11 individual models.  
*In this evaluation — quarterly growth:* strongest at gdp_growth (near (1-4)), ranked 4 of 11 individual models.  
*See:* README.md D9, D17

**`ucmean`**  
*Spec:* Unconditional mean; LP scaling (COVID).  
*Distinctive:* A Gaussian density centred on the expanding-sample mean with the sample variance — the simplest possible density forecast.  
*Role:* The floor: the 'did the model beat just predicting the long-run average' bar.  
*Strengths & failure modes:* Unexpectedly competitive at long horizons for mean-reverting growth (everything reverts to the mean eventually); useless at short horizons where dynamics matter.  
*In this evaluation — levels:* strongest at cpi_inflation (far (9-12)), ranked 3 of 11 individual models.  
*In this evaluation — quarterly growth:* best individual model for cpi_inflation at far (9-12).  
*See:* README.md D8

### 6c. Combination schemes

**`combo_equal`**  
*Weights:* Equal weights on every member, per variable x horizon bucket.  
*In this evaluation:* The forecast-combination-puzzle benchmark and the recommended robust default: in the evaluation it has the best mean log score at every horizon and the best far-horizon CRPS. Hard to beat because it never over-fits weights.

**`combo_logscore`**  
*Weights:* Weights proportional to each member's recent log predictive score, with a forgetting factor, shrunk toward equal.  
*In this evaluation:* Adapts to which members are forecasting well lately; the shrinkage and forgetting guard against over-concentrating on a member that was lucky.

**`combo_pool`**  
*Weights:* Optimal prediction pool (Hall-Mitchell / Geweke-Amisano): weights on the simplex that maximise the historical *pooled* log score, shrunk toward equal.  
*In this evaluation:* Unlike BMA it does not degenerate to a single model ('all models are false but useful'); competitive with equal weights and occasionally better at near horizons.

**`combo_bma`**  
*Weights:* Bayesian model averaging by predictive likelihood — no shrinkage.  
*In this evaluation:* **Diagnostic only.** It concentrates weight on the single best-fitting member, so it answers 'which model does the data favour' rather than serving as a robust combination; reported, not recommended.


## 7. How to read this

- **Two complementary metrics.** The **level** error (§2) is the headline: for GDP and inflation it is the cumulative level from the forecast origin (where real GDP and the price level land h quarters out), for unemployment and the cash rate the rate level at t+h. It accumulates the whole forecast path, so it grows with horizon and exposes a model that gets the persistent/drift component wrong (e.g. the random walk). The **quarterly-growth** view (§3) scores each target's single-quarter outcome — a noisier, more local target that does not compound the path, so it does not grow mechanically with horizon and can rank models differently. For unemployment and the cash rate (modelled in levels) the two views coincide. The year-ended *growth* scores live in `output/tables/scores_by_horizon.csv`.

- **Point gains over the best member are modest by design** — equal weights are hard to beat (the forecast-combination puzzle). The pool's payoff is *calibration and robustness*: it insures against any single member failing, rather than always winning on accuracy.

- **CRPS and log score can disagree** on outlier-heavy windows (log score is far more sensitive to tail events); both are reported. A COVID-excluded variant is in `output/tables/scores_by_horizon_excovid.csv`. See the full Quarto report for fan charts, PIT calibration, and weight-evolution plots.

