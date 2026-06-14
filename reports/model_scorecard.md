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
| gdp_growth | medium_minn (0.988) | small_sv (1.693) | small_sv (2.156) |
| cpi_inflation | rw (0.386) | rw (1.358) | ar4 (2.693) |
| unemp_rate | medium_minn (0.318) | small_ss (0.631) | small_ss (0.721) |
| cash_rate | small_sv (0.307) | combo_equal (0.893) | ar4 (1.447) |


**Density accuracy by variable and horizon (CRPS, lower better)**

**Real GDP growth (qtr %) (`gdp_growth`) — cumulative level from the origin**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_sv | 0.611 | **1.267** | **1.954** | 2.206 |
| medium_minn | **0.610** | 1.275 | 1.995 | 2.314 |
| combo_equal | 0.634 | 1.305 | 2.030 | 2.237 |
| small_tight | 0.682 | 1.366 | 2.101 | **2.165** |
| ucsv | 0.634 | 1.304 | 2.020 | 2.424 |
| medium_conj | 0.672 | 1.383 | 2.139 | 2.238 |
| small_ss | 0.664 | 1.377 | 2.169 | 2.240 |
| small_loose_p5 | 0.627 | 1.376 | 2.173 | 2.313 |
| ar4 | 0.668 | 1.342 | 2.156 | 2.470 |
| small_minn | 0.674 | 1.420 | 2.259 | 2.337 |
| ucmean | 0.680 | 1.380 | 2.273 | 2.662 |
| combo_pool | 0.642 | 1.425 | 2.651 | 3.685 |
| combo_logscore | 0.754 | 1.789 | 4.352 | 3.995 |
| combo_bma | 0.842 | 2.203 | 6.319 | 5.314 |
| rw | 0.941 | 3.264 | 6.711 | 11.084 |

**Trimmed-mean CPI inflation (qtr %) (`cpi_inflation`) — cumulative level from the origin**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| rw | **0.131** | **0.674** | **1.832** | 3.356 |
| ar4 | 0.162 | 0.844 | 1.973 | **3.108** |
| combo_equal | 0.161 | 0.770 | 1.919 | 3.340 |
| combo_pool | 0.147 | 0.761 | 2.043 | 3.523 |
| combo_logscore | 0.137 | 0.729 | 2.075 | 3.592 |
| small_minn | 0.170 | 0.827 | 2.113 | 3.714 |
| small_ss | 0.172 | 0.849 | 2.171 | 3.663 |
| small_sv | 0.170 | 0.849 | 2.129 | 3.761 |
| ucmean | 0.230 | 1.056 | 2.269 | 3.405 |
| medium_minn | 0.165 | 0.814 | 2.141 | 4.002 |
| combo_bma | 0.145 | 0.800 | 2.363 | 3.936 |
| small_tight | 0.198 | 0.925 | 2.241 | 3.910 |
| medium_conj | 0.191 | 0.895 | 2.234 | 4.029 |
| ucsv | 0.175 | 0.920 | 2.331 | 4.056 |
| small_loose_p5 | 0.192 | 0.935 | 2.436 | 4.117 |

**Unemployment rate (%) (`unemp_rate`) — rate level at t+h**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_ss | **0.146** | 0.500 | **0.688** | **0.747** |
| medium_minn | 0.146 | **0.466** | 0.722 | 0.796 |
| combo_equal | 0.159 | 0.490 | 0.699 | 0.782 |
| small_sv | 0.154 | 0.501 | 0.758 | 0.782 |
| small_minn | 0.152 | 0.529 | 0.753 | 0.787 |
| small_loose_p5 | 0.155 | 0.546 | 0.780 | 0.808 |
| small_tight | 0.158 | 0.549 | 0.779 | 0.811 |
| combo_pool | 0.174 | 0.490 | 0.765 | 0.871 |
| ucsv | 0.265 | 0.531 | 0.737 | 0.832 |
| combo_logscore | 0.232 | 0.532 | 0.774 | 0.848 |
| rw | 0.169 | 0.545 | 0.777 | 0.913 |
| medium_conj | 0.163 | 0.545 | 0.820 | 0.896 |
| ar4 | 0.189 | 0.534 | 0.760 | 0.955 |
| combo_bma | 0.148 | 0.475 | 0.935 | 0.886 |
| ucmean | 0.726 | 0.774 | 0.905 | 1.034 |

**Cash rate (%) (`cash_rate`) — rate level at t+h**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_equal | 0.147 | **0.514** | **1.126** | 1.668 |
| ar4 | 0.167 | 0.589 | 1.141 | **1.593** |
| small_ss | 0.139 | 0.525 | 1.151 | 1.676 |
| rw | 0.183 | 0.590 | 1.136 | 1.670 |
| small_minn | 0.145 | 0.527 | 1.152 | 1.785 |
| combo_pool | 0.130 | 0.547 | 1.223 | 1.724 |
| combo_logscore | 0.122 | 0.532 | 1.223 | 1.758 |
| ucsv | 0.170 | 0.639 | 1.266 | 1.832 |
| combo_bma | 0.124 | 0.576 | 1.393 | 1.850 |
| small_tight | 0.157 | 0.598 | 1.286 | 1.935 |
| medium_conj | 0.192 | 0.572 | 1.239 | 2.053 |
| small_loose_p5 | 0.162 | 0.609 | 1.276 | 2.076 |
| small_sv | 0.115 | 0.533 | 1.422 | 2.111 |
| medium_minn | **0.113** | 0.559 | 1.462 | 2.467 |
| ucmean | 1.563 | 1.654 | 1.671 | 1.615 |


**Point accuracy by variable and horizon (RMSE, lower better)**

**`gdp_growth`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_sv | 1.583 | **2.554** | 3.713 | **3.739** |
| small_tight | 1.614 | 2.623 | 3.654 | 3.830 |
| small_ss | 1.603 | 2.599 | 3.682 | 3.845 |
| ucsv | 1.599 | 2.556 | **3.620** | 4.045 |
| small_minn | 1.612 | 2.636 | 3.735 | 3.938 |
| ar4 | 1.632 | 2.632 | 3.631 | 4.036 |
| small_loose_p5 | **1.566** | 2.608 | 3.691 | 4.079 |
| ucmean | 1.619 | 2.634 | 3.685 | 4.152 |
| medium_conj | 1.616 | 2.690 | 3.774 | 4.030 |
| medium_minn | 1.578 | 2.644 | 3.803 | 4.319 |
| combo_equal | 1.620 | 2.818 | 4.130 | 4.938 |
| combo_pool | 1.688 | 3.361 | 6.337 | 8.822 |
| combo_logscore | 2.054 | 5.603 | 11.629 | 8.818 |
| combo_bma | 2.306 | 6.971 | 15.316 | 11.126 |
| rw | 2.377 | 7.807 | 15.379 | 24.092 |

**`cpi_inflation`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| ar4 | 0.315 | 1.433 | **2.905** | **4.138** |
| ucmean | 0.434 | 1.648 | 3.084 | 4.269 |
| combo_equal | 0.326 | 1.501 | 3.322 | 5.216 |
| combo_pool | 0.285 | 1.443 | 3.361 | 5.454 |
| combo_logscore | 0.261 | 1.388 | 3.452 | 5.668 |
| small_ss | 0.332 | 1.589 | 3.549 | 5.427 |
| small_minn | 0.333 | 1.587 | 3.543 | 5.511 |
| rw | **0.235** | **1.301** | 3.486 | 6.042 |
| small_tight | 0.402 | 1.736 | 3.655 | 5.611 |
| small_sv | 0.357 | 1.663 | 3.728 | 5.929 |
| combo_bma | 0.277 | 1.520 | 3.860 | 6.162 |
| medium_conj | 0.400 | 1.750 | 3.777 | 5.951 |
| ucsv | 0.357 | 1.694 | 3.819 | 6.028 |
| medium_minn | 0.341 | 1.599 | 3.722 | 6.257 |
| small_loose_p5 | 0.370 | 1.732 | 3.932 | 6.210 |

**`unemp_rate`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_ss | 0.337 | **0.891** | 1.205 | **1.223** |
| ar4 | 0.413 | 0.915 | **1.174** | 1.300 |
| combo_equal | 0.369 | 0.902 | 1.262 | 1.321 |
| small_minn | 0.349 | 0.948 | 1.338 | 1.366 |
| combo_pool | 0.426 | 0.919 | 1.260 | 1.400 |
| small_tight | 0.362 | 0.953 | 1.351 | 1.371 |
| small_sv | 0.346 | 0.946 | 1.395 | 1.360 |
| medium_minn | **0.336** | 0.918 | 1.377 | 1.466 |
| ucsv | 0.606 | 0.974 | 1.259 | 1.334 |
| combo_bma | 0.339 | 0.940 | 1.507 | 1.416 |
| small_loose_p5 | 0.354 | 1.005 | 1.395 | 1.470 |
| rw | 0.391 | 0.984 | 1.391 | 1.509 |
| combo_logscore | 0.583 | 1.047 | 1.307 | 1.385 |
| medium_conj | 0.375 | 1.050 | 1.487 | 1.574 |
| ucmean | 1.260 | 1.326 | 1.435 | 1.542 |

**`cash_rate`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_ss | 0.278 | 0.990 | **1.908** | 2.579 |
| ar4 | 0.307 | 1.103 | 1.931 | **2.467** |
| combo_equal | 0.342 | **0.969** | 1.923 | 2.783 |
| combo_logscore | 0.264 | 0.975 | 2.141 | 2.869 |
| small_minn | 0.292 | 1.045 | 2.033 | 2.895 |
| ucsv | 0.354 | 1.190 | 2.058 | 2.684 |
| rw | 0.352 | 1.194 | 2.062 | 2.684 |
| combo_pool | 0.279 | 1.008 | 2.161 | 2.856 |
| combo_bma | 0.282 | 1.049 | 2.397 | 2.967 |
| small_tight | 0.313 | 1.190 | 2.230 | 3.078 |
| small_sv | 0.256 | 0.978 | 2.324 | 3.425 |
| small_loose_p5 | 0.355 | 1.130 | 2.221 | 3.433 |
| medium_conj | 0.479 | 1.218 | 2.280 | 3.418 |
| medium_minn | **0.248** | 1.033 | 2.525 | 4.098 |
| ucmean | 2.576 | 2.658 | 2.649 | 2.635 |


**Density calibration — mean log predictive density (higher better)**

Averaged across the 4 targets. **−∞** marks at least one origin where the realization fell outside a member's predictive support; the combinations, which always assign positive density, never do.

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_equal | -1.852 | **-1.573** | **-2.325** | **-2.751** |
| combo_logscore | -1.840 | -1.640 | -2.567 | -2.929 |
| combo_pool | **-1.837** | -1.686 | -2.536 | -2.926 |
| combo_bma | -2.056 | -1.974 | -3.147 | -3.196 |
| medium_minn | -3.127 | -2.740 | -3.785 | -3.335 |
| small_sv | -4.600 | -3.141 | -4.015 | -3.512 |
| ucsv | -5.404 | -6.356 | -8.782 | -7.300 |
| small_loose_p5 | -6.527 | -6.058 | -9.956 | -6.763 |
| small_minn | -6.573 | -6.665 | -11.009 | -10.964 |
| ar4 | −∞ | -4.390 | -5.271 | -6.282 |
| medium_conj | −∞ | -4.830 | -10.309 | -10.904 |
| rw | −∞ | -4.619 | -4.431 | -4.180 |
| small_ss | −∞ | -6.440 | -10.932 | -8.827 |
| small_tight | −∞ | -9.674 | -15.438 | -14.541 |
| ucmean | −∞ | -5.804 | -6.761 | -10.352 |


## 3. Forecast performance — quarterly growth

The same scoring on each target's **single-quarter outcome** at t+h: the quarterly growth rate for real GDP and CPI inflation, and — since they are modelled in levels — the rate level for unemployment and the cash rate (identical to §2 for those two). Quarterly growth is a noisier, more local target than the accumulated level: it does not compound the path, so errors do not grow mechanically with horizon and the rankings can differ from the level view. Same conventions as §2.

**Best single model by variable and horizon bucket (CRPS)** — lowest mean CRPS in the bucket; value in parentheses:

| Variable | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| gdp_growth | small_sv (0.683) | small_tight (0.773) | small_tight (0.846) |
| cpi_inflation | rw (0.189) | ar4 (0.273) | ucmean (0.302) |
| unemp_rate | medium_minn (0.318) | small_ss (0.631) | small_ss (0.721) |
| cash_rate | small_sv (0.307) | combo_equal (0.893) | ar4 (1.447) |


**Density accuracy by variable and horizon (CRPS, lower better)**

**Real GDP growth (qtr %) (`gdp_growth`) — quarterly growth rate (%)**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_loose_p5 | 0.627 | 0.741 | 0.808 | 0.883 |
| combo_equal | 0.634 | 0.734 | 0.817 | 0.887 |
| ar4 | 0.668 | **0.727** | 0.797 | 0.898 |
| small_tight | 0.682 | 0.743 | **0.794** | 0.878 |
| small_sv | 0.611 | 0.735 | 0.850 | 0.901 |
| small_ss | 0.664 | 0.743 | 0.804 | 0.888 |
| small_minn | 0.674 | 0.744 | 0.802 | 0.880 |
| medium_conj | 0.672 | 0.757 | 0.800 | **0.876** |
| medium_minn | **0.610** | 0.743 | 0.840 | 0.925 |
| ucmean | 0.680 | 0.732 | 0.816 | 0.908 |
| ucsv | 0.634 | 0.767 | 0.900 | 0.981 |
| combo_pool | 0.642 | 0.758 | 0.931 | 1.018 |
| combo_logscore | 0.754 | 0.860 | 1.134 | 1.042 |
| combo_bma | 0.842 | 0.946 | 1.409 | 1.164 |
| rw | 0.941 | 1.164 | 1.387 | 1.447 |

**Trimmed-mean CPI inflation (qtr %) (`cpi_inflation`) — quarterly growth rate (%)**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| ar4 | 0.162 | **0.246** | 0.288 | 0.314 |
| ucmean | 0.230 | 0.250 | **0.286** | **0.312** |
| combo_equal | 0.161 | 0.246 | 0.323 | 0.364 |
| combo_pool | 0.147 | 0.257 | 0.354 | 0.391 |
| small_ss | 0.172 | 0.269 | 0.356 | 0.371 |
| rw | **0.131** | 0.248 | 0.361 | 0.431 |
| combo_logscore | 0.137 | 0.263 | 0.365 | 0.408 |
| small_minn | 0.170 | 0.262 | 0.356 | 0.387 |
| small_sv | 0.170 | 0.269 | 0.356 | 0.407 |
| small_tight | 0.198 | 0.272 | 0.350 | 0.395 |
| medium_minn | 0.165 | 0.264 | 0.364 | 0.444 |
| medium_conj | 0.191 | 0.272 | 0.367 | 0.426 |
| combo_bma | 0.145 | 0.288 | 0.401 | 0.439 |
| ucsv | 0.175 | 0.283 | 0.387 | 0.430 |
| small_loose_p5 | 0.192 | 0.295 | 0.415 | 0.422 |

**Unemployment rate (%) (`unemp_rate`) — rate level at t+h (as §2)**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_ss | **0.146** | 0.500 | **0.688** | **0.747** |
| medium_minn | 0.146 | **0.466** | 0.722 | 0.796 |
| combo_equal | 0.159 | 0.490 | 0.699 | 0.782 |
| small_sv | 0.154 | 0.501 | 0.758 | 0.782 |
| small_minn | 0.152 | 0.529 | 0.753 | 0.787 |
| small_loose_p5 | 0.155 | 0.546 | 0.780 | 0.808 |
| small_tight | 0.158 | 0.549 | 0.779 | 0.811 |
| combo_pool | 0.174 | 0.490 | 0.765 | 0.871 |
| ucsv | 0.265 | 0.531 | 0.737 | 0.832 |
| combo_logscore | 0.232 | 0.532 | 0.774 | 0.848 |
| rw | 0.169 | 0.545 | 0.777 | 0.913 |
| medium_conj | 0.163 | 0.545 | 0.820 | 0.896 |
| ar4 | 0.189 | 0.534 | 0.760 | 0.955 |
| combo_bma | 0.148 | 0.475 | 0.935 | 0.886 |
| ucmean | 0.726 | 0.774 | 0.905 | 1.034 |

**Cash rate (%) (`cash_rate`) — rate level at t+h (as §2)**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_equal | 0.147 | **0.514** | **1.126** | 1.668 |
| ar4 | 0.167 | 0.589 | 1.141 | **1.593** |
| small_ss | 0.139 | 0.525 | 1.151 | 1.676 |
| rw | 0.183 | 0.590 | 1.136 | 1.670 |
| small_minn | 0.145 | 0.527 | 1.152 | 1.785 |
| combo_pool | 0.130 | 0.547 | 1.223 | 1.724 |
| combo_logscore | 0.122 | 0.532 | 1.223 | 1.758 |
| ucsv | 0.170 | 0.639 | 1.266 | 1.832 |
| combo_bma | 0.124 | 0.576 | 1.393 | 1.850 |
| small_tight | 0.157 | 0.598 | 1.286 | 1.935 |
| medium_conj | 0.192 | 0.572 | 1.239 | 2.053 |
| small_loose_p5 | 0.162 | 0.609 | 1.276 | 2.076 |
| small_sv | 0.115 | 0.533 | 1.422 | 2.111 |
| medium_minn | **0.113** | 0.559 | 1.462 | 2.467 |
| ucmean | 1.563 | 1.654 | 1.671 | 1.615 |


**Point accuracy by variable and horizon (RMSE, lower better)**

**`gdp_growth`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_sv | 1.583 | 1.700 | 1.812 | **1.893** |
| small_tight | 1.614 | 1.687 | 1.784 | 1.907 |
| small_loose_p5 | **1.566** | 1.702 | 1.799 | 1.925 |
| small_minn | 1.612 | 1.683 | 1.793 | 1.915 |
| small_ss | 1.603 | 1.684 | 1.797 | 1.922 |
| medium_conj | 1.616 | 1.704 | 1.794 | 1.904 |
| medium_minn | 1.578 | 1.701 | 1.810 | 1.930 |
| ar4 | 1.632 | **1.674** | **1.784** | 1.934 |
| ucsv | 1.599 | 1.693 | 1.812 | 1.922 |
| ucmean | 1.619 | 1.684 | 1.801 | 1.929 |
| combo_equal | 1.620 | 1.689 | 1.814 | 1.922 |
| combo_pool | 1.688 | 1.724 | 1.921 | 2.012 |
| combo_logscore | 2.054 | 1.970 | 2.254 | 2.011 |
| combo_bma | 2.306 | 2.182 | 2.570 | 2.084 |
| rw | 2.377 | 2.326 | 2.609 | 2.661 |

**`cpi_inflation`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| ar4 | 0.315 | **0.441** | **0.478** | **0.504** |
| ucmean | 0.434 | 0.449 | 0.480 | 0.507 |
| combo_equal | 0.326 | 0.471 | 0.581 | 0.603 |
| combo_pool | 0.285 | 0.482 | 0.602 | 0.628 |
| small_ss | 0.332 | 0.502 | 0.594 | 0.588 |
| combo_logscore | 0.261 | 0.487 | 0.625 | 0.654 |
| small_minn | 0.333 | 0.498 | 0.595 | 0.608 |
| small_tight | 0.402 | 0.508 | 0.586 | 0.612 |
| rw | **0.235** | 0.481 | 0.666 | 0.726 |
| small_sv | 0.357 | 0.516 | 0.632 | 0.651 |
| combo_bma | 0.277 | 0.526 | 0.675 | 0.686 |
| ucsv | 0.357 | 0.513 | 0.646 | 0.666 |
| medium_conj | 0.400 | 0.518 | 0.623 | 0.655 |
| medium_minn | 0.341 | 0.503 | 0.658 | 0.733 |
| small_loose_p5 | 0.370 | 0.542 | 0.687 | 0.668 |

**`unemp_rate`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_ss | 0.337 | **0.891** | 1.205 | **1.223** |
| ar4 | 0.413 | 0.915 | **1.174** | 1.300 |
| combo_equal | 0.369 | 0.902 | 1.262 | 1.321 |
| small_minn | 0.349 | 0.948 | 1.338 | 1.366 |
| combo_pool | 0.426 | 0.919 | 1.260 | 1.400 |
| small_tight | 0.362 | 0.953 | 1.351 | 1.371 |
| small_sv | 0.346 | 0.946 | 1.395 | 1.360 |
| medium_minn | **0.336** | 0.918 | 1.377 | 1.466 |
| ucsv | 0.606 | 0.974 | 1.259 | 1.334 |
| combo_bma | 0.339 | 0.940 | 1.507 | 1.416 |
| small_loose_p5 | 0.354 | 1.005 | 1.395 | 1.470 |
| rw | 0.391 | 0.984 | 1.391 | 1.509 |
| combo_logscore | 0.583 | 1.047 | 1.307 | 1.385 |
| medium_conj | 0.375 | 1.050 | 1.487 | 1.574 |
| ucmean | 1.260 | 1.326 | 1.435 | 1.542 |

**`cash_rate`**

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| small_ss | 0.278 | 0.990 | **1.908** | 2.579 |
| ar4 | 0.307 | 1.103 | 1.931 | **2.467** |
| combo_equal | 0.342 | **0.969** | 1.923 | 2.783 |
| combo_logscore | 0.264 | 0.975 | 2.141 | 2.869 |
| small_minn | 0.292 | 1.045 | 2.033 | 2.895 |
| ucsv | 0.354 | 1.190 | 2.058 | 2.684 |
| rw | 0.352 | 1.194 | 2.062 | 2.684 |
| combo_pool | 0.279 | 1.008 | 2.161 | 2.856 |
| combo_bma | 0.282 | 1.049 | 2.397 | 2.967 |
| small_tight | 0.313 | 1.190 | 2.230 | 3.078 |
| small_sv | 0.256 | 0.978 | 2.324 | 3.425 |
| small_loose_p5 | 0.355 | 1.130 | 2.221 | 3.433 |
| medium_conj | 0.479 | 1.218 | 2.280 | 3.418 |
| medium_minn | **0.248** | 1.033 | 2.525 | 4.098 |
| ucmean | 2.576 | 2.658 | 2.649 | 2.635 |


**Density calibration — mean log predictive density (higher better)**

Averaged across the 4 targets. **−∞** marks at least one origin where the realization fell outside a member's predictive support; the combinations, which always assign positive density, never do.

| Model | h=1 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|
| combo_equal | -1.852 | **-1.524** | **-1.757** | **-1.872** |
| combo_pool | **-1.837** | -1.670 | -1.954 | -2.028 |
| combo_logscore | -1.840 | -1.649 | -2.013 | -2.042 |
| combo_bma | -2.056 | -1.986 | -2.379 | -2.192 |
| medium_minn | -3.127 | -4.551 | -7.785 | -3.681 |
| small_sv | -4.600 | -6.722 | -6.114 | -3.030 |
| ucsv | -5.404 | -6.487 | -10.739 | -9.902 |
| small_loose_p5 | -6.527 | -9.128 | -9.524 | -8.624 |
| ar4 | −∞ | -7.791 | -10.428 | −∞ |
| medium_conj | −∞ | −∞ | −∞ | -13.455 |
| rw | −∞ | -2.978 | -2.956 | -2.692 |
| small_minn | -6.573 | −∞ | −∞ | −∞ |
| small_ss | −∞ | −∞ | −∞ | -12.302 |
| small_tight | −∞ | -10.973 | −∞ | -12.616 |
| ucmean | −∞ | −∞ | −∞ | −∞ |


## 4. Do the combinations beat the best single model?

Mean CRPS over all 4 targets, by horizon bucket. The honest test of a pool is whether it beats both equal weights and the best individual member. (Buckets average errors of differing scale across horizons, so use this within a bucket to rank models, not to compare buckets.)


**4a. Levels**

| Model | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| combo_equal | 0.533 | **1.189** | **1.812** |
| ar4 | 0.578 | 1.259 | 1.849 |
| small_ss | 0.554 | 1.275 | 1.903 |
| small_sv | 0.534 | 1.267 | 1.992 |
| small_minn | 0.564 | 1.294 | 1.957 |
| medium_minn | **0.525** | 1.264 | 2.094 |
| small_tight | 0.590 | 1.327 | 1.994 |
| ucsv | 0.590 | 1.304 | 2.033 |
| medium_conj | 0.588 | 1.323 | 2.052 |
| small_loose_p5 | 0.583 | 1.364 | 2.098 |
| combo_pool | 0.547 | 1.359 | 2.188 |
| combo_logscore | 0.610 | 1.687 | 2.252 |
| ucmean | 1.027 | 1.576 | 2.042 |
| combo_bma | 0.667 | 2.191 | 2.637 |
| rw | 0.812 | 2.063 | 3.629 |

**4b. Quarterly growth**

| Model | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| combo_equal | 0.389 | **0.650** | **0.860** |
| small_ss | 0.396 | 0.665 | 0.862 |
| ar4 | 0.414 | 0.668 | 0.870 |
| small_minn | 0.401 | 0.678 | 0.894 |
| combo_pool | 0.397 | 0.716 | 0.933 |
| small_tight | 0.421 | 0.707 | 0.936 |
| small_sv | 0.388 | 0.721 | 0.988 |
| medium_conj | 0.423 | 0.706 | 0.968 |
| ucsv | 0.438 | 0.725 | 0.946 |
| small_loose_p5 | 0.419 | 0.724 | 0.967 |
| combo_logscore | 0.428 | 0.765 | 0.944 |
| medium_minn | **0.386** | 0.718 | 1.042 |
| rw | 0.510 | 0.802 | 1.042 |
| combo_bma | 0.444 | 0.904 | 1.014 |
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

- `combo_equal` beats `ar4` on **unemp_rate** at h=12 (DM -3.37***)
- `combo_equal` beats `ar4` on **cash_rate** at h=4 (DM -2.18**)
- `combo_logscore` beats `ar4` on **unemp_rate** at h=12 (DM -2.10**)
- `combo_logscore` beats `rw` on **gdp_growth** at h=4 (DM -2.05**)

**5b. Quarterly growth**

| Combination | beats ar4 | beats rw |
|:--|:--|:--|
| combo_bma | 0 / 12 | 1 / 12 |
| combo_equal | 2 / 12 | 1 / 12 |
| combo_logscore | 1 / 12 | 1 / 12 |
| combo_pool | 0 / 12 | 0 / 12 |

Strongest results (significant at 5% or better):

- `combo_equal` beats `ar4` on **unemp_rate** at h=12 (DM -3.37***)
- `combo_equal` beats `ar4` on **cash_rate** at h=4 (DM -2.18**)
- `combo_logscore` beats `ar4` on **unemp_rate** at h=12 (DM -2.10**)

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
*In this evaluation — quarterly growth:* strongest at gdp_growth (near (1-4)), ranked 2 of 11 individual models.  
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
*In this evaluation — quarterly growth:* strongest at unemp_rate (medium (5-8)), ranked 3 of 11 individual models.  
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

