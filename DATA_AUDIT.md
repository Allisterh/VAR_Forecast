# DATA_AUDIT.md — full audit of data sources, 2026-06-13

Method: every configured series was pulled from its provider and its **metadata**
(series description, table, series type, unit, frequency, coverage) checked
against the intended economic concept; every dlog-transformed series was
additionally tested for residual quarterly seasonality (F-test of quarter
dummies on quarterly log changes); constructed series were round-trip tested.
Sections: verified series → corrections made → proxies and constructed data →
synthetic data inventory → simplifying assumptions.

## 1. Verified correct (provider metadata, as of this audit)

| Variable | Series | Verified metadata |
|---|---|---|
| `gdp_growth` | ABS A2304402X (5206.0 Tab 1) | "Gross domestic product: **Chain volume measures**", **Seasonally Adjusted**, $m, quarterly, 1959Q3– |
| `cons_growth` | ABS A2304081W (5206.0 Tab 2) | "Households; Final consumption expenditure", **chain volume**, **SA**, $m, quarterly |
| `tot` | ABS A2304200A (5206.0 Tab 1) | "Terms of trade: Index", **SA**, quarterly (see §5 note on residual Q1 pattern) |
| `unemp_rate` | ABS A84423050A (6202.0 Tab 1) | "Unemployment rate; Persons", **SA**, %, monthly 1978– |
| `emp_growth` | ABS A84423043C (6202.0 Tab 1) | "Employed total; Persons", **SA**, '000, monthly |
| `wpi_growth` | ABS **A2713849C** (6345.0 Tab 1) | "Total hourly rates of pay excluding bonuses; Australia; Private and Public; All industries", Quarterly Index, **SA** — *corrected in this audit, see §2* |
| `cpi_inflation` | RBA GCPIOCPMTMQP (G1) | "Consumer price index; **Trimmed mean**; Quarterly change (per cent)" — the RBA/ABS trimmed mean is computed from seasonally adjusted CPI components by construction |
| `cash_rate` | RBA FIRMMCRTD (F1) | "**Cash Rate Target**", daily |
| `rtwi` | RBA FRERTWI (F15) | "**Real** trade-weighted index … adjusted for relative consumer price levels", **quarterly native** |
| `bond10y` | RBA FCMYGBAG10D (F2) | "Australian Government 10 year bond" yield, interpolated, daily |
| `f_comm` | RBA GRCPAISDR (I2) | "Index of commodity prices; All items; **SDR**", monthly, 1982– |
| `f_rate` | DBnomics FED/H15/RIFSPFF_N.M | "Federal funds – Overnight", effective rate, monthly |
| `f_act` | DBnomics IMF/IFS/Q.US.AIP_IX | "United States – Economic Activity, Industrial Production, Index", quarterly (IFS mirrors the Fed's seasonally adjusted G.17 index; empirical seasonality test p = 0.83 — clean) |

Empirical seasonality F-tests on quarterly log changes (real data): all clean
(p > 0.05) except `wpi_growth` before the §2 correction (p ≈ 1e-10, the
Original-series fingerprint) and `tot` (p = 0.003; see §5).

## 2. Corrections made during this audit

1. **WPI was the wrong series type.** The original config used `A2603609J`,
   which is the **Original (non-seasonally-adjusted)** quarterly WPI index —
   confirmed both by ABS metadata (`series_type: Original`) and a glaring
   seasonal pattern in its quarterly log changes. Replaced with `A2713849C`
   (identical concept — total hourly rates of pay ex bonuses, Australia,
   private & public, all industries — **Seasonally Adjusted**).
2. *(Previous audits, recorded here for completeness.)* The commodity index
   was originally `GRCPAISAD` — the **A$-denominated, bulk-spot-prices**
   variant that (a) only starts in 2009, silently truncating the panel, and
   (b) embeds the endogenous AUD exchange rate, which contradicts the
   foreign-block price-taker assumption. Replaced with `GRCPAISDR` (all
   items, **SDR-denominated**, 1982–): SDR denomination keeps the foreign
   block a world price. The terms-of-trade ID `A2303731T` did not resolve and
   was replaced with the verified SA index `A2304200A`.

## 3. Proxies and constructed series (real-data mode)

These are **real published data**, but not literally the named concept:

1. **`f_act` "World activity" is US industrial production.** A deliberate,
   documented proxy (DECISIONS.md D1): the domain brief itself flags that the
   US is an imperfect proxy for Australia's Asia-weighted trading partners.
   The intended trade-weighted alternative (`f_act_tw`, OECD G20 GDP via
   DBnomics) is **stale on DBnomics (ends 2023Q3)** and is therefore parked in
   the unused `alt_foreign` set rather than the default suite. No fresh,
   key-free, machine-readable trade-weighted partner GDP series was found
   (the RBA does not publish its trading-partner GDP aggregate as a
   statistical table series).
2. **`cpi_inflation` index level is constructed.** The RBA publishes the
   trimmed mean as a quarterly % change; the data layer cumulates it to an
   index (`100·cumprod(1+q/100)`) purely so the uniform dlog transform
   applies. Round-trip error ≈ 8e-14 (exact to float precision). The model
   forecasts log changes, ~0.2bp below the published arithmetic % change at
   typical inflation rates — negligible and consistent across all dlog series.
3. **`f_act` freshness.** The key-free IFS series lags ~12–18 months
   (currently ends 2024Q4) and, because the panel is trimmed to balance,
   **the entire real-data panel ends 2024Q4** — the "current" forecast origin
   is ~6 quarters stale. With `FRED_API_KEY` set, `f_act` switches to FRED
   INDPRO and the panel extends to the ABS/RBA frontier (2026Q1). This is the
   single biggest practical limitation of the key-free real-data mode.

## 4. Synthetic data inventory

1. **`data.source: synthetic` (the default)** simulates the *entire* panel
   from a block-exogenous DGP. This exists so the pipeline runs offline
   (institutional requirement) and so the block-exogeneity diagnostic has a
   known ground truth. It is clearly labelled in every output (the report
   states the data source and warns that synthetic results validate
   machinery, not Australian dynamics). **No synthetic values are ever mixed
   into the real-data path** — if a real series cannot be obtained, the
   pipeline stops with an explicit error rather than substituting.
2. **Steady-state prior anchors** (`ss_mean`: inflation 2.5 % p.a. target
   midpoint, potential growth 2.8 %, NAIRU 4.5 %, neutral cash rate 3.5 %)
   are judgment numbers entering through priors, not data; sd's are
   configured and documented (DECISIONS.md D5).
3. **External forecasts hook** (`read_external_forecasts`) returns NULL
   unless the user supplies a real file — no placeholder data.

## 5. Simplifying assumptions (accepted and documented)

1. **Quarterly aggregation = within-quarter mean** for daily/monthly series
   (cash rate, bond yield, fed funds, commodity index, unemployment,
   employment, monthly IP). Quarterly-average rates are the standard VAR
   convention; end-of-quarter is the main alternative and would slightly
   change rate dynamics. Native-quarterly series are unaffected.
2. **Cash rate = the target, not the realized interbank rate.** Identical in
   practice over the modelled sample (1997Q3–), where the target is announced
   and the interbank rate tracks it within a basis point.
3. **Final-vintage data, not real-time vintages.** The pseudo-real-time
   evaluation uses today's published history at every origin. Revisions to
   GDP/employment are material in Australia; rankings sensitive to revisions
   (especially GDP at h=1–2) should be read with care. Stated loudly in the
   report; building a vintage database is the natural extension.
4. **Balanced-panel trim.** The common sample is the intersection across all
   active variables: it starts 1997Q3 (WPI's first observation) and ends at
   the stalest series (see §3.3). Alternative: drop WPI from the medium set
   and gain ~7 years of history — rejected for now because wages are central
   to the inflation block, but it is a one-line config change.
5. **`tot` residual Q1 pattern (p = 0.003).** The series is the ABS SA index;
   the pattern (Q1 mean +1.5pp vs ~0, sd 3.4, lag-1 autocorr 0.21) reflects
   bulk-commodity contract repricing historically clustered in Q1 — a price
   phenomenon, not a removable seasonal in the SA sense, and the simple
   F-test overstates significance under autocorrelation. Kept as is.
6. **The default `start: 1990` is aspirational for the real panel** — the
   actual start is determined by the balanced trim (1997Q3). Config start
   only truncates from the left.
7. **No seasonal adjustment is performed in-pipeline.** All series are sourced
   already-SA (or are concepts that need no SA: rates, SDR commodity prices,
   real TWI). If a future series is only available Original, adjust at source
   or add an SA step — do not dlog an Original series (the WPI bug in §2 is
   the cautionary tale).
