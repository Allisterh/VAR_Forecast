# DATA_AUDIT.md — full audit of data sources, updated 2026-06-13

Method: every configured series was pulled from its provider and its **metadata**
(series description, table, series type, unit, frequency, coverage) checked
against the intended economic concept; every dlog-transformed series was
additionally tested for residual quarterly seasonality (F-test of quarter
dummies on quarterly log changes); constructed series were round-trip tested.
Sections: verified series → corrections made → proxies and constructed data →
synthetic data inventory → simplifying assumptions.

**Production configuration (this revision).** `data.source: real` is the
default. The foreign block is sourced from **FRED** (key required, in a
gitignored `.Renviron`): `f_act` = US **real GDP** (GDPC1), `f_rate` = US
**effective fed funds** (FEDFUNDS). Australian data comes key-free from RBA
(`readrba`) and ABS (`readabs`); commodity prices from RBA. The balanced panel
runs **1997Q4–2026Q1** (114 quarters, all 13 series real). The synthetic
generator is retained as an explicit opt-in (`data.source: synthetic`) for
tests/CI and the block-exogeneity ground-truth diagnostic only.

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
| `f_comm` | RBA GRCPAISDR (I2) | "Index of commodity prices; All items; **SDR**", monthly, 1982– (SDR strips the endogenous AUD — correct for the price-taker foreign block) |
| `f_rate` | **FRED FEDFUNDS** | "Effective Federal Funds Rate", %, monthly, NSA (no seasonal concept for a rate); used as level, delta=1 |
| `f_act` | **FRED GDPC1** | "Real Gross Domestic Product", **chain-volume**, **Seasonally Adjusted**, Bil. Chained 2017 $, **quarterly-native** (1947–); used as dlog (qtr % growth), delta=0 — the same real-GDP concept as the domestic `gdp_growth`, so the foreign/domestic activity comovement is consistent |

Empirical seasonality F-tests on quarterly log changes (real data): all clean
(p > 0.05) except `tot` (p = 0.003; see §5). `wpi_growth` was failing
(p ≈ 1e-10) before the §2 correction. `f_act` = GDPC1 is published already
seasonally adjusted and quarterly-native, so the seasonality test is trivially
clean and within-quarter aggregation is a no-op for it.

## 2. Corrections made

1. **Foreign block moved to FRED (this revision).** `f_act` changed from
   DBnomics IMF/IFS US **industrial production** to **FRED GDPC1 (US real
   GDP)** — a better world-activity proxy for a quarterly macro VAR: it is
   quarterly-native (no frequency conversion) and the *same concept* (real
   GDP) as the domestic activity variable, rather than the narrower
   industrial-production index. `f_rate` changed from DBnomics
   `FED/H15/RIFSPFF_N.M` to **FRED FEDFUNDS** (same effective-rate concept,
   primary not fallback, fresher). This removes the previous "panel ends
   2024Q4" staleness (the IFS series lagged ~a year): the panel now reaches
   2026Q1. DBnomics remains only for the parked `alt_foreign` trade-weighted
   variant.
2. **WPI was the wrong series type.** The config used `A2603609J`, the
   **Original (non-seasonally-adjusted)** quarterly WPI index — confirmed by
   ABS metadata (`series_type: Original`) and a glaring seasonal in its log
   changes (F-test p ≈ 1e-10). Replaced with `A2713849C` (identical concept —
   total hourly rates of pay ex bonuses, Australia, private & public, all
   industries — **Seasonally Adjusted**).
3. **Commodity and terms-of-trade IDs.** Commodity index `GRCPAISAD` (A$/bulk
   -spot, starts 2009, embeds the endogenous AUD) → **`GRCPAISDR`**
   (SDR-denominated world price, 1982–). ToT `A2303731T` (did not resolve) →
   verified SA index `A2304200A`.

## 3. Proxies and constructed series (real-data mode)

These are **real published data**, but not literally the named concept:

1. **`f_act` "World activity" is US real GDP (a US-only proxy).** GDPC1 is the
   right *concept* (real GDP, the same as the domestic activity variable) but
   is US-only, a deliberate, documented proxy (DECISIONS.md D1): the domain
   brief flags that the US is an imperfect proxy for Australia's Asia-weighted
   trading partners. The intended trade-weighted alternative (`f_act_tw`, OECD
   G20 GDP via DBnomics) is **stale on DBnomics (ends 2023Q3)** and is parked
   in the unused `alt_foreign` set, not the default suite — no fresh, machine
   -readable trade-weighted partner-GDP series was found (the RBA does not
   publish its trading-partner GDP aggregate as a statistical-table series).
   This US-only-world proxy is the principal remaining modelling simplification
   and is acceptable for production but flagged.
2. **`cpi_inflation` index level is constructed.** The RBA publishes the
   trimmed mean as a quarterly % change; the data layer cumulates it to an
   index (`100·cumprod(1+q/100)`) purely so the uniform dlog transform
   applies. Round-trip error ≈ 8e-14 (exact to float precision). The model
   forecasts log changes, ~0.2bp below the published arithmetic % change at
   typical inflation rates — negligible and consistent across all dlog series.
3. **Frontier and publication lags.** With the foreign block on FRED (GDPC1
   published ~1 month after the quarter; FEDFUNDS ~1 month) and ABS/RBA at
   their normal cadence, the balanced panel reaches **2026Q1** — the binding
   (stalest) series are the quarterly-native ones (ABS national accounts, the
   RBA trimmed mean, FRERTWI, GDPC1), all of which have a 2026Q1 print. The
   `end: 2026-03-01` config pin and the coverage-aware quarterly aggregation
   (§5.7) together ensure no partial (sub-quarter) frontier observation enters
   the panel. The old key-free "ends 2024Q4 / a year behind" limitation no
   longer applies.

## 4. Synthetic data inventory

1. **`data.source: synthetic` (an explicit opt-in, no longer the default)**
   simulates the *entire* panel from a block-exogenous DGP. It is retained so
   the pipeline can run fully offline (tests/CI) and so the block-exogeneity
   diagnostic has a known ground truth — it is still exercised by the test
   suite (test-pipeline.R, test-covid.R). It is clearly labelled in every
   output (the report states the data source and warns that synthetic results
   validate machinery, not Australian dynamics). **No synthetic values are
   ever mixed into the real-data path** — `get_raw_data()` branches strictly
   on `data.source`, and if a real series cannot be obtained the pipeline
   stops with an explicit error rather than substituting.
2. **Steady-state prior anchors** (`ss_mean`: inflation 2.5 % p.a. target
   midpoint, potential growth 2.8 %, NAIRU 4.5 %, neutral cash rate 3.5 %)
   are judgment numbers entering through priors, not data; sd's are
   configured and documented (DECISIONS.md D5).
3. **External forecasts hook** (`read_external_forecasts`) returns NULL
   unless the user supplies a real file — no placeholder data.

## 5. Simplifying assumptions (accepted and documented)

1. **Quarterly aggregation = within-quarter mean** for daily/monthly series
   (cash rate, bond yield, fed funds, commodity index, unemployment,
   employment). Quarterly-average rates are the standard VAR convention;
   end-of-quarter is the main alternative and would slightly change rate
   dynamics. Native-quarterly series (US GDP/GDPC1, AU GDP/consumption/ToT,
   trimmed-mean CPI, WPI, real TWI) are unaffected.
2. **Cash rate = the target, not the realized interbank rate.** Identical in
   practice over the modelled sample (1997Q4–), where the target is announced
   and the interbank rate tracks it within a basis point.
3. **Final-vintage data, not real-time vintages.** The pseudo-real-time
   evaluation uses today's published history at every origin. Revisions to
   GDP/employment are material in Australia; rankings sensitive to revisions
   (especially GDP at h=1–2) should be read with care. Stated loudly in the
   report; building a vintage database is the natural extension.
4. **Balanced-panel trim.** The common sample is the intersection across all
   active variables: it starts **1997Q4** (the first dlog of WPI, whose level
   starts 1997Q3 — the binding series; GDPC1 from 1947 and FEDFUNDS from 1954
   do not bind) and ends at the stalest series (2026Q1, see §3.3). Alternative:
   drop WPI from the medium set and gain ~7 years of history — rejected for now
   because wages are central to the inflation block, but it is a one-line
   config change.
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
8. **Coverage-aware quarterly aggregation (this revision).** `to_quarterly()`
   drops (sets NA) any quarter with materially fewer source observations than
   a full quarter (< 80% of the fullest quarter's count), so a partial
   frontier quarter becomes NA and is trimmed rather than entering as a
   1–2-month average. Combined with the loud interior-NA error in
   `transform_data()` and the pct_change interior-NA guard, the real path
   never silently substitutes or shortens: a gap either trims cleanly
   (trailing) or errors (interior). The cache is keyed by provider+series_id
   so a changed series ID re-downloads rather than serving stale data.
