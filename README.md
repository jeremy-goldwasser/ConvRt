# ConvRt — Generation-Interval R<sub>t</sub>

ConvRt estimates the time-varying reproduction number R<sub>t</sub> with a
regularized identity-link Poisson model: it fits a smoothing-spline R<sub>t</sub>
through the **renewal equation** and predicts the **observed counts**. The only
epidemiological inputs are a **generation interval** `g` and a **reporting delay**
`pi_EY` — each a single PMF you can take from the literature or an external
estimate.

It provides retrospective (end-of-season) and real-time (right-censored) fits,
with pointwise, simultaneous, and split-conformal confidence intervals, plus
weekly-aggregated and trend-filtered variants.

---

## Install

```r
# install.packages("remotes")
remotes::install_github("jeremy-goldwasser/ConvRt")
```

This pulls ConvRt's required dependencies (dplyr, tibble, tidyr; splines/stats/utils
ship with R). To also install what the demos and optional features use, run
[`install_packages.R`](install_packages.R). Two dependencies are **optional** and
not installed by default:

- **mgcv** — only for `gi_deconvolve_exposures(link = "log")`.
- **CVXR** — only for the trend-filter variant (`fit_convrt_tf_*`). It is a heavy,
  rarely-needed solver stack; the TF functions print an install hint if you call
  them without it. Everything else works without CVXR.

---

## Quick start

```r
library(ConvRt)

# 1. Epidemiological inputs: two discretized-Gamma PMFs (use any literature PMFs)
g <- gi_discrete_gamma_delay(mean = 3.5, sd = 1.8)$pmf   # generation interval
# reporting delay is passed as its mean/sd (mean_EY, sd_EY) below

# 2. Retrospective (end-of-season) fit
retro <- fit_convrt_retrospective(
  obs_inc = y, dates = dates, g = g,
  mean_EY = 5.7, sd_EY = 2.3, severity = 0.015,
  first_rt_date         = as.Date("2022-07-01"),   # start estimating R_t here (burn-in anchor)
  likelihood_start_date = as.Date("2022-07-22"))   # first day entering the likelihood
retro$rt_df    # tibble: date, day, Rt_mean, Rt_lo, Rt_hi

# 3. Real-time (right-censored vintage) fit
rt <- fit_convrt_realtime(
  obs_inc = y_vintage, dates = dates_vintage, g = g,
  mean_EY = 5.7, sd_EY = 2.3, severity = 0.015,
  first_rt_date         = as.Date("2022-07-01"),
  likelihood_start_date = as.Date("2022-07-22"))

# 4. Quick plot of any fit (needs ggplot2)
plot_convrt(retro)

# 5. Real-time nowcast with split-conformal edge CIs (refits at many past
#    vintages internally, so it is slower):
nc <- fit_convrt_realtime_conformal(
  obs_inc = y_vintage, dates = dates_vintage, g = g,
  mean_EY = 5.7, sd_EY = 2.3, severity = 0.015,
  first_rt_date         = as.Date("2022-07-01"),
  likelihood_start_date = as.Date("2022-07-22"),
  n_calib = 45)
nc$rt_df    # trailing-edge Rt_lo/Rt_hi are conformal; ci_source marks each row
```

Pass `verbose = TRUE` to either fit wrapper to print CV/tuning progress (silent by
default).

---

## Model

```
Y_t ~ Poisson( mu_t ),   mu_t = rho_t * omega_(t mod 7) * sum_{s<t} R_s * X_{t,s}
X_{t,s} = Lambda_s * pi_EY[t-s],     Lambda_s = sum_{k>=1} g_k * X_exposure[s-k]
```

- `X_exposure` — exposure (infection) incidence, recovered by deconvolving the
  observations against the **reporting delay** `pi_EY` only.
- `Lambda_s` — the **renewal force**: trailing exposures convolved with the
  generation interval `g`. The renewal equation is `E[X_s] = R_s * Lambda_s`.
- `pi_EY` — exposure→observed-outcome (e.g. exposure→hospitalization) delay,
  parameterized by `mean_EY`/`sd_EY` (discretized Gamma).
- `rho_t` — **severity** (reporting/ascertainment rate); scalar or per-day vector.
- `omega` — optional multiplicative **day-of-week** reporting effect (product 1).
- `R_t = (B theta)_t` — a natural cubic smoothing spline.

---

## Methodological choices & defaults

Several steps involve a defensible-but-not-unique modeling choice. They are
collected here, with the rationale and how to change each. **All defaults below
are the values used by the `fit_convrt_*` wrappers.**

### 1. Recovering exposures (deconvolution)

| Choice | Default | Notes |
| --- | --- | --- |
| Link | `identity` | Matches the data-generating model (counts are a *linear* convolution of exposures × severity). `link = "log"` is available (needs **mgcv**) for log-spline deconvolution. |
| Penalty tuning | **GCV** | The exposure curve is a *nuisance* pre-step, so its smoothing penalty is chosen by fast generalized cross-validation rather than K-fold CV (which is reserved for the penalty that actually matters — the R<sub>t</sub> spline). |
| Burn-in | 30 days | Early days have too little trailing data to deconvolve; excluded from the deconvolution likelihood. |

### 2. R<sub>t</sub> spline basis & roughness penalty

| Choice | Default | Notes |
| --- | --- | --- |
| Basis | natural cubic spline (linear tails) | `R_t = B theta`. Linear tails avoid wild extrapolation at the edges. |
| Knot spacing | `knot_step = 5` days | Denser knots → more flexible R<sub>t</sub> but more reliance on the penalty. |
| Penalty | integrated squared 2nd derivative | The standard smoothing-spline roughness penalty; strength `lambda` is tuned (below). |

### 3. Choosing the smoothing strength `lambda`

| Choice | Default | Notes |
| --- | --- | --- |
| Method | **K-fold cross-validation** | On the observation-level Poisson likelihood. |
| Folds (`nfold`) | 5 | Folds are **deterministic** (interleaved: fold *k* = rows *k, k+5, k+10, …*), so λ-selection is **reproducible without setting a seed**. |
| Fold loss (`error_measure`) | `deviance` | Poisson deviance. Alternatives: `mse`, `mae`. |
| Selection rule (`cv_select_rule`) | `min` | Smallest-CV-error λ. `1se` (largest λ within 1 SE of the min — smoother, more conservative) is also available. *(Note: the low-level `gi_select_lambda_cv` defaults to `1se`, but the wrappers override to `min`.)* |
| λ grid (`lam_grid`) | `10^seq(-2, 8, length.out = 30)` | Widen if the chosen λ lands at a grid endpoint. |

### 4. Real-time edge taper (`gamma`)

The real-time fit adds a CDF-tapered penalty on the right edge (where data is
still arriving) to stabilize the nowcast.

| Choice | Default | Notes |
| --- | --- | --- |
| Strength selection | **forward validation (FV)** | One-step-ahead prediction error: refit through day *s*, predict `Y_{s+1}`, for the last `n_fv` days. This is the honest real-time analogue of CV. |
| FV folds (`n_fv`) | 7 (daily) / 4 (weekly) | Number of trailing one-step-ahead holdouts. |
| Selection rule | min FV error | Picks the `gamma` with lowest mean one-step-ahead error. |
| `gamma` grid | `10^seq(-4, 5, length.out = 18)` | |

### 5. Confidence intervals

| Method | Where | Notes |
| --- | --- | --- |
| Pointwise (Laplace/Wald) | `fit_convrt_*$rt_df` (`Rt_lo`/`Rt_hi`) | Sandwich covariance `M_pen^{-1} Z'WZ M_pen^{-1}`. `level = 0.95`. |
| Quasi-Poisson overdispersion | `overdispersion = TRUE` (default) | Inflates SEs by `sqrt(phi)`, with `phi` the Pearson statistic over effective df. Point estimates are unchanged. Set `FALSE` for pure Poisson SEs. |
| Simultaneous band | `gi_extract_rt_simband()` | Sup-norm band drawn from the Laplace covariance (`n_sim = 5000`). |
| Real-time **split-conformal** | `fit_convrt_realtime_conformal()` (one call) or `gi_apply_conformal_to_rt_df()` (manual) | Calibrates the trailing-edge band from how past vintages' real-time estimates differed from the later settled values. `ci_level = 0.90`, `buffer_days = 14` (a date is "settled" once this far behind an edge), `max_lag = 13` (corrects the trailing 14 days), `n_calib = 45` prior vintages. Refits the model many times — slow. |

### 6. Day-of-week effects

Pass `dow_dates` to model multiplicative weekday reporting effects (product
constrained to 1). For DoW-affected daily data you can tune λ with the
DoW-aware selector `gi_select_lambda_cv_dow()` directly.

---

## Demos

Both demos are fully self-contained: they simulate a two-wave outbreak from the
model (no external data) and recover R<sub>t</sub>.

- [`examples/demo.R`](examples/demo.R) — retrospective + a single real-time fit
  with pointwise CIs. Fast (seconds). Writes `examples/demo_rt.pdf`.
- [`examples/demo_realtime_ci.R`](examples/demo_realtime_ci.R) — real-time
  **split-conformal** CIs. Refits the model at ~50 past vintages to calibrate the
  trailing-edge band, so it is much slower (~20 s). Writes
  `examples/demo_realtime_ci.pdf`.

```sh
Rscript examples/demo.R
Rscript examples/demo_realtime_ci.R
```

---

## Practical guidance

- **Low counts ruin the edges.** Below ~5 observed cases/day the convolution
  likelihood has too little signal to resolve R<sub>t</sub>; estimates near such
  troughs (and in the first weeks, before the reporting delay has filled) are
  unreliable. Start `first_rt_date`/`likelihood_start_date` a few weeks into the
  series, and don't read too much into the extreme tails. The real-time CI demo
  deliberately places its target vintage in a healthy-count region for exactly
  this reason.
- **Generation interval must be a PMF.** `g` should sum to 1; `build_gi_design`
  warns otherwise.
- **Reproducibility.** λ-selection is deterministic (no seed needed). The
  trend-filter CV (`fit_convrt_tf_*`) uses seeded random folds (`seed = 1`).

---

## API reference

**Inputs / delays**
- `gi_discrete_gamma_delay(mean, sd)` — discretized Gamma PMF.
- `gi_from_compartmental(pi_lat, pi_IR, ...)` — *optional* helper: the generation
  interval implied by SEIR-style latent + infectious-period PMFs. (Not needed if
  you supply `g` directly.)

**Main fits**
- `fit_convrt_retrospective(...)` — end-of-season R<sub>t</sub> + pointwise CIs.
- `fit_convrt_realtime(...)` — right-censored R<sub>t</sub> + tapered tail + CIs.
- `fit_convrt_weekly_retrospective(...)`, `fit_convrt_weekly_realtime(...)` —
  weekly-aggregated variants.
- `fit_convrt_tf_retrospective(...)` — trend-filter variant (needs **CVXR**).

**Confidence intervals**
- `fit_convrt_realtime_conformal(...)` — one-call real-time nowcast with
  split-conformal edge CIs (refits at `n_calib` past vintages internally).
- `gi_extract_rt_simband(...)` — simultaneous (sup-norm) band.
- `gi_apply_conformal_to_rt_df(...)`, `gi_compute_conformal_q(...)`,
  `gi_load_daily_vintage_cache(...)` — lower-level conformal building blocks.

**Plotting**
- `plot_convrt(fit, truth = NULL)` — ggplot of the R<sub>t</sub> curve + band.

**Building blocks** (call directly for custom pipelines): `gi_deconvolve_exposures`,
`gi_renewal_force`, `build_gi_design`, `gi_solve`, `gi_select_lambda_cv`
(+`_dow`), `gi_extract_rt`, `gi_estimate_dispersion`, `gi_build_tapered_penalty`,
`gi_tune_gamma_fv`, `gi_apply_severity_to_design`, `gi_enforce_likelihood_start`,
`gi_rt_df_from_fit`.

---

## Layout

```
DESCRIPTION / NAMESPACE   package metadata
R/
  delays.R          delay PMFs (gi_discrete_gamma_delay, gi_from_compartmental)
  deconvolution.R   gi_deconvolve_exposures (reporting delay only)
  renewal.R         gi_renewal_force
  design.R          build_gi_design (renewal design)
  wrappers.R        fit_convrt_retrospective / fit_convrt_realtime
  solve.R           gi_solve (penalized IRLS + KKT)
  lambda_select.R   gi_select_lambda_cv (+ _dow)
  extract.R         gi_extract_rt, gi_estimate_dispersion, simband
  taper.R           gi_build_tapered_penalty, gi_tune_gamma_fv
  realtime.R        severity / likelihood-start / rt_df helpers
  conformal.R       split-conformal real-time CIs
  weekly.R          weekly aggregation + weekly wrappers
  tf.R              trend-filter variant (optional; needs CVXR)
examples/
  demo.R                 retrospective + real-time, pointwise CIs
  demo_realtime_ci.R     real-time split-conformal CIs
```
