# examples/demo_realtime_ci.R
#
# Real-time *conformal* confidence intervals for ConvRt.
#
# The pointwise (Wald) CI from a single real-time fit understates the true
# uncertainty at the right edge: the last few days are still being revised as
# data arrives.  Split-conformal calibration fixes this by measuring, on past
# vintages, how far the real-time edge estimate landed from the (later, stable)
# retrospective value -- and inflating the edge band to match.
#
# This requires REFITTING the model at many past vintages, so it is much slower
# than examples/demo.R (here: ~50 real-time fits, a handful of seconds each).
#
# Self-contained: same simulated two-wave outbreak as demo.R, no external data.
# Depends only on the installed ConvRt package.

suppressPackageStartupMessages({
  library(ConvRt); library(ggplot2); library(dplyr); library(tibble)
})

.this <- (function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f) else NULL
})()
convrt_root <- if (!is.null(.this)) dirname(dirname(.this)) else getwd()

set.seed(2)

# ---- Truth + simulation (identical DGP to demo.R) ---------------------------
n_show   <- 170L
burn_in  <- 120L
n_total  <- n_show + burn_in
tt_all   <- seq.int(-burn_in + 1L, n_show)
rt_fn <- function(t) {
  1.0 +
  0.50 * exp(-((t -  55) / 25)^2) -
  0.28 * exp(-((t - 120) / 22)^2) +
  0.30 * exp(-((t - 155) / 18)^2)
}
Rt_true <- rt_fn(tt_all)

g       <- gi_discrete_gamma_delay(mean = 3.5, sd = 1.8)$pmf
mean_EY <- 5.7; sd_EY <- 2.3; pi_EY <- gi_discrete_gamma_delay(mean_EY, sd_EY)$pmf
rho     <- 0.035

X <- numeric(n_total); X[seq_len(min(5L, length(g)))] <- 30
for (t in (min(5L, length(g)) + 1L):n_total) {
  k   <- 1:min(length(g), t - 1L)
  X[t] <- rpois(1L, Rt_true[t] * sum(g[k] * X[t - k]))
}
muY <- vapply(seq_len(n_total), function(t) {
  k <- 1:min(length(pi_EY), t - 1L); if (t == 1L) 0 else rho * sum(pi_EY[k] * X[t - k])
}, numeric(1))
Y <- rpois(n_total, pmax(muY, 0))

keep    <- (burn_in + 1L):n_total
dates   <- as.Date("2023-01-01") + 0:(n_show - 1L)
Y_obs   <- Y[keep]
Rt_show <- Rt_true[keep]
truth   <- tibble(date = dates, Rt_true = Rt_show)

# ---- Vintage schedule -------------------------------------------------------
# Target vintage (the "now" we report) is placed at day 108 -- past the first
# peak but where counts are still healthy (~180/day).  We deliberately do NOT
# push the edge into the post-peak trough (days >120, counts < 100): there the
# real-time signal is too thin and every method's edge estimate degrades, which
# would conflate conformal calibration with genuine non-identifiability.
W_target  <- 108L                  # target vintage = "today"
W_start   <- 60L                   # earliest calibration vintage (counts already healthy)
fit_days  <- 21L                   # first_rt_date offset (3 weeks of burn-in)
lik_days  <- 35L                   # likelihood_start offset (let reporting delay fill)
ci_level    <- 0.90                # nominal coverage of the produced interval
buffer_days <- 14L                 # a date is "settled" once it is this far behind a vintage edge
max_lag     <- 13L                 # conformal corrects the trailing 14 days (d_to_edge 0..13)

vintages <- W_start:W_target
cat(sprintf("Refitting ConvRt at %d daily vintages (day %d .. %d). This is the slow part...\n",
            length(vintages), W_start, W_target))

# One real-time fit at vintage W -> its rt_df, tagged with the vintage date.
# (fit_convrt_realtime is silent unless verbose = TRUE.)
fit_vintage <- function(W) {
  rt <- tryCatch(
    fit_convrt_realtime(
      obs_inc = Y_obs[1:W], dates = dates[1:W], g = g,
      mean_EY = mean_EY, sd_EY = sd_EY, severity = rho,
      first_rt_date         = dates[fit_days],
      likelihood_start_date = dates[lik_days],
      knot_step             = 5L),
    error = function(e) { message(sprintf("  vintage %d failed: %s", W, conditionMessage(e))); NULL })
  if (is.null(rt)) return(NULL)
  rt$rt_df$vintage <- dates[W]
  rt$rt_df
}

t0 <- Sys.time()
vintage_dfs <- lapply(vintages, function(W) {
  df <- fit_vintage(W)
  if (W %% 10L == 0L) cat(sprintf("  ...vintage %d/%d\n", W, W_target))
  df
})
daily_vintages <- bind_rows(Filter(Negate(is.null), vintage_dfs))
cat(sprintf("done: %d vintage fits in %.0fs\n",
            length(unique(daily_vintages$vintage)), as.numeric(Sys.time() - t0, units = "secs")))

# ---- Target fit + conformal calibration -------------------------------------
target_W   <- dates[W_target]
target_rt  <- daily_vintages |> filter(vintage == target_W)

conf <- gi_apply_conformal_to_rt_df(
  rt_df          = target_rt,
  daily_vintages = daily_vintages,
  target_W       = target_W,
  ci_level       = ci_level,
  buffer_days    = buffer_days,
  max_lag        = max_lag)

# ---- Metrics: coverage of the trailing edge, conformal vs Wald --------------
edge <- conf |>
  inner_join(truth, by = "date") |>
  filter(date >= target_W - max_lag, date <= target_W)
cov_conf <- mean(edge$Rt_true >= edge$Rt_lo      & edge$Rt_true <= edge$Rt_hi)
cov_wald <- mean(edge$Rt_true >= edge$Rt_lo_wald & edge$Rt_true <= edge$Rt_hi_wald)
w_conf   <- mean(edge$Rt_hi - edge$Rt_lo)
w_wald   <- mean(edge$Rt_hi_wald - edge$Rt_lo_wald)
n_cal    <- max(edge$conformal_n_cal, na.rm = TRUE)
cat(sprintf("\nTrailing %d days @ %d%% nominal:\n", max_lag + 1L, round(100 * ci_level)))
cat(sprintf("  Wald      : coverage %.2f  mean width %.3f\n", cov_wald, w_wald))
cat(sprintf("  conformal : coverage %.2f  mean width %.3f  (n_cal up to %d)\n",
            cov_conf, w_conf, n_cal))

# ---- Plot -------------------------------------------------------------------
plt <- conf |> filter(date >= target_W - 45)
truth_plt <- truth |> filter(date >= target_W - 45, date <= target_W)

p <- ggplot(plt, aes(date)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey55") +
  geom_ribbon(aes(ymin = Rt_lo_wald, ymax = Rt_hi_wald, fill = "Wald"), alpha = 0.30) +
  geom_ribbon(aes(ymin = Rt_lo,      ymax = Rt_hi,      fill = "conformal"), alpha = 0.30) +
  geom_line(data = truth_plt, aes(date, Rt_true), colour = "grey10", linewidth = 1.2) +
  geom_line(aes(y = Rt_mean), colour = "#D55E00", linewidth = 1.1) +
  geom_point(data = tail(plt, 1L), aes(y = Rt_mean), colour = "#D55E00", size = 2.6) +
  geom_vline(xintercept = target_W, linetype = "dashed", colour = "#D55E00", alpha = 0.5) +
  scale_fill_manual(values = c(Wald = "#999999", conformal = "#0072B2"), name = "edge CI") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(title = "ConvRt real-time conformal CIs",
       subtitle = sprintf(
         "black = truth, orange = real-time nowcast (vintage %s); conformal band widens the trailing %d days",
         format(target_W), max_lag + 1L),
       x = NULL, y = expression(R[t])) +
  theme_bw(base_size = 14) +
  theme(plot.subtitle = element_text(colour = "grey30"), legend.position = "top")

out_pdf <- file.path(convrt_root, "examples", "demo_realtime_ci.pdf")
ggsave(out_pdf, p, width = 10, height = 6)
cat(sprintf("wrote %s\n", out_pdf))
