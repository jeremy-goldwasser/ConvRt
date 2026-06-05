# Shared simulator: the same two-wave outbreak DGP used in examples/demo.R,
# returned as a tidy list so tests can fit and check recovery.
sim_outbreak <- function(seed = 2L, n_show = 170L, burn_in = 120L, rho = 0.035) {
  set.seed(seed)
  n_total <- n_show + burn_in
  tt_all  <- seq.int(-burn_in + 1L, n_show)
  rt_fn <- function(t) {
    1.0 +
    0.50 * exp(-((t -  55) / 25)^2) -
    0.28 * exp(-((t - 120) / 22)^2) +
    0.30 * exp(-((t - 155) / 18)^2)
  }
  Rt_true <- rt_fn(tt_all)
  g       <- gi_discrete_gamma_delay(mean = 3.5, sd = 1.8)$pmf
  mean_EY <- 5.7; sd_EY <- 2.3
  pi_EY   <- gi_discrete_gamma_delay(mean_EY, sd_EY)$pmf

  X <- numeric(n_total); X[seq_len(min(5L, length(g)))] <- 30
  for (t in (min(5L, length(g)) + 1L):n_total) {
    k <- 1:min(length(g), t - 1L)
    X[t] <- rpois(1L, Rt_true[t] * sum(g[k] * X[t - k]))
  }
  muY <- vapply(seq_len(n_total), function(t) {
    k <- 1:min(length(pi_EY), t - 1L); if (t == 1L) 0 else rho * sum(pi_EY[k] * X[t - k])
  }, numeric(1))
  Y <- rpois(n_total, pmax(muY, 0))
  keep <- (burn_in + 1L):n_total

  list(dates   = as.Date("2023-01-01") + 0:(n_show - 1L),
       Y_obs   = Y[keep],
       Rt_show = Rt_true[keep],
       g = g, mean_EY = mean_EY, sd_EY = sd_EY, rho = rho)
}
