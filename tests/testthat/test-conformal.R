test_that("fit_convrt_realtime_conformal widens the trailing edge band", {
  skip_on_cran()  # refits the model dozens of times
  s <- sim_outbreak()
  # Truncate to a healthy-count vintage (day 108) so the edge is well-identified.
  W <- 108L
  res <- fit_convrt_realtime_conformal(
    obs_inc = s$Y_obs[1:W], dates = s$dates[1:W], g = s$g,
    mean_EY = s$mean_EY, sd_EY = s$sd_EY, severity = s$rho,
    first_rt_date         = s$dates[21],
    likelihood_start_date = s$dates[35],
    n_calib = 40L)

  expect_true(res$n_vintages > 20L)
  expect_true(all(c("Rt_lo_wald", "Rt_hi_wald", "ci_source", "conformal_q")
                  %in% names(res$rt_df)))

  # At least some trailing rows should have been conformalized.
  expect_true(any(res$rt_df$ci_source == "conformal"))

  # On the conformalized rows the band should be at least as wide as the Wald band
  # (conformal corrects real-time edge under-coverage by widening).
  cf <- res$rt_df[res$rt_df$ci_source == "conformal", ]
  w_conf <- cf$Rt_hi      - cf$Rt_lo
  w_wald <- cf$Rt_hi_wald - cf$Rt_lo_wald
  expect_true(mean(w_conf) >= mean(w_wald))
})
