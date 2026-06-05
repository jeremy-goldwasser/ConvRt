test_that("retrospective fit recovers R_t on simulated data", {
  s <- sim_outbreak()
  fit <- fit_convrt_retrospective(
    obs_inc = s$Y_obs, dates = s$dates, g = s$g,
    mean_EY = s$mean_EY, sd_EY = s$sd_EY, severity = s$rho,
    first_rt_date         = s$dates[21],
    likelihood_start_date = s$dates[35],
    knot_step             = 5L)

  expect_named(fit$rt_df, c("date", "day", "Rt_mean", "Rt_lo", "Rt_hi"))
  expect_true(all(is.finite(fit$rt_df$Rt_mean)))
  expect_true(all(fit$rt_df$Rt_lo <= fit$rt_df$Rt_hi))

  # Accuracy + coverage over the well-identified interior.
  ev <- merge(fit$rt_df, data.frame(date = s$dates, Rt_true = s$Rt_show), by = "date")
  ev <- ev[ev$date >= s$dates[35] & ev$date <= s$dates[length(s$dates) - 3L], ]
  mae <- mean(abs(ev$Rt_mean - ev$Rt_true))
  cov <- mean(ev$Rt_true >= ev$Rt_lo & ev$Rt_true <= ev$Rt_hi)
  expect_lt(mae, 0.10)
  expect_gt(cov, 0.80)
})

test_that("fit is reproducible without setting a seed (deterministic CV folds)", {
  s <- sim_outbreak()
  args <- list(obs_inc = s$Y_obs, dates = s$dates, g = s$g,
               mean_EY = s$mean_EY, sd_EY = s$sd_EY, severity = s$rho,
               first_rt_date = s$dates[21], likelihood_start_date = s$dates[35])
  a <- do.call(fit_convrt_retrospective, args)
  b <- do.call(fit_convrt_retrospective, args)
  expect_equal(a$lam, b$lam)
  expect_equal(a$rt_df$Rt_mean, b$rt_df$Rt_mean)
})

test_that("real-time fit returns a sane edge estimate", {
  s <- sim_outbreak()
  W  <- 90L
  rt <- fit_convrt_realtime(
    obs_inc = s$Y_obs[1:W], dates = s$dates[1:W], g = s$g,
    mean_EY = s$mean_EY, sd_EY = s$sd_EY, severity = s$rho,
    first_rt_date         = s$dates[21],
    likelihood_start_date = s$dates[35])
  expect_true(rt$gamma >= 0)
  edge <- tail(rt$rt_df$Rt_mean, 1L)
  expect_true(is.finite(edge) && edge > 0 && edge < 5)
})

test_that("fits run silently by default and talk when verbose = TRUE", {
  s <- sim_outbreak()
  expect_silent(
    fit_convrt_retrospective(
      obs_inc = s$Y_obs, dates = s$dates, g = s$g,
      mean_EY = s$mean_EY, sd_EY = s$sd_EY, severity = s$rho,
      first_rt_date = s$dates[21], likelihood_start_date = s$dates[35]))
  expect_output(
    fit_convrt_retrospective(
      obs_inc = s$Y_obs, dates = s$dates, g = s$g,
      mean_EY = s$mean_EY, sd_EY = s$sd_EY, severity = s$rho,
      first_rt_date = s$dates[21], likelihood_start_date = s$dates[35],
      verbose = TRUE),
    "Best lambda")
})
