test_that("plot_convrt returns a ggplot for a fit and for a bare rt_df", {
  skip_if_not_installed("ggplot2")
  s <- sim_outbreak()
  fit <- fit_convrt_retrospective(
    obs_inc = s$Y_obs, dates = s$dates, g = s$g,
    mean_EY = s$mean_EY, sd_EY = s$sd_EY, severity = s$rho,
    first_rt_date = s$dates[21], likelihood_start_date = s$dates[35])

  expect_s3_class(plot_convrt(fit), "ggplot")
  expect_s3_class(plot_convrt(fit$rt_df), "ggplot")
  expect_s3_class(
    plot_convrt(fit, truth = data.frame(date = s$dates, Rt_true = s$Rt_show)),
    "ggplot")
})

test_that("plot_convrt rejects input without the required columns", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_convrt(data.frame(a = 1, b = 2)), "rt_df")
})
