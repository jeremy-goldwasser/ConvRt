test_that("gi_discrete_gamma_delay returns a normalized pmf", {
  d <- gi_discrete_gamma_delay(mean = 3.5, sd = 1.8)
  expect_named(d, c("pmf", "shape", "rate"))
  expect_equal(sum(d$pmf), 1, tolerance = 1e-8)
  expect_true(all(d$pmf >= 0))
  # discretized mean should be close to the requested mean
  k <- seq_along(d$pmf)
  expect_equal(sum(k * d$pmf), 3.5, tolerance = 0.3)
})

test_that("gi_from_compartmental yields a valid generation-interval pmf", {
  g <- gi_from_compartmental(
    pi_lat = gi_discrete_gamma_delay(2.0, 1.2)$pmf,
    pi_IR  = gi_discrete_gamma_delay(2.75, 1.0)$pmf)$g
  expect_true(all(g >= 0))
  expect_equal(sum(g), 1, tolerance = 1e-6)
})
