test_that("weibull_median: mu is the median", {
  for (mu in c(0.5, 2, 7.3)) {
    for (k in c(0.5, 1, 3)) {
      expect_eps(
        stats::integrate(dweibull_median, 0, mu, mu = mu, k = k)$value,
        0.5,
        eps = 1e-6
      )
    }
  }
  set.seed(2323)
  expect_eps(
    median(rweibull_median(1e5, mu = 2, k = 1)),
    2,
    eps = 0.02,
    relative = TRUE
  )
})

test_that("rng_lookup('weibull') matches brms's mean parameterization", {
  set.seed(2324)
  draws <- rng_lookup("weibull")(1e5, mu = 2, shape = 1.5)
  expect_eps(mean(draws), 2, eps = 0.02, relative = TRUE)
})
