test_that("frechet_median: mu is the median", {
  for (mu in c(0.5, 1, 4)) {
    # a finite median does not require a finite mean, so nu <= 1 is allowed
    for (nu in c(0.5, 1, 2, 5)) {
      expect_eps(
        stats::integrate(dfrechet_median, 0, mu, mu = mu, nu = nu)$value,
        0.5,
        eps = 1e-6
      )
    }
  }
  set.seed(1123)
  expect_eps(
    median(rfrechet_median(1e5, mu = 2, nu = 0.5)),
    2,
    eps = 0.03,
    relative = TRUE
  )
})

test_that("rng_lookup('frechet') matches brms's mean parameterization", {
  set.seed(1124)
  draws <- rng_lookup("frechet")(1e5, mu = 2, nu = 5)
  expect_eps(mean(draws), 2, eps = 0.02, relative = TRUE)
})
