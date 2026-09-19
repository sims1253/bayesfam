test_that("frechet_median density matches brms::dfrechet with median scale", {
  mu <- 2
  nu <- 3
  x <- seq(from = 0.1, to = 40, length.out = 100)
  expect_eps(
    dfrechet_median(x, mu = mu, nu = nu),
    brms::dfrechet(x, loc = 0, scale = mu * log(2)^(1 / nu), shape = nu),
    eps = 1e-12,
    relative = TRUE
  )
})

test_that("frechet_median reference CDF at mu equals 0.5", {
  for (mu in c(0.5, 1, 4)) {
    # a finite median does not require a finite mean: nu <= 1 included
    for (nu in c(0.5, 1, 2, 5)) {
      expect_eps(
        stats::integrate(dfrechet_median, lower = 0, upper = mu,
                         mu = mu, nu = nu)$value,
        0.5,
        eps = 1e-6
      )
    }
  }
})

test_that("frechet_median density integrates to one", {
  for (nu in c(1, 3)) {
    expect_eps(
      stats::integrate(dfrechet_median, lower = 0, upper = Inf,
                       mu = 2, nu = nu)$value,
      1,
      eps = 1e-6
    )
  }
})

test_that("frechet_median quantile checks", {
  # closed form quantile under median parameterization
  qfrechet_median <- function(p, mu, nu) {
    mu * log(2)^(1 / nu) * (-log(p))^(-1 / nu)
  }
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  for (mu in c(0.5, 2)) {
    for (nu in c(0.5, 2, 5)) {
      expect_eps(
        qfrechet_median(p, mu = mu, nu = nu),
        brms::qfrechet(p, loc = 0, scale = mu * log(2)^(1 / nu), shape = nu),
        eps = 1e-12,
        relative = TRUE
      )
    }
  }
})

test_that("frechet_median RNG recovers median and mean", {
  set.seed(1123)
  n <- 1e5
  for (mu in c(0.5, 2)) {
    for (nu in c(0.5, 2, 5)) {
      draws <- rfrechet_median(n, mu = mu, nu = nu)
      # median is mu, also for nu <= 1 (no finite mean required)
      expect_eps(median(draws), mu, eps = 0.03, relative = TRUE)
    }
    for (nu in c(2, 5)) {
      draws <- rfrechet_median(n, mu = mu, nu = nu)
      # for nu > 1 the mean is sigma * Gamma(1 - 1 / nu)
      true_mean <- mu * log(2)^(1 / nu) * gamma(1 - 1 / nu)
      expect_eps(mean(draws), true_mean, eps = 0.05, relative = TRUE)
    }
  }
})

test_that("frechet_median RNG recovers quantiles", {
  set.seed(772)
  mu <- 2
  nu <- 2.5
  draws <- rfrechet_median(2e5, mu = mu, nu = nu)
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  true_quantiles <- mu * log(2)^(1 / nu) * (-log(p))^(-1 / nu)
  expect_eps(
    true_quantiles,
    stats::quantile(draws, probs = p),
    eps = 0.03,
    r = 0.2,
    relative = TRUE
  )
})

test_that("frechet_median argument checks", {
  expect_error(dfrechet_median(-1, mu = 1, nu = 2))
  expect_error(dfrechet_median(1, mu = 0, nu = 2))
  expect_error(dfrechet_median(1, mu = 1, nu = 0))
  # nu <= 1 is now allowed (finite median without finite mean)
  expect_silent(dfrechet_median(1, mu = 1, nu = 0.5))
  expect_error(rfrechet_median(10, mu = 0, nu = 2))
  expect_error(rfrechet_median(10, mu = 1, nu = 0))
  expect_silent(rfrechet_median(10, mu = 1, nu = 0.5))
})
