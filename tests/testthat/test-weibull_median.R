test_that("weibull_median density matches stats::dweibull with median scale", {
  mu <- 2
  k <- 3
  x <- seq(from = 0.01, to = 20, length.out = 100)
  expect_eps(
    dweibull_median(x, mu = mu, k = k),
    stats::dweibull(x, shape = k, scale = mu / log(2)^(1 / k)),
    eps = 1e-12,
    relative = TRUE
  )
  expect_eps(
    dweibull_median(x, mu = mu, k = k, log = TRUE),
    stats::dweibull(x, shape = k, scale = mu / log(2)^(1 / k), log = TRUE),
    eps = 1e-12,
    relative = TRUE
  )
})

test_that("weibull_median reference CDF at mu equals 0.5", {
  for (mu in c(0.5, 1, 2, 7.3)) {
    for (k in c(0.5, 1, 2, 5)) {
      expect_eps(
        stats::integrate(dweibull_median, lower = 0, upper = mu,
                         mu = mu, k = k)$value,
        0.5,
        eps = 1e-6
      )
    }
  }
})

test_that("weibull_median density integrates to one", {
  for (k in c(0.7, 1, 3)) {
    expect_eps(
      stats::integrate(dweibull_median, lower = 0, upper = Inf,
                       mu = 2, k = k)$value,
      1,
      eps = 1e-6
    )
  }
})

test_that("weibull_median quantile checks", {
  # closed form quantile under median parameterization
  qweibull_median <- function(p, mu, k) {
    mu * (-log1p(-p) / log(2))^(1 / k)
  }
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  for (mu in c(0.5, 2)) {
    for (k in c(0.7, 1, 3)) {
      expect_eps(
        qweibull_median(p, mu = mu, k = k),
        stats::qweibull(p, shape = k, scale = mu / log(2)^(1 / k)),
        eps = 1e-12,
        relative = TRUE
      )
    }
  }
})

test_that("weibull_median RNG recovers median and mean", {
  set.seed(2323)
  n <- 1e5
  for (mu in c(0.5, 2)) {
    for (k in c(0.7, 1, 3)) {
      draws <- rweibull_median(n, mu = mu, k = k)
      # median is mu
      expect_eps(median(draws), mu, eps = 0.02, relative = TRUE)
      # mean is sigma * Gamma(1 + 1 / k)
      true_mean <- (mu / log(2)^(1 / k)) * gamma(1 + 1 / k)
      expect_eps(mean(draws), true_mean, eps = 0.03, relative = TRUE)
    }
  }
})

test_that("weibull_median RNG recovers quantiles", {
  set.seed(913)
  mu <- 2
  k <- 1.5
  draws <- rweibull_median(2e5, mu = mu, k = k)
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  true_quantiles <- mu * (-log1p(-p) / log(2))^(1 / k)
  expect_eps(
    true_quantiles,
    stats::quantile(draws, probs = p),
    eps = 0.03,
    r = 0.2,
    relative = TRUE
  )
})

test_that("weibull_median argument checks", {
  expect_error(dweibull_median(-1, mu = 1, k = 1))
  expect_error(dweibull_median(1, mu = 0, k = 1))
  expect_error(dweibull_median(1, mu = 1, k = 0))
  expect_error(rweibull_median(10, mu = 0, k = 1))
  expect_error(rweibull_median(10, mu = 1, k = 0))
})
