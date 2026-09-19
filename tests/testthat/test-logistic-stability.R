# Issue #32: dlogistic computed -z - log(sigma) - 2*log1p(exp(-z)) which
# overflowed to -Inf in the left tail even though the log density is
# representable there (stats::dlogis returns it).

test_that("dlogistic equals stats::dlogis on a wide grid, both tails (#32)", {
  z <- c(
    -1e5,
    -3000,
    -1000,
    -745,
    -100,
    -10,
    -1,
    0,
    1,
    10,
    100,
    745,
    1000,
    3000,
    1e5
  )
  for (mu in c(-50, 0, 3.5)) {
    for (s in c(0.1, 1, 7)) {
      x <- mu + s * z
      expect_equal(
        dlogistic(x, mu = mu, sigma = s, log = TRUE),
        stats::dlogis(x, location = mu, scale = s, log = TRUE)
      )
      expect_equal(
        dlogistic(x, mu = mu, sigma = s),
        stats::dlogis(x, location = mu, scale = s)
      )
    }
  }
})

test_that("dlogistic left tail is finite (#32)", {
  expect_equal(dlogistic(-1000, mu = 0, sigma = 1, log = TRUE), -1000)
  expect_equal(dlogistic(1000, mu = 0, sigma = 1, log = TRUE), -1000)
  # asymmetric case from the issue: only the left tail used to break
  expect_true(is.finite(dlogistic(-1e5, mu = -50, sigma = 0.1, log = TRUE)))
})

test_that("dlogistic is symmetric around mu (#32)", {
  mu <- 2
  d <- c(0.5, 1, 10, 100, 800, 5000)
  expect_equal(
    dlogistic(mu + d, mu = mu, sigma = 3, log = TRUE),
    dlogistic(mu - d, mu = mu, sigma = 3, log = TRUE)
  )
  expect_equal(
    dlogistic(mu + d, mu = mu, sigma = 3),
    dlogistic(mu - d, mu = mu, sigma = 3)
  )
})

test_that("dlogistic matches stats::dlogis on a dense random grid (#32)", {
  set.seed(20241001)
  x <- rnorm(500, 0, 50)
  for (mu in c(-3, 0, 4)) {
    for (s in c(0.5, 2, 20)) {
      expect_equal(
        dlogistic(x, mu = mu, sigma = s, log = TRUE),
        stats::dlogis(x, location = mu, scale = s, log = TRUE)
      )
    }
  }
})
