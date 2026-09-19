test_that("kumaraswamy d/p/q/r share the strict open-interval mu policy", {
  # mu outside the open interval (0, 1) errors in every function
  bad_mu <- c(0, 1, -0.1, 1.1)
  for (mu in bad_mu) {
    expect_error(dkumaraswamy(0.3, mu = mu, p = 2), regexp = "median must be")
    expect_error(qkumaraswamy(0.3, mu = mu, p = 2), regexp = "median must be")
    expect_error(pkumaraswamy(0.3, mu = mu, p = 2), regexp = "median must be")
    expect_error(rkumaraswamy(10, mu = mu, p = 2), regexp = "median must be")
  }
})

test_that("kumaraswamy no longer silently clamps mu near the boundary", {
  # values inside the open interval are used as-is, so density, CDF and
  # quantile agree with the unclamped reference (and hence with the
  # Stan kumaraswamy_lpdf, which also uses mu as-is)
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  p <- c(0.05, 0.3, 0.5, 0.7, 0.95)
  for (mu in 1 - 10^-(6:12)) {
    q <- -log(2) / log1p(-mu^2)
    # unclamped density reference
    expect_eps(
      dkumaraswamy(x, mu = mu, p = 2, log = TRUE),
      log(2) + log(log(2)) - log(-(log1p(-mu^2))) + log(x) +
        (q - 1) * log1p(-x^2),
      eps = 1e-12,
      relative = TRUE
    )
    expect_eps(
      dkumaraswamy(x, mu = mu, p = 2),
      extraDistr::dkumar(x, a = 2, b = q),
      eps = 1e-10,
      relative = TRUE
    )
    # q stays consistent with the unclamped density at the same extreme mu
    expect_eps(
      qkumaraswamy(p, mu = mu, p = 2),
      extraDistr::qkumar(p, a = 2, b = q),
      eps = 1e-10,
      relative = TRUE
    )
  }
  # the previously documented clamped value at mu = 1 - 1e-10 must now
  # differ from the clamped result (i.e. clamping is really gone)
  mu <- 1 - 1e-10
  expect_eps(
    dkumaraswamy(0.3, mu = mu, p = 2, log = TRUE),
    -3.892007,
    eps = 1e-5,
    relative = TRUE
  )
})

test_that("kumaraswamy RNG uses unclamped mu", {
  set.seed(42)
  mu <- 1 - 1e-8
  draws <- rkumaraswamy(2e5, mu = mu, p = 2)
  # sampled median recovers the extreme mu (median parameterization)
  expect_eps(median(draws), mu, eps = 1e-4, relative = TRUE)
})
