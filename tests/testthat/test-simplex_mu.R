test_that("simplex d/r share the strict open-interval mu policy", {
  bad_mu <- c(0, 1, -0.1, 1.1)
  for (mu in bad_mu) {
    expect_error(dsimplex(0.5, mu = mu, sigma = 2), regexp = "mean must be")
    expect_error(rsimplex(10, mu = mu, sigma = 2), regexp = "mean must be")
  }
})

test_that("simplex no longer silently clamps mu near the boundary", {
  # values inside the open interval are used as-is, matching the Stan
  # simplex_lpdf, which also uses mu as-is
  x <- c(0.2, 0.5, 0.8)
  unclamped_logpdf <- function(x, mu, sigma) {
    (-0.5) *
      (log(2) + log(pi) + 2 * log(sigma) + 3 * (log(x) + log1p(-x))) +
      ((-1 / (2 * sigma^2)) *
        (((x - mu)^2) / (x * (1 - x) * mu^2 * (1 - mu)^2)))
  }
  for (mu in 1 - 10^-(6:12)) {
    expect_eps(
      dsimplex(x, mu = mu, sigma = 2, log = TRUE),
      unclamped_logpdf(x, mu, 2),
      eps = 1e-12,
      relative = TRUE
    )
  }
  # reference value from issue #42 (unclamped Stan result at mu = 1 - 1e-10)
  expect_eps(
    dsimplex(0.5, mu = 1 - 1e-10, sigma = 2, log = TRUE),
    -1.25e19,
    eps = 0.1,
    relative = TRUE
  )
})
