# Issue #39: posterior_epred_symlognormal called warning() and returned the
# warning text instead of predictions. It now evaluates the closed-form mean
# E[Y] = exp(mu + s^2/2) Phi((mu + s^2)/s)
#      - exp(-mu + s^2/2) Phi((-mu + s^2)/s) + 1 - 2 Phi(mu/s)
# with stable log-space tail arithmetic.

quad_ref <- function(mu, sigma) {
  integrate(
    function(z) sign(z) * expm1(abs(z)) * dnorm(z, mu, sigma),
    mu - 40 * sigma,
    mu + 40 * sigma,
    rel.tol = 1e-12
  )$value
}

make_prep <- function(mu_obs, sigma_obs, ndraws) {
  mu <- matrix(rep(mu_obs, each = ndraws), nrow = ndraws)
  sigma <- matrix(rep(sigma_obs, each = ndraws), nrow = ndraws)
  structure(
    list(dpars = list(mu = mu, sigma = sigma), ndraws = ndraws),
    class = "brmsprep"
  )
}

test_that("symlognormal mean agrees with quadrature on a moderate grid (#39)", {
  mu_values <- c(-5, -2, -1, 0, 1, 2, 5)
  sigma_values <- c(0.01, 0.1, 0.5, 1, 2, 3)
  for (mu in mu_values) {
    for (sigma in sigma_values) {
      ref <- quad_ref(mu, sigma)
      expect_equal(symlognormal_mean(mu, sigma), ref, tolerance = 1e-9)
    }
  }
  # example from the issue
  expect_equal(symlognormal_mean(1, 2), 17.04121, tolerance = 1e-4)
})

test_that("symlognormal mean is 0 at mu = 0 and symmetric (#39)", {
  for (sigma in c(0.1, 0.5, 2, 10, 50)) {
    expect_equal(symlognormal_mean(0, sigma), 0)
  }
  for (mu in c(-5, -1, 0.5, 2)) {
    for (sigma in c(0.1, 1, 3)) {
      expect_equal(
        symlognormal_mean(-mu, sigma),
        -symlognormal_mean(mu, sigma),
        tolerance = 1e-12
      )
    }
  }
})

test_that("symlognormal mean approaches the degenerate limit as sigma -> 0 (#39)", {
  expect_equal(symlognormal_mean(2, 1e-6), exp(2) - 1, tolerance = 1e-12)
  expect_equal(symlognormal_mean(-2, 1e-6), 1 - exp(2), tolerance = 1e-12)
  expect_equal(symlognormal_mean(0, 1e-6), 0, tolerance = 1e-12)
})

test_that("symlognormal mean stays finite or signed-Inf, never NaN (#39)", {
  # extreme sigma with asymmetric mu: the two closed-form terms each overflow
  # double precision, but their log-space difference must not produce NaN
  expect_false(is.nan(symlognormal_mean(0.5, 50)))
  expect_false(is.nan(symlognormal_mean(-0.5, 50)))
  # exactly symmetric case is exactly zero despite overflowing terms
  expect_equal(symlognormal_mean(0, 50), 0)
})

test_that("posterior_epred_symlognormal returns a numeric prediction matrix (#39)", {
  ndraws <- 4
  mu_obs <- c(-1, 0, 1.5)
  sigma_obs <- c(0.5, 2, 1)
  prep <- make_prep(mu_obs, sigma_obs, ndraws)
  out <- posterior_epred_symlognormal(prep)
  expect_silent(posterior_epred_symlognormal(prep))
  expect_type(out, "double")
  expect_equal(dim(out), c(ndraws, length(mu_obs)))
  for (j in seq_along(mu_obs)) {
    expect_equal(
      out[, j],
      rep(symlognormal_mean(mu_obs[j], sigma_obs[j]), ndraws),
      tolerance = 1e-12
    )
  }
})
