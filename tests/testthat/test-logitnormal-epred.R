# Issue #34: posterior_epred_logitnormal returned plogis(mu) (the median) and
# ignored sigma. It now computes E[plogis(Z)], Z ~ N(mu, sigma), via
# deterministic quadrature, keeping draw-by-observation dimensions.

quad_ref <- function(mu, sigma) {
  integrate(
    function(t) plogis(mu + sigma * t) * dnorm(t),
    -Inf,
    Inf,
    rel.tol = 1e-13
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

test_that("posterior_epred_logitnormal computes the mean, not the median (#34)", {
  # issue example: actual mean 0.6477264, median/old value 0.7310586
  expect_equal(logitnormal_mean(1, 2), 0.6477264, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(
    logitnormal_mean(1, 2),
    plogis(1)
  )))
})

test_that("logitnormal mean agrees with quadrature for several (mu, sigma) (#34)", {
  mu_values <- c(-5, -2, -1, 0, 1, 2, 5, 10)
  sigma_values <- c(0.01, 0.1, 0.5, 1, 2, 5, 10, 20)
  for (mu in mu_values) {
    for (sigma in sigma_values) {
      ref <- quad_ref(mu, sigma)
      expect_equal(logitnormal_mean(mu, sigma), ref, tolerance = 1e-9)
    }
  }
})

test_that("logitnormal mean agrees with Monte Carlo (#34)", {
  set.seed(20241003)
  n <- 2e6
  for (par in list(c(0, 1), c(1, 2), c(-2, 0.5), c(3, 5))) {
    draws <- rlogitnormal(n, mu = par[1], sigma = par[2])
    # MC standard error of a mean over (0, 1) is below 3e-4
    expect_equal(
      logitnormal_mean(par[1], par[2]),
      mean(draws),
      tolerance = 3e-3
    )
  }
})

test_that("logitnormal mean depends on sigma correctly (#34)", {
  mu <- 1
  median <- plogis(mu)
  small <- logitnormal_mean(mu, 0.05)
  medium <- logitnormal_mean(mu, 0.5)
  large <- logitnormal_mean(mu, 2)
  # mean is pulled from the median towards 0.5 as sigma grows
  expect_lt(small, median)
  expect_lt(medium, small)
  expect_lt(large, medium)
  expect_gt(large, 0.5)
  # sigma -> 0 recovers plogis(mu)
  expect_equal(logitnormal_mean(mu, 1e-9), median)
})

test_that("logitnormal mean is symmetric at mu = 0 (#34)", {
  for (sigma in c(0.1, 1, 5, 20)) {
    expect_equal(logitnormal_mean(0, sigma), 0.5, tolerance = 1e-12)
    expect_equal(
      logitnormal_mean(1.5, sigma),
      1 - logitnormal_mean(-1.5, sigma),
      tolerance = 1e-12
    )
  }
})

test_that("posterior_epred_logitnormal keeps draw-by-observation dimensions (#34)", {
  ndraws <- 5
  mu_obs <- c(-2, 0, 1)
  sigma_obs <- c(0.5, 2, 10)
  prep <- make_prep(mu_obs, sigma_obs, ndraws)
  out <- expect_silent(posterior_epred_logitnormal(prep))
  expect_type(out, "double")
  expect_equal(dim(out), c(ndraws, length(mu_obs)))
  for (j in seq_along(mu_obs)) {
    expect_equal(
      out[, j],
      rep(logitnormal_mean(mu_obs[j], sigma_obs[j]), ndraws),
      tolerance = 1e-12
    )
  }
})
