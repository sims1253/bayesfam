test_that("Generalized Gamma", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- exp(seq(from = eps, to = 200, length.out = n))
  unit <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = -1, to = 100, length.out = n_small)
  sigma_list <- seq(from = 0.1, to = 100, length.out = n_small)
  Q_list <- seq(from = -1, to = 100, length.out = n_small)
  Q_list[2] <- 0

  accepted_relative_error <- 1e-6

  # Check lengths
  expect_equal(n, length(dgeneralized_gamma(x, mu = 1, sigma = 2, Q = 2)))
  expect_equal(n, length(rgeneralized_gamma(n, mu = 1, sigma = 2, Q = 2)))

  # Test density function

  for (mu in mu_list) {
    for (sigma in sigma_list) {
      for (Q in Q_list) {
        expect_eps(
          dgeneralized_gamma(x, mu = mu, sigma = sigma, Q = Q),
          flexsurv::dgengamma(x, mu = mu, sigma = sigma, Q = Q),
          eps = accepted_relative_error,
          relative = TRUE
        )
      }
    }
  }

  # Currently not checking the rng

  # Check density function for errors
  expect_error(dgeneralized_gamma(1, 2)) # to few arguments
  expect_error(dgeneralized_gamma(1, 2, 3, 4, 5, 6)) # to many arguments
  expect_error(dgeneralized_gamma(-1, mu = 2, sigma = 2, Q = 2)) # x is not allowed to be smaller 0
  expect_error(dgeneralized_gamma(-1, mu = 2, sigma = -2, Q = 2)) # sigma is not allowed to be 0 or smaller

  # Check rng for errors
  expect_error(rgeneralized_gamma(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rgeneralized_gamma(-1, mu = 2, sigma = 2, Q = 2)) # number of drawn samples cannot be smaller 0
  # expect_warning(expect_error(generalized_gamma("r", mu = 2, phi = 2))) # non-numeric arguments are disallowed
  expect_error(rgeneralized_gamma(100, mu = 2, sigma = 0, Q = 2)) # sigma is not allowed to be 0 or smaller

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = 5,
    aux_par = 2,
    aux2_par = 3,
    ref_intercept = 5,
    rng_link = identity,
    parameter_link = log,
    family = generalized_gamma,
    rng = rgeneralized_gamma,
    aux_name = "sigma",
    aux2_name = "Q"
  )
})

test_that("dgeneralized_gamma is vectorized in Q (#23)", {
  x <- c(0.5, 1.2, 3.7, 10)
  mu <- c(0, 0.2, -0.5, 1.1)
  sigma <- c(1, 0.7, 2.3, 0.5)
  for (Q in list(c(0, 0, 0, 0), c(-1.5, 0.3, 2, 7.7), c(0, 0.9, 0, -2.2))) {
    expect_equal(
      dgeneralized_gamma(x, mu, sigma, Q = Q, log = TRUE),
      mapply(dgeneralized_gamma, x, mu, sigma, Q, log = TRUE)
    )
  }

  # log_lik callback: scalar y with length-S draws of Q
  Q <- c(0, 1.1, 0, -0.7, 2.3)
  prep <- structure(
    list(
      dpars = list(mu = rep(0.3, 5), sigma = rep(1.2, 5), Q = Q),
      data = list(Y = c(1.7, 2.5))
    ),
    class = "brmsprep"
  )
  expect_equal(
    log_lik_generalized_gamma(2, prep),
    mapply(dgeneralized_gamma, 2.5, 0.3, 1.2, Q, log = TRUE)
  )
})

# Reference mean without gamma functions: with G ~ Gamma(1 / Q^2, 1) and
# v = log(Q^2 G) / Q, Y = exp(mu + sigma v), and v has a density
# proportional to exp(-v^2 h(Q v)), h(x) = (exp(x) - 1 - x) / x^2.
generalized_gamma_reference_mean <- function(mu, sigma, Q) {
  h <- function(x) {
    out <- (expm1(x) - x) / x^2
    small <- abs(x) < 0.1
    # series sum_j x^j / (j + 2)!
    term <- 1 / 2
    total <- 0
    for (j in 0:14) {
      total <- total + term
      term <- term * x[small] / (j + 3)
    }
    out[small] <- total
    out
  }
  both_sides <- function(f, at) {
    integrate(f, -Inf, at, rel.tol = 1e-13, subdivisions = 5000L)$value +
      integrate(f, at, Inf, rel.tol = 1e-13, subdivisions = 5000L)$value
  }
  exp(mu) *
    both_sides(function(v) exp(sigma * v - v^2 * h(Q * v)), sigma) /
    both_sides(function(v) exp(-v^2 * h(Q * v)), 0)
}

test_that("the generalized gamma mean stays accurate as Q approaches 0", {
  Q <- c(1e-200, 1e-9, 1e-7, 1e-5, 1e-3, 0.05, 0.3, 0.32, 0.8)
  Q <- c(Q, -Q)
  for (p in list(c(mu = 0, sigma = 1), c(mu = 0.3, sigma = 0.5))) {
    for (q in Q) {
      expect_equal(
        generalized_gamma_mean(p[["mu"]], p[["sigma"]], q),
        generalized_gamma_reference_mean(p[["mu"]], p[["sigma"]], q),
        tolerance = 1e-12,
        label = sprintf("mean at mu = %g, sigma = %g, Q = %g", p[["mu"]], p[["sigma"]], q)
      )
    }
  }
})

test_that("the generalized gamma mean is Inf exactly where Q * sigma <= -1", {
  # Q * sigma is exactly -1 here; 1 / Q^2 + sigma / Q rounds to 4e-16
  expect_identical(generalized_gamma_mean(0, 1 / 0.7, -0.7), Inf)
  expect_true(is.finite(generalized_gamma_mean(0, 1.4285, -0.7)))
})
