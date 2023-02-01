test_that("test-custom-beta", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-2
  x <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = eps, to = 1 - eps, length.out = n_small)
  phi_list <- seq(from = 0.001, to = 10, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.15
  warning("current accepted rng_error is high at 0.15")
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(dbeta_mean(x, mu = 0.8, phi = 2)))
  expect_equal(n, length(rbeta_mean(n, mu = 0.8, phi = 2)))

  # Compare density and to built-in
  for (mu in mu_list) {
    for (phi in phi_list) {
      expect_eps(
        dbeta_mean(x, mu = mu, phi = phi),
        dbeta(x, shape1 = mu * phi, shape2 = (1 - mu) * phi),
        eps = accepted_relative_error,
        relative = TRUE
      )
    }
  }

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rbeta_mean,
    metric_mu = mean,
    n = 5 * n,
    mu_list = mu_list,
    aux_list = phi_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check if the RNG can recover the quantiles
  test_rng_quantiles(
    rng_fun = rbeta_mean,
    quantile_fun = qbeta_mean,
    n = 50 * n,
    mu_list = mu_list,
    aux_list = phi_list,
    eps = accepted_rng_error * 3,
    quantiles = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),
    p_acceptable_failures = accepred_rng_failures * 2,
    relative = TRUE
  )


  # Check density function for errors
  expect_error(dbeta_mean(1, 0.8)) # to few arguments
  expect_error(dbeta_mean(1, 0.8, 3, 4, 5)) # to many arguments
  expect_error(dbeta_mean(-1, mu = 0.8, phi = 2)) # x is not allowed to be smaller 0
  expect_error(dbeta_mean("r", mu = 0.8, phi = 2)) # non-numeric arguments are disallowed
  expect_error(dbeta_mean(1, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  expect_error(dbeta_mean(1, mu = 0.8, phi = 2)) # mu is not allowed to be 1 or bigger
  expect_error(dbeta_mean(1, mu = 0.8, phi = 0)) # p is not allowed to be 1 or smaller

  # Check rng for errors
  expect_error(rbeta_mean(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rbeta_mean(-1, mu = 0.8, phi = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(rbeta_mean("r", mu = 0.8, phi = 2))) # non-numeric arguments are disallowed
  expect_error(rbeta_mean(100, mu = 1, phi = 2)) # mu must be between 0 and 1
  expect_error(rbeta_mean(100, mu = 0.8, phi = -1)) # phi is not allowed to be negative
})
