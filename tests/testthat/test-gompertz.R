test_that("custom-gompertz", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
  unit <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = 0.1, to = 50, length.out = n_small) # mu can't be too close to 0 for numerical reasons
  beta_list <- seq(from = eps, to = 1, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_failure_rate <- 5e-4
  accepted_rng_error <- 0.02
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(dgompertz(x, mu = 10, b = 1)))
  expect_equal(n, length(qgompertz(unit, mu = 10, b = 1)))
  expect_equal(n, length(rgompertz(n, 1, 3)))

  # Compare density and quantile functions to extraDistr
  for (mu in mu_list) {
    for (beta in beta_list) {
      expect_eps(
        dgompertz(x, mu = mu, beta = beta),
        extraDistr::dgompertz(x, a = -(beta * log(0.5)) / (exp(mu * beta) - 1), b = beta),
        eps = accepted_relative_error,
        r = accepted_failure_rate,
        relative = TRUE
      )
      expect_eps(
        qgompertz(unit, mu = mu, beta = beta),
        extraDistr::qgompertz(unit, a = -(beta * log(0.5)) / (exp(mu * beta) - 1), b = beta),
        eps = accepted_relative_error,
        r = accepted_failure_rate,
        relative = TRUE
      )
    }
  }

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rgompertz,
    metric_mu = median,
    n = 5 * n,
    mu_list = mu_list,
    aux_list = beta_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check if the RNG can recover the quantiles
  test_rng_quantiles(
    rng_fun = rgompertz,
    quantile_fun = qgompertz,
    n = 50 * n,
    mu_list = mu_list,
    aux_list = beta_list,
    eps = accepted_rng_error,
    quantiles = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check density function for errors
  expect_error(dgompertz(1, 2)) # to few arguments
  expect_error(dgompertz(1, 2, 3, 4, 5)) # to many arguments
  expect_error(dgompertz(-1, mu = 2, beta = 2)) # x is not allowed to be smaller 0
  expect_error(dgompertz(1, mu = 0, beta = 2)) # mu is not allowed to be 0 or smaller
  expect_error(dgompertz(1, mu = 1, beta = 0)) # beta is not allowed to be 0 or smaller
  expect_error(dgompertz("r", mu = 2, beta = 2)) # non-numeric argumen

  # Check quantile function for errors
  expect_error(qgompertz(1, 2)) # to few arguments
  expect_error(qgompertz(1, 2, 3, 4, 5)) # to many arguments
  expect_error(qgompertz(1, mu = 0, beta = 2)) # mu is not allowed to be 0 or smaller
  expect_error(qgompertz(1, mu = 1, beta = 0)) # beta is not allowed to be 0 or smaller
  expect_error(qgompertz(c(-1, 2), mu = 2, beta = 2)) # q is not allowed to be outside [0, 1]
  expect_error(qgompertz("r", mu = 2, beta = 2)) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(rgompertz(100, 2)) # to few arguments
  expect_error(rgompertz(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rgompertz(-1, mu = 2, b = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(rgompertz("r", mu = 2, b = 2))) # non-numeric arguments are disallowed
  expect_error(rgompertz(100, mu = 0, beta = 2)) # mu is not allowed to be 0 or smaller
  expect_error(rgompertz(100, mu = 1, beta = -1)) # beta is not allowed to be 0 or smaller

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = 5,
    aux_par = 2,
    ref_intercept = 5,
    parameter_link = log,
    rng_link = identity,
    family = gompertz,
    rng = rgompertz,
    aux_name = "beta"
  )
})
