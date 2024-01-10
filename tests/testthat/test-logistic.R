test_that("custom-logistic", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- exp(seq(from = -100, to = 100, length.out = n))
  unit <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = -50, to = 50, length.out = n_small)
  aux_list <- seq(from = eps, to = 50, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.085
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(dlogistic(x, mu = 1, sigma = 2)))
  expect_equal(n, length(qlogistic(unit, mu = 1, sigma = 2)))
  expect_equal(n, length(rlogistic(n, mu = 1, sigma = 2)))

  # Compare density and quantile functions to extraDistr
  for (mu in mu_list) {
    for (aux in aux_list) {
      expect_eps(
        dlogistic(x, mu = mu, sigma = aux),
        stats::dlogis(x, mu, aux),
        eps = accepted_relative_error,
        relative = TRUE
      )
      expect_eps(
        qlogistic(unit, mu = mu, sigma = aux),
        stats::qlogis(unit, mu, aux),
        eps = accepted_relative_error,
        relative = TRUE
      )
    }
  }


  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rlogistic,
    metric_mu = mean,
    n = 5 * n,
    mu_list = mu_list,
    aux_list = aux_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check if the RNG can recover the quantiles
  test_rng_quantiles(
    rng_fun = rlogistic,
    quantile_fun = qlogistic,
    n = 5 * n,
    mu_list = mu_list,
    aux_list = aux_list,
    eps = accepted_rng_error,
    quantiles = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check density function for errors
  expect_error(dlogistic(1, 2)) # to few arguments
  expect_error(dlogistic(1, 2, 3, 4, 5)) # to many arguments
  expect_error(dlogistic(1, mu = 1, sigma = 0)) # aux is not allowed to be 0 or smaller
  # expect_error(dlogistic("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed

  # Check quantile function for errors
  expect_error(qlogistic(1, 2)) # to few arguments
  expect_error(qlogistic(1, 2, 3, 4, 5)) # to many arguments
  expect_error(qlogistic(1, mu = 1, sigma = 0)) # aux is not allowed to be 0 or smaller
  expect_error(qlogistic(c(-1, 2), mu = 2, sigma = 2)) # q is not allowed to be outside [0, 1]
  # expect_error(qlogistic("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(rlogistic(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rlogistic(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  # expect_error(rlogistic("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed
  expect_error(rlogistic(100, mu = 1, sigma = 0)) # phi is not allowed to be 0 or smaller


  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = 5,
    aux_par = 2,
    ref_intercept = 5,
    rng_link = identity,
    parameter_link = identity,
    family = logistic,
    rng = rlogistic,
    aux_name = "sigma"
  )
})
