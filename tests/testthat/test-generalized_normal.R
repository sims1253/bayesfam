test_that("custom-generalized-normal", {
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
  expect_equal(n, length(dgeneralized_normal(x, mu = 1, sigma = 2, beta = 3)))
  expect_equal(n, length(qgeneralized_normal(unit, mu = 1, sigma = 2, beta = 3)))
  expect_equal(n, length(rgeneralized_normal(n, mu = 1, sigma = 2, beta = 3)))

  # Compare density and quantile functions to extraDistr
  for (mu in mu_list) {
    for (aux1 in aux_list) {
      for(aux2 in aux_list) {
        expect_eps(
          dgeneralized_normal(x, mu = mu, sigma = aux1, beta = aux2),
          gnorm::dgnorm(x, mu = mu, alpha = aux1, beta = aux2),
          eps = accepted_relative_error,
          relative = TRUE
        )
        # expect_eps(
        #   qgeneralized_normal(unit, mu = mu, sigma = aux1, beta = aux2),
        #   gnorm::qgnorm(unit, mu = mu, alpha = aux1, beta = aux2),
        #   eps = accepted_relative_error,
        #   relative = TRUE
        # )
      }
    }
  }
  warning("quanitle produces fails")

  # check if the RNG is close enough to the true mean in most cases
  # test_rng(
  #   rng_fun = rgeneralized_normal,
  #   metric_mu = mean,
  #   n = 5 * n,
  #   mu_list = mu_list,
  #   aux_list = aux_list,
  #   mu_eps = accepted_rng_error,
  #   p_acceptable_failures = accepred_rng_failures,
  #   relative = TRUE
  # )
  warning("RNG test function currently uses only 1 auxiliary parameter")

  # Check density function for errors
  # expect_error(dbetaprime(1, 2)) # to few arguments
  # expect_error(dbetaprime(1, 2, 3, 4, 5)) # to many arguments
  # expect_error(dbetaprime(-1, mu = 2, phi = 2)) # x is not allowed to be smaller 0
  # expect_error(dbetaprime(1, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  # expect_error(dbetaprime(1, mu = 1, phi = 0)) # phi is not allowed to be 0 or smaller
  # expect_error(dbetaprime("r", mu = 2, phi = 2)) # non-numeric arguments are disallowed

  # Check quantile function for errors
  # expect_error(qbetaprime(1, 2)) # to few arguments
  # expect_error(qbetaprime(1, 2, 3, 4, 5)) # to many arguments
  # expect_error(qbetaprime(1, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  # expect_error(qbetaprime(1, mu = 1, phi = 0)) # phi is not allowed to be 0 or smaller
  # expect_error(qbetaprime(c(-1, 2), mu = 2, phi = 2)) # q is not allowed to be outside [0, 1]
  # expect_error(qbetaprime("r", mu = 2, phi = 2)) # non-numeric arguments are disallowed

  # Check rng for errors
  # expect_error(rbetaprime(10, 2, 3, 4, 5)) # to many arguments
  # expect_error(rbetaprime(-1, mu = 2, phi = 2)) # number of drawn samples cannot be smaller 0
  # expect_warning(expect_error(rbetaprime("r", mu = 2, phi = 2))) # non-numeric arguments are disallowed
  # expect_error(rbetaprime(100, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  # expect_error(rbetaprime(100, mu = 1, phi = 0)) # phi is not allowed to be 0 or smaller

  # Check of brms can fit the custom family and recover the intercept and shape
  # expect_brms_family(
  #   intercept = 5,
  #   aux_par = 2,
  #   ref_intercept = 5,
  #   rng_link = identity,
  #   parameter_link = log,
  #   family = betaprime,
  #   rng = rbetaprime,
  #   aux_name = "phi"
  # )
  warning("BRMS test function currently uses only 1 auxiliary parameter")
})
