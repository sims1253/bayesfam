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
