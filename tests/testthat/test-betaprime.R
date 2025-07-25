test_that("custom-betaprime", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
  unit <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = eps, to = 100, length.out = n_small)
  phi_list <- seq(from = eps, to = 50, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.085
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(dbetaprime(x, mu = 1, phi = 2)))
  expect_equal(n, length(qbetaprime(unit, mu = 1, phi = 2)))
  expect_equal(n, length(rbetaprime(n, 2, 3)))

  # Compare density and quantile functions to extraDistr
  for (mu in mu_list) {
    for (phi in phi_list) {
      expect_eps(
        dbetaprime(x, mu = mu, phi = phi),
        extraDistr::dbetapr(x, mu * (phi + 1), phi + 2),
        eps = accepted_relative_error,
        relative = TRUE
      )
      expect_eps(
        qbetaprime(unit, mu = mu, phi = phi),
        extraDistr::qbetapr(unit, mu * (phi + 1), phi + 2),
        eps = accepted_relative_error,
        relative = TRUE
      )
    }
  }

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rbetaprime,
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
    rng_fun = rbetaprime,
    quantile_fun = qbetaprime,
    n = 5 * n,
    mu_list = mu_list,
    aux_list = phi_list,
    eps = accepted_rng_error,
    quantiles = c(
      0.01,
      0.05,
      0.1,
      0.2,
      0.3,
      0.4,
      0.5,
      0.6,
      0.7,
      0.8,
      0.9,
      0.95,
      0.99
    ),
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check density function for errors
  expect_error(dbetaprime(1, 2)) # to few arguments
  expect_error(dbetaprime(1, 2, 3, 4, 5)) # to many arguments
  expect_error(dbetaprime(-1, mu = 2, phi = 2)) # x is not allowed to be smaller 0
  expect_error(dbetaprime(1, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  expect_error(dbetaprime(1, mu = 1, phi = 0)) # phi is not allowed to be 0 or smaller
  # expect_error(dbetaprime("r", mu = 2, phi = 2)) # non-numeric arguments are disallowed

  # Check quantile function for errors
  expect_error(qbetaprime(1, 2)) # to few arguments
  expect_error(qbetaprime(1, 2, 3, 4, 5)) # to many arguments
  expect_error(qbetaprime(1, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  expect_error(qbetaprime(1, mu = 1, phi = 0)) # phi is not allowed to be 0 or smaller
  expect_error(qbetaprime(c(-1, 2), mu = 2, phi = 2)) # q is not allowed to be outside [0, 1]
  # expect_error(qbetaprime("r", mu = 2, phi = 2)) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(rbetaprime(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rbetaprime(-1, mu = 2, phi = 2)) # number of drawn samples cannot be smaller 0
  # expect_warning(expect_error(rbetaprime("r", mu = 2, phi = 2))) # non-numeric arguments are disallowed
  expect_error(rbetaprime(100, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  expect_error(rbetaprime(100, mu = 1, phi = 0)) # phi is not allowed to be 0 or smaller

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = 5,
    aux_par = 2,
    ref_intercept = 5,
    rng_link = identity,
    parameter_link = log,
    family = betaprime,
    rng = rbetaprime,
    aux_name = "phi"
  )
})
