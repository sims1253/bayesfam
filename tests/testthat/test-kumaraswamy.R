test_that("custom-kumaraswamy", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- seq(from = eps, to = 1 - eps, length.out = n)
  unit <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = eps + 0.075, to = 1 - eps, length.out = n_small)
  p_list <- seq(from = 0.1, to = 10, length.out = n_small)
  accepted_relative_error <- 1e-5
  accepted_rng_error <- 0.1
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(dkumaraswamy(x, mu = 0.5, p = 2)))
  expect_equal(n, length(qkumaraswamy(unit, mu = 0.5, p = 2)))
  expect_equal(n, length(rkumaraswamy(n, mu = 0.5, p = 2)))

  # Compare density and quantile functions to extraDistr
  for (mu in mu_list) {
    for (p in p_list) {
      expect_eps(
        dkumaraswamy(x, mu = mu, p = p),
        extraDistr::dkumar(x, a = p, b = -log(2) / log1p(-mu^p)),
        eps = accepted_relative_error,
        relative = TRUE
      )
      expect_eps(
        qkumaraswamy(unit, mu = mu, p = p),
        extraDistr::qkumar(unit, a = p, b = -log(2) / log1p(-mu^p)),
        eps = accepted_relative_error,
        relative = TRUE
      )
    }
  }

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rkumaraswamy,
    metric_mu = median,
    n = 5 * n,
    mu_list = mu_list,
    aux_list = p_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check if the RNG can recover the quantiles
  test_rng_quantiles(
    rng_fun = rkumaraswamy,
    quantile_fun = qkumaraswamy,
    n = 50 * n,
    mu_list = mu_list,
    aux_list = p_list,
    eps = accepted_rng_error,
    quantiles = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check density function for errors
  expect_error(dkumaraswamy(1, 0.8)) # to few arguments
  expect_error(dkumaraswamy(1, 0.8, 3, 4, 5)) # to many arguments
  expect_error(dkumaraswamy(-1, mu = 0.8, p = 2)) # x is not allowed to be smaller 0
  expect_error(dkumaraswamy("r", mu = 0.8, p = 2)) # non-numeric arguments are disallowed
  expect_error(dkumaraswamy(1, mu = 0, p = 2)) # mu is not allowed to be 0 or smaller
  expect_error(dkumaraswamy(1, mu = 0.8, p = 2)) # mu is not allowed to be 1 or bigger
  expect_error(dkumaraswamy(1, mu = 0.8, p = 0)) # p is not allowed to be 1 or smaller

  # Check quantile function for errors
  expect_error(qkumaraswamy(1, 0.8)) # to few arguments
  expect_error(qkumaraswamy(1, 0.8, 3, 4, 5)) # to many arguments
  expect_error(qkumaraswamy("r", mu = 0.8, p = 2)) # non-numeric arguments are disallowed
  expect_error(qkumaraswamy(1, mu = 0, p = 2)) # mu is not allowed to be 0 or smaller
  expect_error(qkumaraswamy(1, mu = 1, p = 2)) # mu is not allowed to be 1 or bigger
  expect_error(qkumaraswamy(1, mu = 0.8, p = 0)) # p is not allowed to be 0 or smaller
  expect_error(qkumaraswamy(c(-1, 2), mu = 2, p = 2)) # q is not allowed to be outside [0, 1]

  # Check rng for errors
  expect_error(rkumaraswamy(100, 0.8)) # to few arguments
  expect_error(rkumaraswamy(100.82, 3, 4, 5)) # to many arguments
  expect_error(rkumaraswamy(-1, mu = 0.8, p = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(rkumaraswamy("r", mu = 0.8, p = 2))) # non-numeric arguments are disallowed
  expect_error(rkumaraswamy(100, mu = 0, p = 2)) # mu is not allowed to be 0 or smaller
  expect_error(rkumaraswamy(100, mu = 1, p = 2)) # mu is not allowed to be 1 or bigger
  expect_error(rkumaraswamy(100, mu = 0.8, p = 0)) # p is not allowed to be 0 or smaller

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = 0.5,
    aux_par = 2,
    ref_intercept = 0.5,
    rng_link = identity,
    parameter_link = logit,
    family = kumaraswamy,
    rng = rkumaraswamy,
    aux_name = "p"
  )
})
