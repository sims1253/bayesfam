test_that("custom-symlognormal", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-3
  x <- exp(seq(from = -100, to = 100, length.out = n))
  unit <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = -50, to = 50, length.out = n_small)
  aux_list <- seq(from = eps, to = 50, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.2
  accepred_rng_failures <- 0.2
  # symlognormal RNG is not that stable

  # Check lengths
  expect_equal(n, length(dsymlognormal(x, mu = 1, sigma = 2)))
  expect_equal(n, length(rsymlognormal(n, mu = 1, sigma = 2)))

  # no reference Density!

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rsymlognormal,
    metric_mu = median,
    n = 5 * n,
    mu_list = mu_list,
    aux_list = aux_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE,
    mu_link = symlog
  )

  # Check density function for errors
  expect_error(dsymlognormal(1, 2)) # to few arguments
  expect_error(dsymlognormal(1, 2, 3, 4, 5)) # to many arguments
  expect_error(dsymlognormal(1, mu = 1, sigma = -1)) # aux is not allowed to be smaller 0
  expect_error(dsymlognormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(rsymlognormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rsymlognormal(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(rsymlognormal("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed
  expect_error(rsymlognormal(100, mu = 1, sigma = -1)) # phi is not allowed to be smaller 0

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = 5,
    aux_par = 2,
    ref_intercept = 5,
    rng_link = identity,
    parameter_link = identity,
    family = symlognormal,
    rng = rsymlognormal,
    aux_name = "sigma"
  )
})
