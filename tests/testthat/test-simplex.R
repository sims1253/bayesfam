test_that("custom-simplex", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = eps, to = 1 - eps, length.out = n_small)
  sigma_list <- seq(from = eps, to = 10, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.05
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(dsimplex(x, mu = 0.5, sigma = 2)))
  expect_equal(n, length(rsimplex(n, mu = 0.5, sigma = 2)))

  # Compare density to reference implementation
  for (mu in mu_list) {
    for (sigma in sigma_list) {
      expect_eps(
        dsimplex(x, mu = mu, sigma = sigma),
        rmutil::dsimplex(x, m = mu, s = sigma^2),
        eps,
        relative = TRUE
      )
    }
  }

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rsimplex,
    metric_mu = mean,
    n = n,
    mu_list = mu_list,
    aux_list = sigma_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check if the RNG can recover the quantiles
  warning("No quantile function available to test rng quantile recovery.")

  # Check density function for errors
  expect_error(dsimplex(1, 0.8)) # to few arguments
  expect_error(dsimplex(1, 0.8, 3, 4, 5)) # to many arguments
  expect_error(dsimplex(-1, mu = 0.8, sigma = 2)) # x is not allowed to be smaller 0
  # expect_error(dsimplex("r", mu = 0.8, sigma = 2)) # non-numeric arguments are disallowed
  expect_error(dsimplex(1, mu = 0, sigma = 2)) # mu is not allowed to be 0 or smaller
  expect_error(dsimplex(1, mu = 1, sigma = 2)) # mu is not allowed to be 1 or bigger
  expect_error(dsimplex(1, mu = 0.8, sigma = 0)) # p is not allowed to be 1 or smaller

  # Check rng for errors
  expect_error(rsimplex(100.82, 3, 4, 5)) # to many arguments
  expect_error(rsimplex(-1, mu = 0.8, sigma = 2)) # number of drawn samples cannot be smaller 0
  # expect_warning(expect_error(rsimplex("r", mu = 0.8, sigma = 2))) # non-numeric arguments are disallowed
  expect_error(rsimplex(100, mu = 0, sigma = 2)) # mu is not allowed to be 0 or smaller
  expect_error(rsimplex(100, mu = 1, sigma = 2)) # mu is not allowed to be 1 or bigger
  expect_error(rsimplex(100, mu = 0.8, sigma = 0)) # p is not allowed to be 0 or smaller

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = 0.5,
    aux_par = 2,
    ref_intercept = 0.5,
    rng_link = identity,
    parameter_link = logit,
    family = simplex,
    rng = rsimplex,
    aux_name = "sigma"
  )
})
