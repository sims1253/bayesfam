test_that("custom-cloglognormal", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = eps, to = 1 - eps, length.out = n_small)
  sigma_list <- seq(from = 0.01, to = 10, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.05
  accepred_rng_failures <- 0.05

  # Check lengths
  expect_equal(n, length(dcloglognormal(x, mu = cloglog(0.5), sigma = 0.4)))
  expect_equal(n, length(rcloglognormal(n, mu = cloglog(0.5), sigma = 0.4)))

  # Compare density to reference implementation
  warning("No reference density available to test against!")

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rcloglognormal,
    metric_mu = median,
    n = 10 * n,
    mu_list = mu_list,
    aux_list = sigma_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE,
    mu_link = cloglog
  )

  # Check if the RNG can recover the quantiles
  warning("No quantile function available to test rng quantile recovery.")

  # Check density function for errors
  expect_error(dcloglognormal(0.5, 2)) # to few arguments
  expect_error(dcloglognormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(dcloglognormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(dcloglognormal(2, mu = 2, sigma = 2)) # x is not allowed to be bigger 1
  expect_error(dcloglognormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  # expect_error(dcloglognormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(rcloglognormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rcloglognormal(-1, mu = cloglog(0.5), sigma = 0.4)) # number of drawn samples cannot be smaller 0
  # expect_warning(expect_error(rcloglognormal("r", mu = cloglog(0.5), sigma = 0.4))) # non-numeric arguments are disallowed
  expect_error(rcloglognormal(100, mu = cloglog(-1), sigma = 0.4)) # mu must be between 0 and 1
  expect_error(rcloglognormal(100, mu = cloglog(0.5), sigma = -1)) # sigma is not allowed to be negative

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = cloglog(0.5),
    aux_par = 0.4,
    ref_intercept = 0.5,
    parameter_link = cloglog,
    rng_link = identity,
    family = cloglognormal,
    rng = rcloglognormal,
    aux_name = "sigma"
  )
})
