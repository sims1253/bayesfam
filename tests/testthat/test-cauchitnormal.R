test_that("custom-cauchitnormal", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-3
  x <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = 0.1, to = 0.9, length.out = 9)
  sigma_list <- seq(from = 0.1, to = 10, length.out = n_small)
  accepted_rng_error <- 0.01
  accepred_rng_failures <- 0.11

  # Check lengths
  expect_equal(n, length(dcauchitnormal(x, mu = cauchit(0.5), sigma = 0.4)))
  expect_equal(n, length(rcauchitnormal(n, mu = cauchit(0.5), sigma = 0.4)))

  # Compare density to reference implementation
  warning("No reference density available to test against!")

  # check if the RNG is close enough to the true mean in most cases
  bayesfam:::test_rng(
    rng_fun = rcauchitnormal,
    metric_mu = median,
    n = 50 * n,
    mu_list = mu_list,
    aux_list = sigma_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE,
    mu_link = cauchit
  )

  # Check if the RNG can recover the quantiles
  warning("No quantile function available to test rng quantile recovery.")

  # Check density function for errors
  expect_error(dcauchitnormal(0.5, 2)) # to few arguments
  expect_error(dcauchitnormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(dcauchitnormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(dcauchitnormal(2, mu = 2, sigma = 2)) # x is not allowed to be bigger 1
  expect_error(dcauchitnormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  # expect_error(dcauchitnormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(rcauchitnormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rcauchitnormal(-1, mu = cauchit(0.5), sigma = 0.4)) # number of drawn samples cannot be smaller 0
  # expect_warning(expect_error(rcauchitnormal("r", mu = cauchit(0.5), sigma = 0.4))) # non-numeric arguments are disallowed
  expect_error(rcauchitnormal(100, mu = cauchit(-1), sigma = 0.4)) # mu must be between 0 and 1
  expect_error(rcauchitnormal(100, mu = cauchit(0.5), sigma = -1)) # sigma is not allowed to be negative

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = cauchit(0.5),
    aux_par = 0.4,
    ref_intercept = 0.5,
    rng_link = identity,
    parameter_link = cauchit,
    family = cauchitnormal,
    rng = rcauchitnormal,
    aux_name = "sigma"
  )
})
