test_that("custom-logitnormal", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = eps, to = 1 - eps, length.out = n_small)
  sigma_list <- seq(from = 0.01, to = 10, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.1
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(dlogitnormal(x, mu = logit(0.5), sigma = 0.4)))
  expect_equal(n, length(rlogitnormal(n, mu = logit(0.5), sigma = 0.4)))

  # Compare density to reference implementation
  warning("No reference density available to test against!")

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rlogitnormal,
    metric_mu = median,
    n = 10 * n,
    mu_list = mu_list,
    aux_list = sigma_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE,
    mu_link = logit
  )

  # Check if the RNG can recover the quantiles
  warning("No quantile function available to test rng quantile recovery.")

  # Check density function for errors
  expect_error(dlogitnormal(0.5, 2)) # to few arguments
  expect_error(dlogitnormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(dlogitnormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(dlogitnormal(2, mu = 2, sigma = 2)) # x is not allowed to be bigger 1
  expect_error(dlogitnormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_error(dlogitnormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed


  # Check rng for errors
  expect_error(rlogitnormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rlogitnormal(-1, mu = logit(0.5), sigma = 0.4)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(rlogitnormal("r", mu = logit(0.5), sigma = 0.4))) # non-numeric arguments are disallowed
  expect_error(rlogitnormal(100, mu = logit(-1), sigma = 0.4)) # mu must be between 0 and 1
  expect_error(rlogitnormal(100, mu = logit(0.5), sigma = -1)) # sigma is not allowed to be negative

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = logit(0.5),
    aux_par = 0.4,
    ref_intercept = 0.5,
    parameter_link = logit,
    rng_link = identity,
    family = logitnormal,
    rng = rlogitnormal,
    aux_name = "sigma"
  )
})
