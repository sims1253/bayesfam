test_that("custom-softplusnormal", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
  n_small <- 10
  mu_list <- seq(from = eps, to = 1 - eps, length.out = n_small)
  sigma_list <- seq(from = 0.01, to = 10, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.03
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(dsoftplusnormal(x, mu = softplus(1), sigma = 2)))
  expect_equal(n, length(rsoftplusnormal(n, mu = softplus(1), sigma = 2)))


  # Compare density and to built-in
  warning("No reference density available to test against.")

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rsoftplusnormal,
    metric_mu = median,
    n = 10 * n,
    mu_list = mu_list,
    aux_list = sigma_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE,
    mu_link = softplus
  )

  # Check if the RNG can recover the quantiles
  warning("No quantile function available to test rng quantile recovery.")

  # Check density function for errors
  expect_error(dsoftplusnormal(0.5, 2)) # to few arguments
  expect_error(dsoftplusnormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(dsoftplusnormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(dsoftplusnormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_error(dsoftplusnormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(rsoftplusnormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rsoftplusnormal(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(rsoftplusnormal("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed
  expect_error(rsoftplusnormal(100, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = softplus(5),
    aux_par = 2,
    ref_intercept = 5,
    parameter_link = softplus,
    rng_link = identity,
    family = softplusnormal,
    rng = rsoftplusnormal,
    aux_name = "sigma"
  )
})
