test_that("lognormal", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
  n_small <- 10
  mu_list <- seq(from = eps, to = 1 - eps, length.out = n_small)
  sigma_list <- seq(from = 0.01, to = 10, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.1
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(dlognormal(x, mu = log(1), sigma = 2)))
  expect_equal(n, length(rlognormal(n, mu = log(1), sigma = 2)))

  # Compare density and to built-in
  for (mu in mu_list) {
    for (sigma in sigma_list) {
      expect_eps(dlognormal(x, mu = log(mu), sigma = sigma),
        dlnorm(x, log(mu), sigma),
        eps = accepted_relative_error,
        relative = TRUE
      )
    }
  }

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rlognormal,
    metric_mu = median,
    n = 10 * n,
    mu_list = mu_list,
    aux_list = sigma_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE,
    mu_link = log
  )

  # Check if the RNG can recover the quantiles
  test_rng_quantiles(
    rng_fun = rlognormal,
    quantile_fun = qlnorm,
    n = 20 * n,
    mu_list = mu_list,
    aux_list = sigma_list,
    eps = accepted_rng_error,
    quantiles = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE,
    mu_link = log
  )

  # Check density function for errors
  expect_error(dlognormal(0.5, 2)) # to few arguments
  expect_error(dlognormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(dlognormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(dlognormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_warning(expect_error(dlognormal("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(rlognormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rlognormal(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_error(rlognormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed
  expect_error(rlognormal(100, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
})
