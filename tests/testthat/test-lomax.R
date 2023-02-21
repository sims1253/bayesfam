test_that("custom-lomax", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
  unit <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = eps, to = 1000, length.out = n_small)
  alpha_list <- seq(from = 1 + eps, to = 50, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.03
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(dlomax(x, mu = 5, alpha = 2)))
  expect_equal(n, length(qlomax(unit, mu = 5, alpha = 2)))
  expect_equal(n, length(rlomax(n, 5, 2)))

  # Compare density and quantile functions to extraDistr
  for (mu in mu_list) {
    for (alpha in alpha_list) {
      expect_eps(
        dlomax(x, mu = mu, alpha = alpha),
        extraDistr::dlomax(x, lambda = 1 / (mu * (alpha - 1)), kappa = alpha),
        eps = accepted_relative_error,
        relative = TRUE
      )
      expect_eps(
        qlomax(unit, mu = mu, alpha = alpha),
        extraDistr::qlomax(unit, lambda = 1 / (mu * (alpha - 1)), kappa = alpha),
        eps = accepted_relative_error,
        relative = TRUE
      )
    }
  }

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = 5 * n,
    mu_list = mu_list,
    aux_list = alpha_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check if the RNG can recover the quantiles
  test_rng_quantiles(
    rng_fun = rlomax,
    quantile_fun = qlomax,
    n = 20 * n,
    mu_list = mu_list,
    aux_list = alpha_list,
    eps = accepted_rng_error,
    quantiles = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check density function for errors
  expect_error(dlomax(1, 2)) # to few arguments
  expect_error(dlomax(1, 2, 3, 4, 5)) # to many arguments
  expect_error(dlomax(-1, mu = 2, alpha = 2)) # x is not allowed to be smaller 0
  expect_error(dlomax("r", mu = 2, alpha = 2)) # non-numeric arguments are disallowed
  expect_error(dlomax(1, mu = 0, alpha = 2)) # mu is not allowed to be 0 or smaller
  expect_error(dlomax(1, mu = 1, alpha = 0)) # alpha is not allowed to be 1 or smaller

  # Check quantile function for errors
  expect_error(qlomax(1, 2)) # to few arguments
  expect_error(qlomax(1, 2, 3, 4, 5)) # to many arguments
  expect_error(qlomax("r", mu = 2, alpha = 2)) # non-numeric arguments are disallowed
  expect_error(qlomax(1, mu = 0, alpha = 2)) # mu is not allowed to be 0 or smaller
  expect_error(qlomax(1, mu = 1, alpha = 0)) # alpha is not allowed to be 0 or smaller
  expect_error(qlomax(c(-1, 2), mu = 2, alpha = 2)) # q is not allowed to be outside [0, 1]

  # Check rng for errors
  expect_error(rlomax(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rlomax(-1, mu = 2, alpha = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(rlomax("r", mu = 2, alpha = 2))) # non-numeric arguments are disallowed
  expect_error(rlomax(100, mu = 0, alpha = 2)) # mu is not allowed to be 0 or smaller
  expect_error(rlomax(100, mu = 1, alpha = 0)) # alpha is not allowed to be 0 or smaller

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    n_data_sampels = 10000,
    intercept = 5,
    aux_par = 2,
    ref_intercept = 5,
    rng_link = identity,
    parameter_link = log,
    family = lomax,
    rng = rlomax,
    aux_name = "alpha"
  )
})
