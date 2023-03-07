test_that("custom-gumbel", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- exp(seq(from = -100, to = 100, length.out = n))
  unit <- seq(from = eps, to = 1 - eps, length.out = n)
  n_small <- 10
  mu_list <- seq(from = -50, to = 50, length.out = n_small)
  aux_list <- seq(from = eps, to = 50, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.085
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(dgumbel_mean(x, mu = 1, sigma = 2)))
  expect_equal(n, length(qgumbel_mean(unit, mu = 1, sigma = 2)))
  expect_equal(n, length(rgumbel_mean(n, mu = 1, sigma = 2)))

  euler_mascheroni <- 0.57721566490153

  # Compare density and quantile functions to extraDistr
  for (mu in mu_list) {
    for (aux in aux_list) {
      loc <- mu - aux * euler_mascheroni
      expect_eps(
        dgumbel_mean(x, mu = mu, sigma = aux),
        evd::dgumbel(x, loc, aux),
        eps = accepted_relative_error,
        relative = TRUE
      )
      expect_eps(
        qgumbel_mean(unit, mu = mu, sigma = aux),
        evd::qgumbel(unit, loc, aux),
        eps = accepted_relative_error,
        relative = TRUE
      )
    }
  }

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rgumbel_mean,
    metric_mu = mean,
    n = 5 * n,
    mu_list = mu_list,
    aux_list = aux_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check if the RNG can recover the quantiles
  test_rng_quantiles(
    rng_fun = rgumbel_mean,
    quantile_fun = qgumbel_mean,
    n = 5 * n,
    mu_list = mu_list,
    aux_list = aux_list,
    eps = accepted_rng_error,
    quantiles = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE
  )

  # Check density function for errors
  expect_error(dgumbel_mean(1, 2)) # to few arguments
  expect_error(dgumbel_mean(1, 2, 3, 4, 5)) # to many arguments
  expect_error(dgumbel_mean(1, mu = 1, sigma = 0)) # aux is not allowed to be 0 or smaller
  expect_error(dgumbel_mean("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed

  # Check quantile function for errors
  expect_error(qgumbel_mean(1, 2)) # to few arguments
  expect_error(qgumbel_mean(1, 2, 3, 4, 5)) # to many arguments
  expect_error(qgumbel_mean(1, mu = 1, sigma = 0)) # aux is not allowed to be 0 or smaller
  expect_error(qgumbel_mean(c(-1, 2), mu = 2, sigma = 2)) # q is not allowed to be outside [0, 1]
  expect_error(qgumbel_mean("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(rgumbel_mean(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rgumbel_mean(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_error(rgumbel_mean("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed
  expect_error(rgumbel_mean(100, mu = 1, sigma = 0)) # phi is not allowed to be 0 or smaller

  #skip("not implemented")

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = 5,
    aux_par = 2,
    ref_intercept = 5,
    rng_link = identity,
    parameter_link = identity,
    family = dgumbel_mean,
    rng = rgumbel_mean,
    aux_name = "sigma"
  )
})
