test_that("lognormal_natural", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
  n_small <- 10
  mu_list <- seq(from = eps, to = 10, length.out = n_small)
  sigma_list <- seq(from = 0.01, to = 10, length.out = n_small)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.1
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(dlognormal_natural(x, mu = 1, sigma = 2)))
  expect_equal(n, length(rlognormal_natural(n, mu = 1, sigma = 2)))

  # Compare density and to built-in
  for (mu in mu_list) {
    for (sigma in sigma_list) {
      common_term <- log(1+sigma^2/mu^2)
      expect_eps(dlognormal_natural(x, mu = mu, sigma = sigma),
                 dlnorm(x, log(mu)-common_term/2, sqrt(common_term)),
                 eps = accepted_relative_error,
                 relative = TRUE
      )
    }
  }

  # check if the RNG is close enough to the true mean in most cases
  test_rng(
    rng_fun = rlognormal_natural,
    metric_mu = mean,
    n = 10 * n,
    mu_list = mu_list,
    aux_list = sigma_list,
    mu_eps = accepted_rng_error,
    p_acceptable_failures = accepred_rng_failures,
    relative = TRUE,
    mu_link = identity
  )

  # Check density function for errors
  expect_error(dlognormal_natural(0.5, 2)) # to few arguments
  expect_error(dlognormal_natural(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(dlognormal_natural(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(dlognormal_natural(0.5, mu = -1, sigma = 1)) # mu has to be > 0
  expect_warning(expect_error(dlognormal_natural("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(rlognormal_natural(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rlognormal_natural(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_error(rlognormal_natural("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed
  expect_error(rlognormal_natural(0.5, mu = -1, sigma = 1)) # mu has to be > 0

  expect_brms_family(
    intercept = 5,
    aux_par = 2,
    ref_intercept = 5,
    rng_link = exp,
    parameter_link = identity,
    family = lognormal_natural,
    rng = rlognormal_natural,
    aux_name = "sigma"
  )
})
