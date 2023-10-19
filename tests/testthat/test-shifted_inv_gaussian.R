test_that("shifted-inverse-gauss", {
  n <- 10000
  n_small <- 10
  eps <- 1e-3
  x <- seq(from = eps, to = 100, length.out = n)
  mu_list <- seq(from = 1, to = 100, length.out = n_small)
  shape_list <- seq(from = 1, to = 10, length.out = n_small)
  shift_list <- seq(from = 1, to = 10, length.out = n_small)
  accepted_relative_error <- 1e-6

  # Check lengths
  expect_equal(n, length(dshifted_inv_gaussian(x + 1, 0.5, 0.5, 0.5)))
  expect_equal(n, length(rshifted_inv_gaussian(n, 0.5, 0.5, 0.5)))

  for (mu in mu_list) {
    for (aux1 in shape_list) {
      for (aux2 in shift_list) {
        expect_eps(
          dshifted_inv_gaussian(x + aux2, mu, aux1, aux2),
          extraDistr::dwald(x, mu, aux1),
          eps = accepted_relative_error,
          relative = TRUE
        )
      }
    }
  }

  warning("RNG test missing, features of unbound dist branch required")

  # Check density function for errors
  expect_error(dshifted_inv_gaussian(0.5, 0.5, 0.5)) # to few arguments
  expect_error(dshifted_inv_gaussian(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)) # to many arguments
  expect_error(dshifted_inv_gaussian(-1, 0.5, 0.5, 0.5)) # x > 0
  expect_error(dshifted_inv_gaussian(0.5, -0.5, 0.5, 0.5)) # mu > 0
  expect_error(dshifted_inv_gaussian(0.5, 0.5, -0.5, 0.5)) # sigma > 0
  expect_error(dshifted_inv_gaussian(0.5, 0.5, 0.5, -0.5)) # shift >= 0
  # expect_true(is.na(dshifted_inv_gaussian(0.5, 0.5, 0.5, 0.5))) # shift results in x_unshifted <= 0
  # expect_warning(expect_error(dshifted_inv_gaussian("r", 0.5, 0.5, -0.5))) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(rshifted_inv_gaussian()) # to few arguments
  expect_error(rshifted_inv_gaussian(10, 0.8, 0.8, 0.8, 0.8)) # to many arguments
  expect_error(rshifted_inv_gaussian(-1, 0.5, 0.5, 0.5)) # number of drawn samples cannot be smaller 0
  # expect_error(rshifted_inv_gaussian("r", 0.5, 0.5, 0.5)) # non-numeric arguments are disallowed
  expect_error(rshifted_inv_gaussian(1, -0.5, 0.5, 0.5)) # mu > 0
  expect_error(rshifted_inv_gaussian(1, 0.5, -0.5, 0.5)) # sigma > 0
  expect_error(rshifted_inv_gaussian(1, 0.5, 0.5, -0.5)) # shift >= 0

  intercept <- 2
  shape <- 2
  shift <- 3
  thresh <- 0.05
  data <- list(y = rshifted_inv_gaussian(n = 1000, mu = exp(intercept), shape = shape, shift = shift))
  posterior_fit <- brms::brm(
    formula = y ~ 1, data = data,
    family = shifted_inv_gaussian(), stanvars = shifted_inv_gaussian()$stanvars,
    refresh = 0, silent = 2
  )

  intercept_recovered <- test_brms_quantile(
    posterior_fit, "b_Intercept", intercept, thresh
  )
  shape_par_recovered <- test_brms_quantile(
    posterior_fit, "shape", shape, thresh
  )
  shift_par_recovered <- test_brms_quantile(
    posterior_fit, "ndt", shift, thresh
  )
  expect_true(intercept_recovered && shape_par_recovered && shift_par_recovered)
})
