test_that("shifted-inverse-gauss", {

  n <- 10
  x <- seq(from=0.1, to=0.9, length.out=n)
  # Check lengths
  expect_equal(n, length(dshifted_inv_gaussian(x, 0.5, 0.5, -0.5)))
  expect_equal(n, length(rshifted_inv_gaussian(n, 0.5, 0.5, 0.5)))

  # Check density function for errors
  expect_error(dshifted_inv_gaussian(0.5, 0.5, 0.5)) # to few arguments
  expect_error(dshifted_inv_gaussian(0.5, 0.5, 0.5, -0.5, 0.5)) # to many arguments
  expect_error(dshifted_inv_gaussian(-1, 0.5, 0.5, -0.5)) # x > 0
  expect_error(dshifted_inv_gaussian(0.5, -0.5, 0.5, -0.5)) # mu > 0
  expect_error(dshifted_inv_gaussian(0.5, 0.5, -0.5, -0.5)) # sigma > 0
  expect_true(is.na(dshifted_inv_gaussian(0.5, 0.5, 0.5, 0.5))) # shift results in x_unshifted <= 0
  expect_warning(expect_error(dshifted_inv_gaussian("r", 0.5, 0.5, -0.5))) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(rshifted_inv_gaussian()) # to few arguments
  expect_error(rshifted_inv_gaussian(10, 0.8, 0.8, 0.8, 0.8)) # to many arguments
  expect_error(rshifted_inv_gaussian(-1, 0.5, 0.5, -0.5)) # number of drawn samples cannot be smaller 0
  expect_error(rshifted_inv_gaussian("r", 0.5, 0.5, -0.5)) # non-numeric arguments are disallowed
  expect_error(rshifted_inv_gaussian(1, -0.5, 0.5, -0.5)) # mu > 0
  expect_error(rshifted_inv_gaussian(1, 0.5, -0.5, -0.5)) # sigma > 0

  data <- list(y = rshifted_inv_gaussian(n = 1000, mu = exp(1), shape = exp(3), shift = -1))
  fit <- brms::brm(formula = y ~ 1, data = data,
                   family = shifted_inv_gaussian(), stanvars = shifted_inv_gaussian()$stanvars,
                   refresh = 0)
})
