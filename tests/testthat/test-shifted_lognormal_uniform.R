test_that("custom-shifted_lognormal_uniform", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
  unit <- seq(from = eps, to = 1 - eps, length.out = n)
  accepted_relative_error <- 1e-6
  accepted_rng_error <- 0.085
  accepred_rng_failures <- 0.1

  # Check lengths
  expect_equal(n, length(
    dshifted_lognormal_uniform(x, meanlog = 1, sdlog = 2,mix = 0.1,
                               max_uniform = 10, shift = 0.9)))
  expect_equal(n, length(
    rshifted_lognormal_uniform(n, meanlog = 1, sdlog = 2,
                               mix = 0.1, max_uniform = 10,shift = 0.9)))

  # Check density function for errors
  expect_error(dshifted_lognormal_uniform(-1, meanlog = 2)) # y is not allowed to be smaller 0
  expect_error(dshifted_lognormal_uniform(1, sdlog = 0)) # sdlog is not allowed to be 0 or smaller
  expect_error(dshifted_lognormal_uniform(1, meanlog = 1, shift = -0.1))


  # Check rng for errors
  expect_error(rshifted_lognormal_uniform(-1)) # number of drawn samples cannot be smaller 0
  expect_error(dshifted_lognormal_uniform(100, sdlog = 0)) # mu is not allowed to be 0 or smaller
  expect_error(dshifted_lognormal_uniform(100, meanlog = 1, shift = -0.1))

  # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = 5,
    aux_par = 2,
    ref_intercept = 5,
    rng_link = identity,
    parameter_link = identity,
    family = shifted_lognormal_uniform,
    rng = rshifted_lognormal_uniform,
    aux_name = "sigma",
    formula = y | vreal(0.01, 100) ~ 1,
    prior = c(brms::prior(constant(0.1), class = "mix"),
              brms::prior(constant(0), class = "shiftprop"))

  )
})
