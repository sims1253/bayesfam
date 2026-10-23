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

  # recovery fit on seeded data, like the other family tests
  expect_brms_family(
    intercept = 2,
    aux_par = 2,
    aux2_par = 3,
    rng_link = exp,
    parameter_link = identity,
    family = shifted_inv_gaussian,
    rng = rshifted_inv_gaussian,
    aux_name = "shape",
    aux2_name = "ndt"
  )
})

test_that("dshifted_inv_gaussian is 0 below the shift", {
  expect_no_warning(d <- dshifted_inv_gaussian(c(0.5, 1, 2), 1, 2, 1))
  expect_equal(d, c(0, 0, brms::dinv_gaussian(1, 1, 2)))
  expect_equal(
    dshifted_inv_gaussian(c(0.5, 2), 1, 2, 1, log = TRUE),
    c(-Inf, brms::dinv_gaussian(1, 1, 2, log = TRUE))
  )
  # arguments of different lengths recycle per observation
  x <- c(0.5, 2)
  mu <- c(1, 2, 3)
  expect_equal(
    dshifted_inv_gaussian(x, mu, rep(2, 6), 1),
    sapply(1:6, function(i) {
      dshifted_inv_gaussian(rep_len(x, 6)[i], rep_len(mu, 6)[i], 2, 1)
    })
  )
})
