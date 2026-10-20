test_that("dshifted_lognormal_uniform returns the density unless log = TRUE", {
  y <- c(0.5, 2, 20)
  d <- dshifted_lognormal_uniform(y, meanlog = 1, sdlog = 0.5, shift = 1)
  expect_equal(
    dshifted_lognormal_uniform(
      y,
      meanlog = 1,
      sdlog = 0.5,
      shift = 1,
      log = TRUE
    ),
    log(d)
  )
  expect_true(all(d > 0 & d < 1))
})

test_that("custom-shifted_lognormal_uniform", {
  # Setup of testing space
  n <- 10000
  eps <- 1e-6
  x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
  unit <- seq(from = eps, to = 1 - eps, length.out = n)

  # Check lengths
  expect_equal(
    n,
    length(
      dshifted_lognormal_uniform(
        x,
        meanlog = 1,
        sdlog = 2,
        mix = 0.1,
        max_uniform = 10,
        shift = 0.9
      )
    )
  )
  expect_equal(
    n,
    length(
      rshifted_lognormal_uniform(
        n,
        meanlog = 1,
        sdlog = 2,
        mix = 0.1,
        max_uniform = 10,
        shift = 0.9
      )
    )
  )

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
    prior = c(
      brms::prior(constant(0.1), class = "mix"),
      brms::prior(constant(0), class = "shiftprop")
    )
  )
})

test_that("dshifted_lognormal_uniform integrates to 1 (#24)", {
  shift <- 0.5
  max_uniform <- 2
  dens <- function(y, mix) {
    dshifted_lognormal_uniform(y, 0, 1, mix, shift, max_uniform)
  }
  # integrate piecewise over the support of each mixture
  for (mix in c(0, 0.5, 1)) {
    breaks <- c(
      if (mix > 0) 0,
      shift,
      max_uniform,
      if (mix < 1) Inf
    )
    total <- sum(vapply(
      seq_len(length(breaks) - 1),
      function(k) integrate(dens, breaks[k], breaks[k + 1], mix = mix)$value,
      numeric(1)
    ))
    expect_equal(total, 1, tolerance = 1e-4)
  }
})

test_that("posterior_epred_shifted_lognormal_uniform uses per-observation bounds (#25)", {
  max_shift <- c(1, 2, 4)
  max_uniform <- c(4, 8, 16)
  for (S in c(2, 3)) {
    # length-S dpar vectors (constant over observations), N = 3
    mu <- seq_len(S) * 0.1
    prep <- structure(
      list(
        ndraws = S,
        dpars = list(
          mu = mu,
          sigma = rep(0.7, S),
          mix = rep(0.1, S),
          shiftprop = rep(0.5, S)
        ),
        data = list(vreal1 = max_shift, vreal2 = max_uniform)
      ),
      class = "brmsprep"
    )
    expected <- t(vapply(
      mu,
      function(m) {
        0.1 * 0.5 * max_uniform + 0.9 * (0.5 * max_shift + exp(m + 0.7^2 / 2))
      },
      numeric(3)
    ))
    expect_equal(posterior_epred_shifted_lognormal_uniform(prep), expected)
  }
})
