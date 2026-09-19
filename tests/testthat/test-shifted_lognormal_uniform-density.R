test_that("dshifted_lognormal_uniform integrates to 1 (#24)", {
  integrate_density <- function(meanlog, sdlog, mix, shift, max_uniform) {
    # split the domain at shift and max_uniform where components switch on/off
    pieces <- integrate(
      function(y) exp(dshifted_lognormal_uniform(
        y, meanlog, sdlog, mix, shift, max_uniform
      )),
      lower = 0, upper = shift
    )$value +
      integrate(
        function(y) exp(dshifted_lognormal_uniform(
          y, meanlog, sdlog, mix, shift, max_uniform
        )),
        lower = shift, upper = max_uniform
      )$value +
      integrate(
        function(y) exp(dshifted_lognormal_uniform(
          y, meanlog, sdlog, mix, shift, max_uniform
        )),
        lower = max_uniform, upper = Inf
      )$value
    pieces
  }

  for (mix in c(0, 0.1, 0.5, 1)) {
    expect_eps(
      integrate_density(0, 1, mix, shift = 0.5, max_uniform = 2),
      1,
      eps = 1e-4
    )
  }
  # also with a nonzero meanlog / larger sdlog
  for (mix in c(0, 0.25, 1)) {
    expect_eps(
      integrate_density(0.7, 1.5, mix, shift = 0.3, max_uniform = 5),
      1,
      eps = 1e-4
    )
  }
})

test_that("dshifted_lognormal_uniform matches the Stan lpdf formula (#24)", {
  meanlog <- 0
  sdlog <- 1
  mix <- 0.1
  shift <- 0.5
  max_uniform <- 2

  # mirror of the shifted_lognormal_uniform_lpdf Stan function
  stan_lpdf <- function(y, meanlog, sdlog, mix, shift, max_uniform) {
    if (y <= shift) {
      log(mix) + stats::dunif(y, 0, max_uniform, log = TRUE)
    } else if (y >= max_uniform) {
      log1p(-mix) + stats::dlnorm(y - shift, meanlog, sdlog, log = TRUE)
    } else {
      uniform_llh <- dunif(y, 0, max_uniform, log = TRUE)
      lognormal_llh <- dlnorm(y - shift, meanlog, sdlog, log = TRUE)
      log(mix * exp(uniform_llh) + (1 - mix) * exp(lognormal_llh))
    }
  }

  # below the shift: only the uniform component
  expect_eps(
    dshifted_lognormal_uniform(0.2, meanlog, sdlog, mix, shift, max_uniform),
    stan_lpdf(0.2, meanlog, sdlog, mix, shift, max_uniform),
    eps = 1e-12
  )
  # in the overlap: both components
  expect_eps(
    dshifted_lognormal_uniform(1, meanlog, sdlog, mix, shift, max_uniform),
    stan_lpdf(1, meanlog, sdlog, mix, shift, max_uniform),
    eps = 1e-12
  )
  # above the uniform bound: only the lognormal component
  expect_eps(
    dshifted_lognormal_uniform(3, meanlog, sdlog, mix, shift, max_uniform),
    stan_lpdf(3, meanlog, sdlog, mix, shift, max_uniform),
    eps = 1e-12
  )

  # the example from the issue: no truncation normalizer anymore
  expect_eps(
    dshifted_lognormal_uniform(2, 0, 1, mix = 0.1, shift = 0, max_uniform = 1),
    log1p(-0.1) + dlnorm(2, 0, 1, log = TRUE),
    eps = 1e-12
  )

  # edge mixture weights keep working
  expect_eps(
    dshifted_lognormal_uniform(1, meanlog, sdlog, mix = 0, shift, max_uniform),
    stan_lpdf(1, meanlog, sdlog, 0, shift, max_uniform),
    eps = 1e-12
  )
  expect_eps(
    dshifted_lognormal_uniform(0.2, meanlog, sdlog, mix = 1, shift, max_uniform),
    stan_lpdf(0.2, meanlog, sdlog, 1, shift, max_uniform),
    eps = 1e-12
  )
  # all-uniform mixture above the bound has zero density
  expect_equal(
    dshifted_lognormal_uniform(3, meanlog, sdlog, mix = 1, shift, max_uniform),
    -Inf
  )
  # all-lognormal mixture below the shift has zero density
  expect_equal(
    dshifted_lognormal_uniform(0.2, meanlog, sdlog, mix = 0, shift, max_uniform),
    -Inf
  )

  # vectorized: scalar y with posterior-draw vectors (log_lik pattern)
  S <- 5
  expect_equal(
    length(dshifted_lognormal_uniform(
      1.3,
      meanlog = rnorm(S, 0, 0.1),
      sdlog = rep(1, S),
      mix = rep(0.1, S),
      shift = rep(0.5, S),
      max_uniform = rep(2, S)
    )),
    S
  )
  # vector y with scalar parameters
  expect_equal(
    length(dshifted_lognormal_uniform(
      c(0.2, 1, 3),
      meanlog = 0,
      sdlog = 1,
      mix = 0.1,
      shift = 0.5,
      max_uniform = 2
    )),
    3
  )
})
