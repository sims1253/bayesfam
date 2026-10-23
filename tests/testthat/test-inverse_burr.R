mu_list <- c(0.5, 5)
tau_list <- c(0.5, 4)
gamma_list <- c(1.5, 6)

test_that("inverse_burr density matches actuar", {
  skip_if_not_installed("actuar")
  x <- exp(seq(from = log(1e-3), to = log(100), length.out = 200))
  for (mu in mu_list) {
    for (tau in tau_list) {
      for (gamma in gamma_list) {
        expect_eps(
          dinverse_burr(x, mu = mu, tau = tau, gamma = gamma),
          actuar::dinvburr(
            x,
            shape1 = tau,
            shape2 = gamma,
            scale = mu * (2^(1 / tau) - 1)^(1 / gamma)
          ),
          eps = 1e-10,
          relative = TRUE
        )
      }
    }
  }
})

test_that("custom-inverse_burr", {
  unit <- seq(from = 0.01, to = 0.99, length.out = 99)
  for (mu in mu_list) {
    for (tau in tau_list) {
      for (gamma in gamma_list) {
        # Density integrates to 1, mu is the median, CDF/quantile round trip
        expect_eps(
          integrate(
            dinverse_burr,
            0,
            Inf,
            mu = mu,
            tau = tau,
            gamma = gamma
          )$value,
          1,
          eps = 1e-6
        )
        expect_eps(pinverse_burr(mu, mu, tau, gamma), 0.5, eps = 1e-12)
        expect_eps(
          pinverse_burr(qinverse_burr(unit, mu, tau, gamma), mu, tau, gamma),
          unit,
          eps = 1e-10,
          relative = TRUE
        )
        # RNG recovers the median
        expect_eps(
          median(rinverse_burr(1e5, mu, tau, gamma)),
          mu,
          eps = 0.03,
          relative = TRUE
        )
      }
    }
  }

  # Check argument checks
  expect_error(dinverse_burr(-1, mu = 1, tau = 1, gamma = 1))
  expect_error(dinverse_burr(1, mu = 1, tau = 1, gamma = 0))
  expect_error(qinverse_burr(c(0.5, 1), mu = 1, tau = 1, gamma = 1))
  expect_error(rinverse_burr(10, mu = 0, tau = 1, gamma = 1))

  # brms callbacks with a mock prep; epred is Inf where gamma <= 1
  ndraws <- 4
  y <- c(0.5, 1, 3)
  prep <- structure(
    list(
      data = list(Y = y),
      ndraws = ndraws,
      dpars = list(
        mu = matrix(2, ndraws, length(y)),
        tau = matrix(2, ndraws, length(y)),
        gamma = matrix(c(3, 3, 1, 0.5), ndraws, length(y))
      )
    ),
    class = "brmsprep"
  )
  expect_equal(
    log_lik_inverse_burr(2, prep),
    dinverse_burr(1, 2, 2, c(3, 3, 1, 0.5), log = TRUE)
  )
  expect_length(posterior_predict_inverse_burr(2, prep), ndraws)
  epred <- posterior_epred_inverse_burr(prep)
  expect_equal(dim(epred), c(ndraws, length(y)))
  expect_true(all(is.finite(epred[1:2, ])))
  expect_true(all(epred[3:4, ] == Inf))

  # The injected Stan code parses and brms accepts the gamma dpar
  expect_no_error(brms::make_stancode(
    y ~ 1,
    data = list(y = y),
    family = inverse_burr(),
    stanvars = inverse_burr()$stanvars,
    parse = TRUE
  ))
})

test_that("inverse_burr brms fit", {
  # Check of brms can fit the custom family and recover the intercept and shapes
  expect_brms_family(
    intercept = 5,
    aux_par = 2,
    aux2_par = 3,
    ref_intercept = 5,
    rng_link = identity,
    parameter_link = log,
    family = inverse_burr,
    rng = rinverse_burr,
    aux_name = "tau",
    aux2_name = "gamma"
  )
})
