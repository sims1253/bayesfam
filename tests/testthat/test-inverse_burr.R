mock_inverse_burr_prep <- function(y, mu, tau, gamma, ndraws = 5) {
  n <- length(y)
  structure(
    list(
      data = list(Y = y),
      ndraws = ndraws,
      family = inverse_burr(),
      dpars = list(
        mu = matrix(mu, nrow = ndraws, ncol = n),
        tau = matrix(tau, nrow = ndraws, ncol = n),
        gamma = matrix(gamma, nrow = ndraws, ncol = n)
      )
    ),
    class = "brmsprep"
  )
}

inverse_burr_mean <- function(mu, tau, gamma) {
  scale <- mu * (2^(1 / tau) - 1)^(1 / gamma)
  scale * exp(
    lgamma(tau + 1 / gamma) +
      lgamma(1 - 1 / gamma) -
      lgamma(tau)
  )
}

test_that("inverse_burr density matches actuar", {
  skip_if_not_installed("actuar")
  x <- exp(seq(from = log(1e-3), to = log(100), length.out = 200))
  for (mu in c(0.5, 2, 5)) {
    for (tau in c(0.5, 2, 4)) {
      for (gamma in c(1.5, 3, 6)) {
        scale <- mu * (2^(1 / tau) - 1)^(1 / gamma)
        expect_eps(
          dinverse_burr(x, mu = mu, tau = tau, gamma = gamma),
          actuar::dinvburr(
            x,
            shape1 = tau,
            shape2 = gamma,
            scale = scale
          ),
          eps = 1e-12,
          relative = TRUE
        )
        expect_eps(
          pinverse_burr(x, mu = mu, tau = tau, gamma = gamma),
          actuar::pinvburr(
            x,
            shape1 = tau,
            shape2 = gamma,
            scale = scale
          ),
          eps = 1e-12,
          relative = TRUE
        )
        p <- seq(from = 0.01, to = 0.99, length.out = 100)
        expect_eps(
          qinverse_burr(p, mu = mu, tau = tau, gamma = gamma),
          actuar::qinvburr(
            p,
            shape1 = tau,
            shape2 = gamma,
            scale = scale
          ),
          eps = 1e-12,
          relative = TRUE
        )
      }
    }
  }
})

test_that("inverse_burr density integrates to one and mu is the median", {
  for (mu in c(0.5, 2)) {
    for (tau in c(0.5, 2)) {
      for (gamma in c(1.5, 3)) {
        expect_eps(
          stats::integrate(
            dinverse_burr,
            lower = 0,
            upper = Inf,
            mu = mu,
            tau = tau,
            gamma = gamma
          )$value,
          1,
          eps = 1e-6
        )
        # reference CDF at mu equals 0.5
        expect_eps(
          pinverse_burr(mu, mu = mu, tau = tau, gamma = gamma),
          0.5,
          eps = 1e-12
        )
        expect_eps(
          stats::integrate(
            dinverse_burr,
            lower = 0,
            upper = mu,
            mu = mu,
            tau = tau,
            gamma = gamma
          )$value,
          0.5,
          eps = 1e-6
        )
      }
    }
  }
})

test_that("inverse_burr CDF/quantile round trip", {
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  for (tau in c(0.5, 1, 4)) {
    for (gamma in c(1.5, 3)) {
      expect_eps(
        pinverse_burr(
          qinverse_burr(p, mu = 2, tau = tau, gamma = gamma),
          mu = 2,
          tau = tau,
          gamma = gamma
        ),
        p,
        eps = 1e-12,
        relative = TRUE
      )
    }
  }
})

test_that("inverse_burr RNG recovers median, mean and quantiles", {
  set.seed(4826)
  n <- 2e5
  for (mu in c(0.5, 2)) {
    for (tau in c(0.5, 2, 4)) {
      for (gamma in c(1.5, 3)) {
        draws <- rinverse_burr(n, mu = mu, tau = tau, gamma = gamma)
        # mu is the median for any shape parameters
        expect_eps(median(draws), mu, eps = 0.02, relative = TRUE)
        # for gamma > 1 the mean exists and matches the closed form
        expect_eps(
          mean(draws),
          inverse_burr_mean(mu, tau, gamma),
          eps = 0.05,
          relative = TRUE
        )
      }
    }
  }
  mu <- 2
  tau <- 2
  gamma <- 3
  draws <- rinverse_burr(n, mu = mu, tau = tau, gamma = gamma)
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  expect_eps(
    qinverse_burr(p, mu = mu, tau = tau, gamma = gamma),
    stats::quantile(draws, probs = p),
    eps = 0.03,
    r = 0.2,
    relative = TRUE
  )
})

test_that("inverse_burr argument checks", {
  expect_error(dinverse_burr(-1, mu = 1, tau = 1, gamma = 1))
  expect_error(dinverse_burr(1, mu = 0, tau = 1, gamma = 1))
  expect_error(dinverse_burr(1, mu = 1, tau = 0, gamma = 1))
  expect_error(dinverse_burr(1, mu = 1, tau = 1, gamma = 0))
  expect_error(pinverse_burr(-1, mu = 1, tau = 1, gamma = 1))
  expect_error(qinverse_burr(c(-0.1, 0.5), mu = 1, tau = 1, gamma = 1))
  expect_error(qinverse_burr(c(0.5, 1), mu = 1, tau = 1, gamma = 1))
  expect_error(qinverse_burr(0.5, mu = 1, tau = 0, gamma = 1))
  expect_error(rinverse_burr(10, mu = 1, tau = 1, gamma = 0))
})

test_that("inverse_burr brms callbacks work with a mock prep", {
  y <- c(0.4, 1.1, 2.9)
  mu <- 1.5
  tau <- 2
  gamma <- 3
  prep <- mock_inverse_burr_prep(y, mu, tau, gamma, ndraws = 5)

  # log_lik returns one log density per draw for observation i
  ll <- log_lik_inverse_burr(2, prep)
  expect_length(ll, 5)
  expect_eps(
    ll,
    rep(
      dinverse_burr(y[2], mu = mu, tau = tau, gamma = gamma,
                    log = TRUE),
      5
    ),
    eps = 1e-12
  )
  # the log_lik matches the injected Stan lpdf body
  scale <- mu * (2^(1 / tau) - 1)^(1 / gamma)
  w <- (y[2] / scale)^gamma
  stan_lpdf <- log(tau) + log(gamma) +
    tau * gamma * (log(y[2]) - log(scale)) -
    log(y[2]) -
    (tau + 1) * log1p(w)
  expect_eps(ll, rep(stan_lpdf, 5), eps = 1e-12)

  # posterior_predict returns ndraws draws on the positive scale
  pp <- posterior_predict_inverse_burr(2, prep)
  expect_length(pp, 5)
  expect_true(all(pp > 0))

  # posterior_epred returns the closed-form mean (exists for gamma > 1)
  pe <- posterior_epred_inverse_burr(prep)
  expect_eps(
    pe,
    matrix(inverse_burr_mean(mu, tau, gamma), nrow = 5, ncol = 3),
    eps = 1e-12,
    relative = TRUE
  )

  # posterior_epred errors for gamma <= 1 where the mean does not exist
  prep_no_mean <- mock_inverse_burr_prep(y, mu, tau, 0.5, ndraws = 5)
  expect_error(
    posterior_epred_inverse_burr(prep_no_mean),
    regexp = "gamma <= 1"
  )
})

test_that("inverse_burr is registered in the family registry", {
  expect_s3_class(brms_family_lookup("inverse_burr"), "brmsfamily")
  expect_equal(brms_family_lookup("inverse_burr"), inverse_burr())
  expect_identical(rng_lookup("inverse_burr"), rinverse_burr)
  expect_equal(
    aux_family_parameters_lookup("inverse_burr"),
    c("tau", "gamma")
  )
  expect_equal(
    aux_limits_lookup("inverse_burr"),
    list(lb = c(0, 0), ub = c(Inf, Inf))
  )
})

test_that("inverse_burr family construction", {
  fam <- inverse_burr()
  expect_s3_class(fam, "brmsfamily")
  expect_equal(fam$name, "inverse_burr")
  expect_equal(fam$dpars, c("mu", "tau", "gamma"))
  expect_equal(c(fam$link, fam$link_tau, fam$link_gamma), c("log", "log", "log"))
  expect_equal(unlist(fam$lb), c(mu = "0", tau = "0", gamma = "0"))
  scode <- fam$stanvars[[1]]$scode
  expect_true(grepl("inverse_burr_lpdf", scode))
  expect_true(grepl("inverse_burr_rng", scode))
})
