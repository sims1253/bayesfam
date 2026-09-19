test_that("dgeneralized_gamma is vectorized in Q (#23)", {
  x <- c(0.5, 1.2, 3.7, 10)
  mu <- c(0, 0.2, -0.5, 1.1)
  sigma <- c(1, 0.7, 2.3, 0.5)
  scalar_call <- function(xi, mi, si, Qi) {
    dgeneralized_gamma(xi, mu = mi, sigma = si, Q = Qi, log = TRUE)
  }

  # all-zero Q
  Q <- c(0, 0, 0, 0)
  expect_equal(
    dgeneralized_gamma(x, mu, sigma, Q = Q, log = TRUE),
    vapply(seq_along(x), function(i) scalar_call(x[i], mu[i], sigma[i], Q[i]), numeric(1))
  )

  # all-nonzero Q
  Q <- c(-1.5, 0.3, 2, 7.7)
  expect_equal(
    dgeneralized_gamma(x, mu, sigma, Q = Q, log = TRUE),
    vapply(seq_along(x), function(i) scalar_call(x[i], mu[i], sigma[i], Q[i]), numeric(1))
  )

  # mixed zero / nonzero Q
  Q <- c(0, 0.9, 0, -2.2)
  expect_equal(
    dgeneralized_gamma(x, mu, sigma, Q = Q, log = TRUE),
    vapply(seq_along(x), function(i) scalar_call(x[i], mu[i], sigma[i], Q[i]), numeric(1))
  )

  # same for the density itself
  Q <- c(0, 0.9, 0, -2.2)
  expect_equal(
    dgeneralized_gamma(x, mu, sigma, Q = Q, log = FALSE),
    exp(vapply(seq_along(x), function(i) scalar_call(x[i], mu[i], sigma[i], Q[i]), numeric(1)))
  )

  # scalars are recycled against longer vectors
  expect_equal(
    dgeneralized_gamma(x, mu = 0.2, sigma = 1.5, Q = 0.4, log = TRUE),
    dgeneralized_gamma(x, mu = rep(0.2, 4), sigma = rep(1.5, 4), Q = rep(0.4, 4), log = TRUE)
  )

  # the log_lik callback pattern: scalar y with length-S posterior draws
  S <- 7
  mu_s <- rnorm(S, 0.3, 0.1)
  sigma_s <- rep(1.2, S)
  Q_s <- c(0, 1.1, 0, -0.7, 2.3, 0, 0.5)
  expect_equal(
    dgeneralized_gamma(1.7, mu = mu_s, sigma = sigma_s, Q = Q_s, log = TRUE),
    vapply(seq_len(S), function(s) scalar_call(1.7, mu_s[s], sigma_s[s], Q_s[s]), numeric(1))
  )
})

test_that("brms log_lik and loo work for a fitted generalized gamma model (#23)", {
  skip_on_cran()
  set.seed(223341)
  n <- 500
  dat <- data.frame(y = rgeneralized_gamma(n, mu = 5, sigma = 2, Q = 3))
  fit <- brms::brm(
    formula = y ~ 1,
    data = dat,
    family = generalized_gamma(),
    stanvars = generalized_gamma()$stanvars,
    chains = 2,
    cores = 2,
    iter = 500,
    warmup = 250,
    seed = 223341,
    refresh = 0,
    silent = 2
  )

  ll <- brms::log_lik(fit)
  expect_true(is.matrix(ll))
  expect_equal(dim(ll), c(posterior::ndraws(fit), n))
  expect_true(all(is.finite(ll)))

  loo_res <- loo::loo(ll)
  expect_s3_class(loo_res, "loo")
  expect_true(all(is.finite(loo_res$estimates)))
})
