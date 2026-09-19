test_that("posterior_predict_lomax samples from the Lomax distribution (#22)", {
  # Deterministic callback test: the callback must delegate to rlomax(),
  # not rgompertz().
  prep <- structure(list(
    ndraws = 1000L,
    dpars = list(mu = rep(1, 1000), alpha = rep(2, 1000))
  ), class = "brmsprep")
  set.seed(101)
  actual <- posterior_predict_lomax(1, prep)
  set.seed(101)
  expected <- rlomax(1000, mu = 1, alpha = 2)
  expect_equal(actual, expected)

  # The two candidate distributions disagree, so the seed test is meaningful.
  set.seed(101)
  expect_false(isTRUE(all.equal(
    posterior_predict_lomax(1, prep),
    rgompertz(1000, mu = 1, beta = 2)
  )))

  # Vectorized dpars are passed through elementwise.
  prep_vec <- structure(list(
    ndraws = 4L,
    dpars = list(mu = c(1, 2, 3, 4), alpha = c(2, 3, 4, 5))
  ), class = "brmsprep")
  set.seed(42)
  actual_vec <- posterior_predict_lomax(1, prep_vec)
  set.seed(42)
  expected_vec <- rlomax(4, mu = c(1, 2, 3, 4), alpha = c(2, 3, 4, 5))
  expect_equal(actual_vec, expected_vec)
})

test_that("posterior_predict on a fitted lomax model draws Lomax values (#22)", {
  skip_on_cran()
  set.seed(213471)
  n <- 500
  mu_true <- 3
  alpha_true <- 3
  dat <- data.frame(y = rlomax(n, mu = mu_true, alpha = alpha_true))
  fit <- brms::brm(
    formula = y ~ 1,
    data = dat,
    family = lomax(),
    stanvars = lomax()$stanvars,
    chains = 2,
    cores = 2,
    iter = 600,
    warmup = 300,
    seed = 213471,
    refresh = 0,
    silent = 2
  )

  pred <- brms::posterior_predict(fit)
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(posterior::ndraws(fit), n))
  expect_true(all(pred >= 0))

  # Predictions must follow the posterior predictive CDF obtained by mixing
  # the Lomax CDF over the posterior draws. A Gompertz RNG (the old bug)
  # fails this check by construction.
  draws <- brms::as_draws_df(fit)
  mu_draws <- exp(draws$Intercept) # log link for mu
  alpha_draws <- draws$alpha # already on the response scale
  lomax_cdf <- function(q) {
    1 - (1 + q / (mu_draws * (alpha_draws - 1)))^(-alpha_draws)
  }
  pp_cdf <- function(q) {
    vapply(
      q,
      function(qi) mean(lomax_cdf(qi)),
      numeric(1)
    )
  }
  ks <- stats::ks.test(as.vector(pred), pp_cdf)
  expect_gt(ks$p.value, 0.01)
})
