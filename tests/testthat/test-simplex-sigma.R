test_that("simplex rejects nonpositive sigma with a clear error (#27)", {
  # density
  expect_error(
    dsimplex(0.5, mu = 0.5, sigma = -1),
    "sigma > 0"
  )
  expect_error(
    dsimplex(0.5, mu = 0.5, sigma = 0),
    "sigma > 0"
  )
  expect_error(
    dsimplex(c(0.3, 0.7), mu = 0.5, sigma = c(1, -2)),
    "sigma > 0"
  )
  # rng
  expect_error(
    rsimplex(5, mu = 0.5, sigma = -1),
    "sigma > 0"
  )
  expect_error(
    rsimplex(5, mu = 0.5, sigma = 0),
    "sigma > 0"
  )
  # positive sigma keeps working
  expect_true(all(is.finite(dsimplex(c(0.2, 0.5), mu = 0.5, sigma = 2, log = TRUE))))
  expect_true(all(rsimplex(5, mu = 0.5, sigma = 2) > 0 &
    rsimplex(5, mu = 0.5, sigma = 2) < 1))
})

test_that("simplex family constrains sigma to be positive (#27)", {
  fam <- simplex()
  expect_equal(fam$dpars, c("mu", "sigma"))
  # API change: the sigma auxiliary parameter now defaults to a log link
  expect_equal(fam$link_sigma, "log")
  # sigma gets a lower bound of 0
  expect_equal(as.numeric(fam$lb$sigma), 0)
  expect_true(is.na(as.numeric(fam$ub$sigma)))
})

test_that("simplex model with a sigma predictor samples on the positive scale (#27)", {
  skip_on_cran()
  set.seed(271124)
  n <- 400
  x <- rnorm(n)
  dat <- data.frame(
    y = rsimplex(n, mu = plogis(0.3 + 0.5 * x), sigma = 1.5),
    x = x
  )
  fit <- brms::brm(
    formula = brms::bf(y ~ x, sigma ~ x),
    data = dat,
    family = simplex(),
    stanvars = simplex()$stanvars,
    chains = 2,
    cores = 2,
    iter = 500,
    warmup = 250,
    seed = 271124,
    refresh = 0,
    silent = 2
  )

  # with the log link, every posterior draw of sigma is strictly positive,
  # so the likelihood never evaluates log(negative)
  sigma_linpred <- brms::posterior_linpred(
    fit,
    dpar = "sigma",
    newdata = data.frame(x = c(0, 1))
  )
  expect_true(all(sigma_linpred > 0))
  # posterior_predict runs through the (now validating) RNG
  pred <- brms::posterior_predict(fit)
  expect_equal(dim(pred), c(posterior::ndraws(fit), n))
  expect_true(all(pred > 0 & pred < 1))
})
