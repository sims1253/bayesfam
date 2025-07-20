test_that("unit-lindley", {
  n <- 10
  x <- seq(from = 0.1, to = 0.9, length.out = n)
  # Check lengths
  expect_equal(n, length(dunit_lindley(x, mu = 0.5)))
  expect_equal(n, length(qunit_lindley(x, mu = 0.5)))
  expect_equal(n, length(runit_lindley(n, mu = 0.5)))

  warning("No reference density to test against")

  # Check density function for errors
  expect_error(dunit_lindley(0.5)) # to few arguments
  expect_error(dunit_lindley(0.5, 0.5, 0.5, 0.5)) # to many arguments
  expect_error(dunit_lindley(-1, mu = 0.5)) # unit x
  expect_error(dunit_lindley(0.5, mu = -1)) # unit mu
  expect_error(dunit_lindley("r", mu = 0.5)) # non-numeric arguments are disallowed

  # Check quantiles for errors
  expect_error(qunit_lindley(0.5)) # to few arguments
  expect_error(qunit_lindley(0.5, 0.5, 0.5, 0.5)) # to many arguments
  expect_error(qunit_lindley(-1, mu = 0.5)) # unit x
  expect_error(qunit_lindley(0.5, mu = -1)) # unit mu
  expect_error(qunit_lindley("r", mu = 0.5)) # non-numeric arguments are disallowed

  # Check rng for errors
  expect_error(runit_lindley()) # to few arguments
  expect_error(runit_lindley(10, 0.8, 0.8, 0.8)) # to many arguments
  expect_error(runit_lindley(-1, mu = 0.8)) # number of drawn samples cannot be smaller 0
  # expect_error(runit_lindley("r", mu = 0.8)) # non-numeric arguments are disallowed
  expect_error(runit_lindley(0.5, mu = -1)) # unit mu

  for (i in seq(from = 0.1, to = 0.9, length.out = 9)) {
    s <- runit_lindley(1000000, i)
    expect_eps(i, mean(s), eps = 0.001)
  }

  data <- list(y = runit_lindley(n = 1000, mu = 0.2))
  fit <- brms::brm(
    formula = y ~ 1,
    data = data,
    family = unit_lindley(),
    stanvars = unit_lindley()$stanvars,
    refresh = 0,
    silent = 2
  )
  qs <- inv_logit(brms::fixef(fit)[c(3, 4)])
  expect_equal(0.2 > qs[1] && 0.2 < qs[2], TRUE)
})
