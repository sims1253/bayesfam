test_that("unit-lindley", {

  n <- 10
  x <- seq(from=0.1, to=0.9, length.out=n)
  # Check lengths
  expect_equal(n, length(dunit_lindley(x, mu = 0.5)))
  expect_equal(n, length(qunit_lindley(x, mu = 0.5)))
  expect_equal(n, length(runit_lindley(n, mu = 0.5)))

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
  expect_error(runit_lindley("r", mu = 0.8)) # non-numeric arguments are disallowed
  expect_error(runit_lindley(0.5, mu = -1)) # unit mu

  warning("RNG test missing")

  data <- list(y = runit_lindley(n = 1000, mu = 0.2))
  fit <- brms::brm(formula = y ~ 1, data = data,
                   family = unit_lindley(), stanvars = unit_lindley()$stanvars,
                   refresh = 0)
})
