# Regression tests for issue #35:
# - normale_difference used to overflow to 0 or NaN for large finite inputs
# - expect_eps dropped NaN comparisons via na.rm, so grossly different
#   inputs could false-pass

test_that("normale_difference handles large finite values without overflow", {
  # squaring the raw values overflows to Inf, which used to yield a difference of 0
  expect_equal(normale_difference(1e200, 2e200), 1 / sqrt(5))
  # the metric has to stay symmetric
  expect_equal(normale_difference(2e200, 1e200), 1 / sqrt(5))
  # opposite signs near the double limit used to produce NaN
  expect_equal(normale_difference(1e308, -1e308), sqrt(2))
  expect_equal(normale_difference(-1e308, 1e308), sqrt(2))
  # one operand completely dwarfs the other
  expect_equal(normale_difference(1e308, 1e-308), 1)
  # vectors work element-wise, including mixed magnitudes
  expect_equal(
    normale_difference(c(1e200, 1e308, 1e150), c(2e200, -1e308, 3e150)),
    c(1 / sqrt(5), sqrt(2), 2 / sqrt(10))
  )
})

test_that("normale_difference handles infinite operands", {
  # equal infinities have zero distance
  expect_equal(normale_difference(Inf, Inf), 0)
  expect_equal(normale_difference(-Inf, -Inf), 0)
  # opposite-sign infinities are maximally apart in the euler metric
  expect_equal(normale_difference(Inf, -Inf), sqrt(2))
  expect_equal(normale_difference(-Inf, Inf), sqrt(2))
  # an infinite value against a finite one has normalized distance 1
  expect_equal(normale_difference(Inf, 5), 1)
  expect_equal(normale_difference(-5, Inf), 1)
  expect_equal(normale_difference(1e-300, -Inf), 1)
  # mixed vectors
  expect_equal(
    normale_difference(c(Inf, -Inf, 1, Inf), c(Inf, -Inf, 1, 2)),
    c(0, 0, 0, 1)
  )
})

test_that("normale_difference near-zero values", {
  expect_equal(normale_difference(0, 0), 0)
  expect_equal(normale_difference(-0.0, 0), 0)
  # denormal-range values normalize safely upwards
  expect_equal(normale_difference(1e-300, 2e-300), 1 / sqrt(5))
  expect_equal(normale_difference(1e-320, 1e-320), 0)
  expect_equal(normale_difference(c(0, 1e-200), c(0, 2e-200)), c(0, 1 / sqrt(5)))
})

test_that("expect_eps counts non-finite comparison results as failures", {
  # Inf - Inf is NaN in absolute mode and used to be dropped via na.rm,
  # silently passing the comparison
  expect_failure(expect_eps(Inf, Inf, 0.1))
  expect_failure(expect_eps(c(1, Inf), c(1, Inf), 0.1))
  expect_failure(expect_eps(c(Inf, 2), c(-Inf, 2), 10))
  # relative mode: opposite-sign infinities have distance sqrt(2) > eps
  expect_failure(expect_eps(Inf, -Inf, 0.99, relative = TRUE))
  # large finite mismatches used to false-pass through overflow
  expect_failure(expect_eps(1e200, 2e200, 0.1, relative = TRUE))
  expect_failure(expect_eps(1e308, -1e308, 0.1, relative = TRUE))
  expect_failure(expect_eps(c(1e200, 1e-200), c(2e200, 1e-200), 0.1, relative = TRUE))
  # a single non-finite entry inside a longer comparison is detected
  expect_failure(expect_eps(c(1e200, 1e200), c(1e200, 2e200), 0.1, relative = TRUE))
  # ... and can be tolerated explicitly via the r argument
  expect_success(
    expect_eps(c(1e200, 1e200), c(1e200, 2e200), 0.1, r = 0.5, relative = TRUE)
  )
  expect_failure(
    expect_eps(c(1e200, 1e200), c(1e200, 2e200), 0.1, r = 0.4, relative = TRUE)
  )
})

test_that("expect_eps still passes genuinely equal large, tiny and infinite inputs", {
  expect_success(expect_eps(1e300, 1e300, 1e-8, relative = TRUE))
  expect_success(expect_eps(1e-300, 1e-300, 1e-8, relative = TRUE))
  expect_success(
    expect_eps(c(1e200, 1e-200), c(1e200, 1e-200), 1e-8, relative = TRUE)
  )
  # ordinary finite values are unaffected by the overflow fix
  expect_success(expect_eps(1, 1.1, 0.2, relative = TRUE))
  expect_failure(expect_eps(1, 2, 0.2, relative = TRUE))
})
