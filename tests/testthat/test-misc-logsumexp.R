# Issue #33: logsumexp subtracted max(x) without special-casing infinities,
# so all-(-Inf) input produced NaN and +Inf input produced NaN. Impossible
# mixture observations in dshifted_lognormal_uniform inherited the NaN.

test_that("logsumexp agrees with the naive definition for finite inputs (#33)", {
  set.seed(20241002)
  for (i in seq_len(200)) {
    x <- rnorm(sample(seq_len(25), 1), mean = -5, sd = 15)
    # all representable on the direct scale here with high probability
    if (max(x) < 500 && length(x) > 1L) {
      expect_equal(logsumexp(x), log(sum(exp(x))), tolerance = 1e-12)
    }
  }
  expect_equal(logsumexp(c(0, 0)), log(2))
  expect_equal(logsumexp(5), 5)
  expect_equal(logsumexp(c(1, 2, 3)), log(exp(1) + exp(2) + exp(3)))
})

test_that("logsumexp handles infinities explicitly (#33)", {
  # all-impossible terms
  expect_equal(logsumexp(c(-Inf, -Inf)), -Inf)
  expect_equal(logsumexp(rep(-Inf, 10)), -Inf)
  # a single possible term dominates correctly
  expect_equal(logsumexp(c(-Inf, 0)), 0)
  expect_equal(logsumexp(c(-Inf, 100)), 100)
  expect_equal(logsumexp(c(-Inf, -Inf, 3)), 3)
  # any +Inf term gives +Inf (not NaN)
  expect_equal(logsumexp(c(Inf, 0)), Inf)
  expect_equal(logsumexp(c(-Inf, Inf)), Inf)
  # mixed finite and -Inf
  expect_equal(logsumexp(c(-Inf, 1, -Inf)), 1)
})

test_that("logsumexp defines NA/NaN and empty behavior (#33)", {
  expect_equal(logsumexp(c(1, NA)), NA_real_)
  expect_equal(logsumexp(c(1, NaN)), NA_real_)
  expect_equal(logsumexp(NA_real_), NA_real_)
  # empty sum is empty: log(0)
  expect_equal(logsumexp(numeric(0)), -Inf)
})

test_that("impossible mixture observations yield -Inf, not NaN (#33)", {
  # observation below the shift with zero uniform weight: impossible under
  # both mixture components
  expect_true(
    is.finite(dshifted_lognormal_uniform(0.5, mix = 0.1, shift = 1))
  )
  expect_equal(
    dshifted_lognormal_uniform(0.5, mix = 0, shift = 1),
    -Inf
  )
  # observation above max_uniform with zero lognormal weight: impossible
  expect_equal(
    dshifted_lognormal_uniform(200, mix = 1, max_uniform = 100),
    -Inf
  )
  # degenerate weights keep the possible component
  possible <- dlnorm(200, meanlog = 0, sdlog = 1, log = TRUE) -
    plnorm(100, meanlog = 0, sdlog = 1, log.p = TRUE)
  expect_equal(
    dshifted_lognormal_uniform(200, mix = 0, max_uniform = 100),
    possible
  )
})

test_that("mixture density stays column-wise vectorized (#33)", {
  y <- c(0.5, 200, 50)
  out <- dshifted_lognormal_uniform(
    y,
    meanlog = c(0, 0, 0),
    sdlog = c(1, 1, 2),
    mix = c(0, 1, 0.1),
    shift = c(1, 0, 0),
    max_uniform = c(100, 100, 100)
  )
  expect_equal(length(out), 3L)
  # column 1: below shift with mix = 0 -> impossible under both components
  expect_equal(out[[1]], -Inf)
  # column 2: above max_uniform with mix = 1 -> impossible under both components
  expect_equal(out[[2]], -Inf)
  # column 3: ordinary mixture observation -> finite density
  expect_true(is.finite(out[[3]]))
  # each column equals the elementwise computation
  for (j in seq_along(y)) {
    expect_equal(
      out[[j]],
      dshifted_lognormal_uniform(
        y[j],
        meanlog = c(0, 0, 0)[j],
        sdlog = c(1, 1, 2)[j],
        mix = c(0, 1, 0.1)[j],
        shift = c(1, 0, 0)[j],
        max_uniform = 100
      )
    )
  }
})
