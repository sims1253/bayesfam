# Exploratory tests for issue #18: assess RNG means via bootstrap cumulative
# means. These only exercise the new helper on well-known RNGs, the regular
# family test files keep using the fixed-n approach for now.

test_that("bootstrap cumulative means cover a well-behaved rng (rnorm)", {
  result <- rng_bootstrap_cummean(
    rng_fun = function(n, mu) stats::rnorm(n, mean = mu),
    mu_list = c(-1, 0, 2.5),
    n = 2000,
    checkpoints = c(100, 2000),
    n_boot = 400,
    seed = 1234
  )

  # one row per combination and checkpoint
  expect_equal(nrow(result), 6)
  expect_setequal(result$n_checkpoint, c(100, 2000))
  expect_setequal(result$mu, c(-1, 0, 2.5))
  expect_true(all(is.na(result$aux)))

  # every combination has to be covered at the final checkpoint
  final <- result[result$n_checkpoint == 2000, ]
  expect_true(all(final$covered))
  # the bootstrap intervals contain their own estimates
  expect_true(all(final$lower <= final$estimate & final$estimate <= final$upper))
})

test_that("bootstrap cumulative means cover a skewed rng (rexp)", {
  # rexp has no location parameter, wrap it so the mean equals mu
  result <- rng_bootstrap_cummean(
    rng_fun = function(n, mu) stats::rexp(n, rate = 1 / mu),
    mu_list = c(1, 2),
    n = 4000,
    checkpoints = c(100, 4000),
    n_boot = 400,
    seed = 42
  )

  final <- result[result$n_checkpoint == 4000, ]
  expect_true(all(final$covered))
})

test_that("bootstrap cumulative means detect a biased rng", {
  # the drawn sample has mean mu + 3, which no interval should cover
  # at a meaningful sample size
  result <- rng_bootstrap_cummean(
    rng_fun = function(n, mu) stats::rexp(n, rate = 1 / (mu + 3)),
    mu_list = c(1, 2),
    n = 2000,
    checkpoints = c(100, 2000),
    n_boot = 400,
    seed = 42
  )

  final <- result[result$n_checkpoint == 2000, ]
  expect_false(any(final$covered))
})

test_that("rng_bootstrap_cummean dispatches auxiliary parameters", {
  result <- rng_bootstrap_cummean(
    rng_fun = function(n, mu, aux) stats::rnorm(n, mean = mu, sd = aux),
    mu_list = c(0, 1),
    aux_list = c(0.5, 5),
    n = 2000,
    checkpoints = 2000,
    n_boot = 200,
    seed = 7
  )

  expect_equal(nrow(result), 4)
  expect_setequal(result$aux, c(0.5, 5))
  expect_true(all(result$covered))
})

test_that("rng_bootstrap_cummean restores the RNG state when seeded", {
  set.seed(11)
  state_before <- .Random.seed

  result <- rng_bootstrap_cummean(
    rng_fun = function(n, mu) stats::rnorm(n, mean = mu),
    mu_list = 1,
    n = 100,
    checkpoints = 50,
    n_boot = 50,
    seed = 99
  )

  expect_equal(nrow(result), 1)
  expect_identical(.Random.seed, state_before)
})

test_that("rng_bootstrap_cummean validates its arguments", {
  expect_error(rng_bootstrap_cummean(rng_fun = 0, mu_list = 1))
  expect_error(rng_bootstrap_cummean(
    rng_fun = stats::rnorm,
    mu_list = 1,
    n = -1
  ))
  expect_error(rng_bootstrap_cummean(
    rng_fun = function(n, mu) stats::rnorm(n, mu),
    mu_list = 1,
    n = 10,
    checkpoints = 11
  ))
  expect_error(rng_bootstrap_cummean(
    rng_fun = function(n, mu) stats::rnorm(n, mu),
    mu_list = 1,
    n_boot = 0
  ))
  expect_error(rng_bootstrap_cummean(
    rng_fun = function(n, mu) stats::rnorm(n, mu),
    mu_list = 1,
    conf_level = 1
  ))
})
