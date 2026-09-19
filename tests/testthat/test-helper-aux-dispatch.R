# Regression tests for issue #36:
# - test_rng() never defined expected_mus in the no-aux branch and pooled all
#   mus into a single rng call
# - construct_brms() called the rng a second time with aux_par = NA for
#   single-parameter distributions (two non-exclusive if branches)
#
# All fits are mocked, no MCMC is required.

mock_test_family <- function() {
  list(name = "mock_test_family", stanvars = NULL)
}

test_that("test_rng works without auxiliary parameters", {
  calls <- list()
  recording_rng <- function(n, mu) {
    calls[[length(calls) + 1]] <<- list(n = n, mu = mu)
    stats::rnorm(n, mean = mu, sd = 0.05)
  }

  expect_success(test_rng(
    rng_fun = recording_rng,
    metric_mu = mean,
    n = 100,
    mu_list = c(0, 1, 2),
    aux_list = NA,
    mu_eps = 0.1,
    p_acceptable_failures = 0
  ))

  # the rng has to be called exactly once per mu with only n and mu arguments
  expect_equal(length(calls), 3)
  expect_equal(vapply(calls, function(call) call$mu, numeric(1)), c(0, 1, 2))
  expect_true(all(vapply(calls, function(call) call$n, numeric(1)) == 100))

  # a biased rng has to be detected in the no-aux branch as well
  expect_failure(test_rng(
    rng_fun = function(n, mu) stats::rnorm(n, mean = mu + 5),
    metric_mu = mean,
    n = 100,
    mu_list = c(0, 1, 2),
    aux_list = NA,
    mu_eps = 0.1,
    p_acceptable_failures = 0,
    debug = FALSE
  ))
})

test_that("test_rng no-aux branch applies mu_link to the rng call only", {
  calls <- list()
  # shifts the linked mu back, so the comparison against the raw mus succeeds
  recording_rng <- function(n, mu) {
    calls[[length(calls) + 1]] <<- list(n = n, mu = mu)
    stats::rnorm(n, mean = mu - 100, sd = 0.05)
  }

  expect_success(test_rng(
    rng_fun = recording_rng,
    metric_mu = mean,
    n = 100,
    mu_list = c(0, 1, 2),
    aux_list = NA,
    mu_eps = 0.1,
    p_acceptable_failures = 0,
    mu_link = function(x) x + 100
  ))

  # the link was applied before calling the rng ...
  expect_equal(vapply(calls, function(call) call$mu, numeric(1)), c(100, 101, 102))
  # ... while expected_mus stayed the raw mu_list (verified by the success above)
})

test_that("test_rng dispatches by number of auxiliary parameters", {
  calls <- list()
  recording_rng <- function(n, mu, aux = NA_real_, aux2 = NA_real_) {
    calls[[length(calls) + 1]] <<- list(mu = mu, aux = aux, aux2 = aux2)
    stats::rnorm(n, mean = mu, sd = 0.05)
  }

  # one auxiliary parameter: aux outer loop, mu inner loop
  expect_success(test_rng(
    rng_fun = recording_rng,
    metric_mu = mean,
    n = 50,
    mu_list = c(1, 2),
    aux_list = c(10, 20),
    mu_eps = 0.1,
    p_acceptable_failures = 0
  ))
  expect_equal(length(calls), 4)
  expect_equal(
    vapply(calls, function(call) call$aux, numeric(1)),
    c(10, 10, 20, 20)
  )
  expect_equal(
    vapply(calls, function(call) call$mu, numeric(1)),
    c(1, 2, 1, 2)
  )

  # two auxiliary parameters: aux, aux2 outer loops, mu innermost
  calls <- list()
  expect_success(test_rng(
    rng_fun = recording_rng,
    metric_mu = mean,
    n = 50,
    mu_list = c(1, 2),
    aux_list = 3,
    aux2_list = c(7, 8),
    mu_eps = 0.1,
    p_acceptable_failures = 0
  ))
  expect_equal(length(calls), 4)
  combos <- vapply(
    calls,
    function(call) paste(call$aux, call$aux2, call$mu),
    character(1)
  )
  expect_equal(combos, c("3 7 1", "3 7 2", "3 8 1", "3 8 2"))

  # zero auxiliary parameters never pass an aux argument to the rng
  calls <- list()
  expect_success(test_rng(
    rng_fun = recording_rng,
    metric_mu = mean,
    n = 50,
    mu_list = c(1, 2),
    aux_list = NA,
    mu_eps = 0.1,
    p_acceptable_failures = 0
  ))
  expect_equal(length(calls), 2)
  expect_true(all(is.na(vapply(calls, function(call) call$aux, numeric(1)))))
  expect_true(all(is.na(vapply(calls, function(call) call$aux2, numeric(1)))))
})

test_that("construct_brms dispatches the rng by auxiliary parameter count", {
  brm_calls <- list()
  local_mocked_bindings(
    brm = function(formula, data, ...) {
      brm_calls[[length(brm_calls) + 1]] <<- list(formula = formula, data = data)
      structure(list(), class = "brmsfit")
    },
    .package = "brms"
  )

  rng_calls <- list()
  recording_rng <- function(...) {
    args <- list(...)
    rng_calls[[length(rng_calls) + 1]] <<- args
    seq_len(args[[1]])
  }

  # zero auxiliary parameters: exactly two arguments, no NA passed through
  construct_brms(
    n_data_sampels = 10,
    intercept = 2,
    aux_par = NA,
    aux2_par = NA,
    rng_link = identity,
    family = mock_test_family,
    rng = recording_rng
  )
  expect_equal(length(rng_calls), 1)
  expect_equal(length(rng_calls[[1]]), 2)
  expect_identical(rng_calls[[1]][[1]], 10)
  expect_identical(rng_calls[[1]][[2]], 2)

  # one auxiliary parameter: rng receives (n, linked intercept, aux_par)
  construct_brms(
    n_data_sampels = 10,
    intercept = 2,
    aux_par = 5,
    aux2_par = NA,
    rng_link = identity,
    family = mock_test_family,
    rng = recording_rng
  )
  expect_equal(length(rng_calls), 2)
  expect_equal(length(rng_calls[[2]]), 3)
  expect_identical(rng_calls[[2]][[3]], 5)

  # two auxiliary parameters: rng receives (n, linked intercept, aux, aux2)
  construct_brms(
    n_data_sampels = 10,
    intercept = 2,
    aux_par = 5,
    aux2_par = 7,
    rng_link = identity,
    family = mock_test_family,
    rng = recording_rng
  )
  expect_equal(length(rng_calls), 3)
  expect_equal(length(rng_calls[[3]]), 4)
  expect_identical(rng_calls[[3]][[3]], 5)
  expect_identical(rng_calls[[3]][[4]], 7)

  # rng_link is applied to the intercept before the rng call
  construct_brms(
    n_data_sampels = 10,
    intercept = 2,
    aux_par = NA,
    aux2_par = NA,
    rng_link = function(x) x + 100,
    family = mock_test_family,
    rng = recording_rng
  )
  expect_equal(length(rng_calls), 4)
  expect_identical(rng_calls[[4]][[2]], 102)

  # the mocked fit was called once per construction with the generated data
  expect_equal(length(brm_calls), 4)
  expect_equal(brm_calls[[1]]$data$y, 1:10)
  expect_s3_class(brm_calls[[1]]$formula, "formula")
})
