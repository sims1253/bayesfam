# Regression tests for issue #37:
# - construct_brms restored .Random.seed via set.seed(old_seed), which
#   initializes a new stream instead of restoring the caller's state
# - .Random.seed was accessed unconditionally and errored in fresh sessions
# - without on.exit, an erroring rng or fit leaked the helper's stream

mock_test_family <- function() {
  list(name = "mock_test_family", stanvars = NULL)
}

test_that("construct_brms restores the full caller RNG state after success", {
  local_mocked_bindings(
    brm = function(...) structure(list(), class = "brmsfit"),
    .package = "brms"
  )

  set.seed(42)
  state_before <- .Random.seed
  expected_next <- runif(3)
  # rewind to the saved state, then let construct_brms consume the stream
  assign(".Random.seed", state_before, envir = globalenv())

  construct_brms(
    n_data_sampels = 5,
    intercept = 1,
    aux_par = NA,
    aux2_par = NA,
    rng_link = identity,
    family = mock_test_family,
    rng = function(n, mu) stats::rnorm(n, mean = mu),
    seed = 99
  )

  # the full state vector is restored, so the caller's stream continues
  expect_identical(.Random.seed, state_before)
  expect_equal(runif(3), expected_next)
})

test_that("construct_brms restores the RNG state after an erroring rng", {
  set.seed(7)
  state_before <- .Random.seed

  expect_error(construct_brms(
    n_data_sampels = 5,
    intercept = 1,
    aux_par = NA,
    aux2_par = NA,
    rng_link = identity,
    family = mock_test_family,
    rng = function(n, mu) stop("rng boom"),
    seed = 99
  ))
  expect_identical(.Random.seed, state_before)
})

test_that("construct_brms restores the RNG state after an erroring fit", {
  local_mocked_bindings(
    brm = function(...) stop("brms boom"),
    .package = "brms"
  )

  set.seed(7)
  state_before <- .Random.seed

  expect_error(construct_brms(
    n_data_sampels = 5,
    intercept = 1,
    aux_par = NA,
    aux2_par = NA,
    rng_link = identity,
    family = mock_test_family,
    rng = function(n, mu) stats::rnorm(n, mean = mu),
    seed = 99
  ))
  expect_identical(.Random.seed, state_before)
})

test_that("construct_brms removes a helper-created seed in fresh sessions", {
  local_mocked_bindings(
    brm = function(...) structure(list(), class = "brmsfit"),
    .package = "brms"
  )

  # simulate a fresh session in which no RNG was used yet
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }

  construct_brms(
    n_data_sampels = 5,
    intercept = 1,
    aux_par = NA,
    aux2_par = NA,
    rng_link = identity,
    family = mock_test_family,
    rng = function(n, mu) stats::rnorm(n, mean = mu),
    seed = 99
  )

  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("construct_brms without a seed leaves the RNG stream untouched", {
  local_mocked_bindings(
    brm = function(...) structure(list(), class = "brmsfit"),
    .package = "brms"
  )

  set.seed(1)
  state_before <- .Random.seed

  construct_brms(
    n_data_sampels = 5,
    intercept = 1,
    aux_par = NA,
    aux2_par = NA,
    rng_link = identity,
    family = mock_test_family,
    rng = function(n, mu) stats::rnorm(n, mean = mu),
    seed = NULL
  )

  # no seeding also means no restoring, the rng draws advanced the stream
  expect_false(identical(.Random.seed, state_before))
})
