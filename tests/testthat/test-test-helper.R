test_that("test custom expect_eps function", {
  # wrong amount of arguments
  expect_error(expect_eps(1, 1))
  expect_error(expect_eps(1, 1, 1, 1, 1))
  # usage of non-numeric types, which is disallowed
  expect_error(expect_warning(expect_eps(c("a", 1), c("b", 1.1), c(2, 0.2))))
  # all scalars
  expect_success(expect_eps(1, 1.1, 0.2))
  expect_success(expect_eps(1, 3, 3))
  expect_success(expect_eps(-1, -1.1, 0.2))
  expect_failure(expect_eps(1, 2, 0.2))
  expect_failure(expect_eps(2, 1, 0.2))
  expect_error(expect_eps(0.1, 0.11, -0.1))
  # a vector, b scalar, eps scalar
  expect_success(expect_eps(c(1, 1), 1.1, 0.2))
  expect_success(expect_eps(c(2, 1, 1), 1.1, eps = 0.2, r = 0.4))
  expect_failure(expect_eps(c(1, 2), 1.1, 0.2))
  expect_failure(expect_eps(c(1, 1.1), 2, 0.2))
  expect_failure(expect_eps(c(2, 1, 1), 1.1, eps = 0.2, r = 0.2))
  # a vector, b vector, eps scalar
  expect_success(expect_eps(c(1, 1), c(1.1, 1.1), 0.2))
  expect_success(expect_eps(c(2, 1, 1), c(1.1, 1.1, 1.1), eps = 0.2, r = 0.4))
  expect_failure(expect_eps(c(1, 2), c(1.1, 1.1), 0.2))
  expect_failure(expect_eps(c(2, 1, 1), c(1.1, 1.1, 1.1), eps = 0.2, r = 0.2))
  expect_error(expect_eps(c(1, 1, 1), c(1.1, 1.1), 0.2))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1), 0.2))
  # a vector, b scalar, eps vector
  expect_success(expect_eps(c(1, 1), 1.1, c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 1, 2), 1.1, eps = c(0.2, 0.3, 0.2), r = 0.4))
  expect_success(expect_eps(c(1, 1), 1.1, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 2), 1.1, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 1), 2, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 1, 2), 1.1, eps = c(0.2, 0.3, 0.2), r = 0.2))
  expect_error(expect_eps(c(1, 1, 1), 1.1, c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 2, 3), 1.1, c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 1), 1.1, c(0.2, -0.3)))
  # all vectors
  expect_success(expect_eps(c(1, 1), c(1.1, 1.1), c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 3), c(1.1, 3.1), c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 5), c(1.1, 7), c(0.2, 3)))
  expect_success(expect_eps(
    c(1, 1, 2),
    c(1.1, 1.1, 1.1),
    eps = c(0.2, 0.3, 0.2),
    r = 0.4
  ))
  expect_failure(expect_eps(c(1, 2), c(1.1, 1.1), c(0.2, 0.3)))
  expect_failure(expect_eps(
    c(1, 1, 2),
    c(1.1, 1.1, 1.1),
    eps = c(0.2, 0.3, 0.2),
    r = 0.2
  ))
  expect_error(expect_eps(c(1, 1, 1), c(1.1, 1.1), c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1), c(0.2, 0.3)))
  # additional r-error-tests
  expect_error(expect_eps(
    c(1, 1, 2),
    c(1.1, 1.1, 1.1),
    eps = c(0.2, 0.3, 0.2),
    r = -0.1
  ))
  expect_error(expect_eps(
    c(1, 1, 2),
    c(1.1, 1.1, 1.1),
    eps = c(0.2, 0.3, 0.2),
    r = 1.1
  ))
  expect_error(expect_eps(
    c(1, 1, 2),
    c(1.1, 1.1, 1.1),
    eps = c(0.2, 0.3, 0.2),
    r = c(0.1, 0.2, 0.3)
  ))
  expect_error(expect_eps(1, 1.1, 0.2, "r"))

  # now all normal use cases, but as relative test
  expect_success(expect_eps(1, 1.1, 0.2, relative = TRUE))
  expect_success(expect_eps(1, 3, 0.8, relative = TRUE))
  expect_success(expect_eps(-1, -1.1, 0.2, relative = TRUE))
  expect_failure(expect_eps(1, 2, 0.2, relative = TRUE))
  expect_failure(expect_eps(2, 1, 0.2, relative = TRUE))
  expect_error(expect_eps(0.1, 0.11, -0.1, relative = TRUE))
  expect_error(expect_eps(0.1, 0.11, 1, relative = TRUE))
  # a vector, b scalar, eps scalar
  expect_success(expect_eps(c(1, 1), 1.1, 0.2, relative = TRUE))
  expect_success(expect_eps(
    c(2, 1, 1),
    1.1,
    eps = 0.2,
    r = 0.4,
    relative = TRUE
  ))
  expect_failure(expect_eps(c(1, 2), 1.1, 0.2, relative = TRUE))
  expect_failure(expect_eps(c(1, 1.1), 2, 0.2, relative = TRUE))
  expect_failure(expect_eps(
    c(2, 1, 1),
    1.1,
    eps = 0.2,
    r = 0.2,
    relative = TRUE
  ))
  # a vector, b vector, eps scalar
  expect_success(expect_eps(c(1, 1), c(1.1, 1.1), 0.2, relative = TRUE))
  expect_success(expect_eps(
    c(2, 1, 1),
    c(1.1, 1.1, 1.1),
    eps = 0.2,
    r = 0.4,
    relative = TRUE
  ))
  expect_failure(expect_eps(c(1, 2), c(1.1, 1.1), 0.2, relative = TRUE))
  expect_failure(expect_eps(
    c(2, 1, 1),
    c(1.1, 1.1, 1.1),
    eps = 0.2,
    r = 0.2,
    relative = TRUE
  ))
  expect_error(expect_eps(c(1, 1, 1), c(1.1, 1.1), 0.2, relative = TRUE))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1), 0.2, relative = TRUE))
  # a vector, b scalar, eps vector
  expect_success(expect_eps(c(1, 1), 1.1, c(0.2, 0.3), relative = TRUE))
  expect_success(expect_eps(
    c(1, 1, 2),
    1.1,
    eps = c(0.2, 0.3, 0.2),
    r = 0.4,
    relative = TRUE
  ))
  expect_success(expect_eps(c(1, 1), 1.1, c(0.2, 0.3), relative = TRUE))
  expect_failure(expect_eps(c(1, 2), 1.1, c(0.2, 0.3), relative = TRUE))
  expect_failure(expect_eps(c(1, 1), 2, c(0.2, 0.3), relative = TRUE))
  expect_failure(expect_eps(
    c(1, 1, 2),
    1.1,
    eps = c(0.2, 0.3, 0.2),
    r = 0.2,
    relative = TRUE
  ))
  expect_error(expect_eps(c(1, 1, 1), 1.1, c(0.2, 0.3), relative = TRUE))
  expect_error(expect_eps(c(1, 2, 3), 1.1, c(0.2, 0.3), relative = TRUE))
  expect_error(expect_eps(c(1, 1), 1.1, c(0.2, -0.3), relative = TRUE))
  # all vectors
  expect_success(expect_eps(c(1, 1), c(1.1, 1.1), c(0.2, 0.3), relative = TRUE))
  expect_success(expect_eps(c(1, 3), c(1.1, 3.1), c(0.2, 0.3), relative = TRUE))
  expect_success(expect_eps(c(1, 5), c(1.1, 7), c(0.1, 0.3), relative = TRUE))
  expect_success(expect_eps(
    c(1, 1, 2),
    c(1.1, 1.1, 1.1),
    eps = c(0.2, 0.3, 0.2),
    r = 0.4,
    relative = TRUE
  ))
  expect_failure(expect_eps(c(1, 2), c(1.1, 1.1), c(0.2, 0.3), relative = TRUE))
  expect_failure(expect_eps(
    c(1, 1, 2),
    c(1.1, 1.1, 1.1),
    eps = c(0.2, 0.3, 0.2),
    r = 0.2,
    relative = TRUE
  ))
  expect_error(expect_eps(
    c(1, 1, 1),
    c(1.1, 1.1),
    c(0.2, 0.3),
    relative = TRUE
  ))
  expect_error(expect_eps(
    c(1, 1, 2),
    c(1.1, 1.1),
    c(0.2, 0.3),
    relative = TRUE
  ))
})

test_that("test normal_difference metric-function", {
  # both scalars
  res1 <- normale_difference(1, 1.5)
  expect_equal(length(res1), 1)
  expect_eps(res1, 0.2773501, eps = 1e-6)
  # first vector, second scalar
  res2 <- normale_difference(c(1, 1.5, 2), 1.5)
  expect_equal(length(res2), 3)
  expect_eps(res2, c(0.2773501, 0, 0.2), eps = 1e-6)
  # now other way round
  res3 <- normale_difference(1.5, c(1, 1.5, 2))
  expect_equal(res2, res3)
  # just exchanging the arguments should not change anything in the results
  res4 <- normale_difference(c(0.5, 1, 1.5), c(1, 1.5, 2))
  expect_equal(length(res4), 3)
  expect_eps(res4, c(0.4472136, 0.2773501, 0.2), eps = 1e-6)

  # so much for the expected use-cases, now check edge cases and errors
  expect_equal(normale_difference(0, 0), 0)
  # 0,0 results in a small exception, where the result is not defined.
  # proof, that the exception is handled
  expect_eps(normale_difference(1e-15, 1e-15), 0, eps = 1e-15)
  # getting close to the same exception makes no trouble
  expect_error(normale_difference(c(2, 3), c(4, 5, 6)))
  # different lengths are dissallowed!
  expect_error(expect_warning(normale_difference(c(1, 2), c(1, "2"))))
  # all data has to be numeric, will produce warning
  expect_error(expect_warning(normale_difference(c(1, 2), c(1, NA))))
  # NAs are dissallowed as well, will produce warning
})

test_that("test the test_rng-wrapper", {
  # Test-RNG constants
  eps <- 1e-6
  n_small <- 10
  n <- 10000
  accepted_means_eps <- 0.35
  p_acceptable_failures <- 0.05
  mus <- seq(from = 1 + eps, to = 10, length.out = n_small)
  alphas_r <- seq(from = 2 + eps, to = 10, length.out = n_small)

  # the one test, if it works (not much point, in checking any other expected state)
  expect_success(test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = n,
    mu_list = mus,
    aux_list = alphas_r,
    mu_eps = accepted_means_eps,
    p_acceptable_failures = p_acceptable_failures
  ))
  # if the margins are too low, the function has a high likelihood to fail
  # for eps=0, will always fail
  expect_failure(test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = n,
    mu_list = mus,
    aux_list = alphas_r,
    mu_eps = 0,
    p_acceptable_failures = 0,
    debug = FALSE
  ))

  # else check all forbidden arguments
  # non-function type function arguments
  expect_error(test_rng(
    0,
    metric_mu = mean,
    n = n,
    mu_list = mus,
    aux_par = alphas_r,
    mu_eps = 0,
    p_acceptable_failures = 0
  ))
  expect_error(test_rng(
    rng_fun = rlomax,
    0,
    n = n,
    mu_list = mus,
    aux_par = alphas_r,
    mu_eps = 0,
    p_acceptable_failures = 0
  ))
  # switched function arguments
  expect_error(test_rng(
    rng_fun = mean,
    metric_mu = rlomax,
    n = n,
    mu_list = mus,
    aux_par = alphas_r,
    mu_eps = 0,
    p_acceptable_failures = 0
  ))
  # non numeric sample number argument
  expect_error(test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = "R",
    mu_list = mus,
    aux_par = alphas_r,
    mu_eps = accepted_means_eps,
    p_acceptable_failures = p_acceptable_failures
  ))
  # vector of numbers to sample
  expect_error(test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = c(42, 73),
    mu_list = mus,
    aux_par = alphas_r,
    mu_eps = accepted_means_eps,
    p_acceptable_failures = p_acceptable_failures
  ))
  # non numeric eps argument
  # expect_warning(expect_error(test_rng(
  #   rng_fun = rlomax, metric_mu = mean, n = n, mu_list = mus, aux_par = alphas_r,
  #   mu_eps = "R", p_acceptable_failures = p_acceptable_failures
  # )))
  expect_error(expect_warning(test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = n,
    mu_list = mus,
    aux_par = alphas_r,
    mu_eps = "R",
    p_acceptable_failures = p_acceptable_failures
  )))
  # negative eps argument
  expect_error(test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = n,
    mu_list = mus,
    aux_par = alphas_r,
    mu_eps = -1,
    p_acceptable_failures = p_acceptable_failures
  ))
  # vector eps argument
  expect_error(test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = n,
    mu_list = mus,
    aux_par = alphas_r,
    mu_eps = c(0.1, 0.2),
    p_acceptable_failures = p_acceptable_failures
  ))
  # non numeric p argument
  expect_error(test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = n,
    mu_list = mus,
    aux_par = alphas_r,
    mu_eps = accepted_means_eps,
    p_acceptable_failures = "R"
  ))
  # too small p argument (smaller 0)
  expect_error(test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = n,
    mu_list = mus,
    aux_par = alphas_r,
    mu_eps = accepted_means_eps,
    p_acceptable_failures = -1
  ))
  # too big p argument (bigger 1)
  expect_error(test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = n,
    mu_list = mus,
    aux_par = alphas_r,
    mu_eps = accepted_means_eps,
    p_acceptable_failures = 2
  ))
  # vector p argument
  expect_error(test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = n,
    mu_list = mus,
    aux_par = alphas_r,
    mu_eps = accepted_means_eps,
    p_acceptable_failures = c(0.1, 0.2)
  ))
})

test_that("test_brms_quantile", {
  # To test this, first construct a brms model with Cloglognormal
  n_brms <- 1000
  intercept <- 0.3
  sigma <- 0.6
  tresh <- 0.05

  # save old seed, to reset it later
  old_seed <- .Random.seed
  # Set predefined seed. Generating correct and "random" RNG data is not part of the brms recovery test.
  set.seed(9001)
  cloglog_data <- rcloglognormal(n_brms, intercept, sigma)
  set.seed(old_seed)
  # Now that the data was generated, reset the old seed (as if nothing ever happened)

  # limit the interval. Cloglognormal brms is very sensitive for data at the boundary.
  eps_brms <- 1e-12
  allowed_interval <- c(eps_brms, 1 - eps_brms)
  cloglog_data <- limit_data(cloglog_data, allowed_interval)

  fit <- brms::brm(
    y ~ 1,
    family = cloglognormal(),
    stanvars = cloglognormal()$stanvars,
    data = list(y = cloglog_data),
    backend = "rstan",
    cores = 2,
    chains = 2,
    silent = 2,
    refresh = 0,
    init = 0.1
  )

  # OK, after all that preamble, now it is getting interesting!
  expect_true(
    test_brms_quantile(fit, "b_Intercept", intercept, tresh) &&
      test_brms_quantile(fit, "sigma", sigma, tresh)
  )
  # This test should be correct (is the almost the same as in the Cloglog Testthat)
  expect_warning(expect_false(test_brms_quantile(fit, "alpha", sigma, tresh)))
  # No alpha in Cloglognormal, which return false and throws a warning
  sigma_data <- posterior::extract_variable_matrix(fit, variable = "sigma")
  median_sigma <- median(sigma_data)
  expect_false(test_brms_quantile(
    fit,
    "sigma",
    2 * tresh + 2 * median_sigma,
    tresh
  ))
  # definitively data not within quantiles
  expect_error(test_brms_quantile())
  # wrong amount of arguments
  expect_error(test_brms_quantile(c(1, 2, 3), "sigma", sigma, tresh))
  # vector is not type brms (which is a R list)
  expect_error(test_brms_quantile(fit, sigma, sigma, tresh))
  # sigma is not a string
  expect_error(test_brms_quantile(fit, c("sigma", "b_Intercept"), sigma, tresh))
  # only one value at a time
  expect_error(test_brms_quantile(fit, "sigma", "sigma", tresh))
  # reference value should be real scalar, not string
  expect_error(test_brms_quantile(fit, "sigma", NA, tresh))
  # reference value may not be NA
  expect_error(test_brms_quantile(fit, "sigma", c(sigma, b_Intercept), tresh))
  # reference value should not be a vector
  expect_error(test_brms_quantile(fit, "sigma", sigma, "tresh"))
  # threshold has to be of type real
  expect_error(test_brms_quantile(fit, "sigma", sigma, -1))
  # threshold has to be in the unit-interval
  expect_error(test_brms_quantile(fit, "sigma", sigma, c(0.1, 0.5, 0.9)))
  # threshold has to have 2 entries at most
  expect_error(test_brms_quantile(fit, "sigma", sigma, c()))
  # no entries for threshold is forbidden as well
  expect_error(test_brms_quantile(fit, "sigma", sigma, NA))
  # threshold cannot be NA
  expect_error(test_brms_quantile(fit, "sigma", sigma, 0.6))
  # for 0.6, the threshold would create bounds, with lowerbound > upperbound
  expect_error(test_brms_quantile(fit, "sigma", sigma, thresh, debug = "TRUE"))
  # debug has to be of type boolean
  expect_error(test_brms_quantile(
    fit,
    "sigma",
    sigma,
    thresh,
    debug = c(TRUE, FALSE)
  ))
  # debug has to be a single boolean
  expect_error(test_brms_quantile(fit, "sigma", sigma, thresh, debug = 0))
  # debug has to be type boolean
})


# Test a few fúrther miscellanios helping functions
test_that("Real of length n isNum_len", {
  # correct lengths and all numerics
  expect_true(isNum_len(0.2))
  expect_true(isNum_len(0.2, 1))
  expect_true(isNum_len(c(0.2, 0.3), 2))
  # incorrect lengths and all numerics
  expect_false(isNum_len(0.2, 2))
  expect_false(isNum_len(c(0.2, 0.3)))
  expect_false(isNum_len(c(0.2, 0.3), 1))
  # correct lengths, but not numerics
  expect_false(isNum_len("r"))
  expect_warning(expect_false(isNum_len(isNum_len)))
  expect_false(isNum_len(c(0.2, NA), 2)) # numeric vectors containing NAs are
  # usually also called numeric, but should not be set as such
  expect_false(isNum_len(c("r", 0.2), 2))
  expect_false(isNum_len(c("r", 0.2), 2))
  # incorrect length and non numerics
  expect_false(isNum_len(c("r", 0.2), 3))
})

test_that("Integer of length n isInt_len", {
  # correct lengths and all integers
  expect_true(isInt_len(1))
  expect_true(isInt_len(1, 1))
  expect_true(isInt_len(c(1, 2), 2))
  # incorrect lengths and all integers
  expect_false(isInt_len(1, 2))
  expect_false(isInt_len(c(1, 2)))
  expect_false(isInt_len(c(1, 2), 1))
  # correct lengths, but not integers
  expect_false(isInt_len("r"))
  expect_false(isInt_len(c("r", 1), 2))
  expect_false(isInt_len(c(1.1, 1), 2))
  expect_warning(expect_false(isInt_len(isInt_len)))
  expect_false(isInt_len(c(1, NA), 2)) # numeric vectors containing NAs are
  # usually also called numeric, but should not be set as such
  expect_false(isInt_len(c("r", 1), 2))
  # incorrect lengths and non integers
  expect_false(isInt_len(c("r", 2), 3))

  # because it is almost similar, just check the edge cases of isNat here
  expect_false(isNat_len(-1))
  expect_true(isInt_len(-1))
  expect_true(isNat_len(0))
  expect_false(isNat_len(c(-1, 0)))
})

test_that("Boolean of length n isLogic_len", {
  # correct lengths and all boolean
  expect_true(isLogic_len(TRUE))
  expect_true(isLogic_len(TRUE, 1))
  expect_true(isLogic_len(c(TRUE, FALSE), 2))
  # incorrect lengths and all boolean
  expect_false(isLogic_len(TRUE, 2))
  expect_false(isLogic_len(c(TRUE, FALSE)))
  expect_false(isLogic_len(c(TRUE, FALSE), 1))
  # correct lengths, but not boolean
  expect_false(isLogic_len("r"))
  expect_warning(expect_false(isLogic_len(isLogic_len)))
  expect_false(isLogic_len(c(TRUE, 0), 2))
  # other languages may recognise integers as boolean, which R does not do (I think)
  expect_false(isLogic_len(c("r", TRUE), 2))
  # incorrect length and non boolean
  expect_false(isLogic_len(c("r", TRUE), 3))
})

test_that("Single string isSingleString", {
  # correct strings
  expect_true(isSingleString("r"))
  expect_true(isSingleString("abc"))
  expect_true(isSingleString(c("abc")))
  # correct strings, but not single string
  expect_false(isSingleString(c("abc", "r")))
  # not a string
  expect_false(isSingleString(1))
  expect_warning(expect_false(isSingleString(isSingleString)))
  expect_false(isSingleString(c(NA, "abc")))
  # multiple inputs containing non strings
  expect_false(isSingleString(c("abc", 1, isSingleString)))
})

test_that("lenEqual length and type_check assertion", {
  # some test-vectors
  va <- c(1, 2, 3)
  vb <- c(NA, 5, 6)
  vc <- c(7, 8)
  vd <- c("r", 9, 10)
  ve <- c(12, 13, 14)
  # and test scalars
  sa <- 11
  sb <- "r"

  # now test it all
  # simple test, only checks lengths with all defaults
  expect_true(lenEqual(list(va, vd)))
  # with NAs included, fails with warning
  expect_warning(expect_false(lenEqual(list(va, vb, vd))))
  # now allow for NAs
  expect_true(lenEqual(list(va, vb, vd), na_allowed = TRUE))
  # now check numeric
  expect_true(lenEqual(list(va, ve), type_check = is.numeric))
  expect_warning(expect_false(lenEqual(
    list(va, vb, vd),
    type_check = is.numeric
  )))
  # now include scalars into the mix
  expect_true(lenEqual(list(va, ve, sa, sb), scalars_allowed = TRUE))
  expect_true(lenEqual(
    list(va, ve, sa),
    scalars_allowed = TRUE,
    type_check = is.numeric
  ))
  # non numerics
  expect_warning(expect_false(lenEqual(
    list(va, ve, vd, sa),
    scalars_allowed = TRUE,
    type_check = is.numeric
  )))
  expect_warning(expect_false(lenEqual(
    list(va, ve, sa, sb),
    scalars_allowed = TRUE,
    type_check = is.numeric
  )))
  expect_true(lenEqual(list(va, ve, sa, sb), scalars_allowed = TRUE))
  # NAs in vector still makes it a numeric
  expect_true(lenEqual(
    list(va, vb, sa),
    scalars_allowed = TRUE,
    type_check = is.numeric,
    na_allowed = TRUE
  ))
  # NA in scalar however will fail the is.numeric. Small peculiar quirk of R
  expect_warning(expect_false(lenEqual(
    list(va, NA, sa),
    scalars_allowed = TRUE,
    type_check = is.numeric,
    na_allowed = TRUE
  )))

  # so far for all the use cases, with correct lengths.
  expect_false(lenEqual(list(va, vc)))
  # usually scalars are dissalowed, hence they will have the wrong length
  expect_false(lenEqual(list(va, sa)))
  # even with scalars allowed, the actual non scalars have to be of same length
  expect_false(lenEqual(list(va, vc, sa), scalars_allowed = TRUE))
  # those two might look trivial, but depending on implementation, one might actually
  # break with incorrect length, before checking the wrong type (hence why no warning)
  expect_warning(expect_false(lenEqual(
    list(va, vc, vd),
    type_check = is.numeric
  )))
  expect_warning(expect_false(lenEqual(
    list(va, vd, vc),
    type_check = is.numeric
  )))
  # but this implementation does check in both cases.

  # I suppose, one might still find about 1000 different relevant permutations.
  # But I think, those should cover about 90% of all really important use-cases.
})

test_that("data limiting function limit_data", {
  input <- c(1, 2, 3, 4)
  result_both <- limit_data(input, c(2, 3))
  # length should not be changed
  expect_true(length(result_both) == 4)
  expect_true(all(2 <= result_both | result_both <= 3))
  # now only change the upper boundary
  result_upper <- limit_data(input, c(NA, 3))
  expect_true(all(1 <= result_upper | result_upper <= 3))
  # also the first input should not be changed now
  expect_true(input[1] == result_upper[1])
  # ... while the last one was set to 3
  expect_true(3 == result_upper[4])
  # should not work for limits beeing wrong way around
  expect_error(limit_data(input, c(3, 2)))
  # wrong amount of limit arguments
  expect_error(limit_data(input, c(1, 2, 3)))
  expect_error(limit_data(input, c(1)))
  # wrong input data
  expect_error(limit_data(c(1, "r"), c(1, 2)))
  expect_error(limit_data(c(1, NA), c(1, 2)))
})
