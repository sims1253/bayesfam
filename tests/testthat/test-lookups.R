custom_families <- c(
  "betaprime",
  "cauchitnormal",
  "cloglognormal",
  "generalized_gamma",
  "generalized_normal",
  "gompertz",
  "gumbel_mean",
  "inverse_burr",
  "kumaraswamy",
  "logistic",
  "logitnormal",
  "lognormal_natural",
  "lomax",
  "shifted_inv_gaussian",
  "shifted_lognormal_uniform",
  "simplex",
  "softplusnormal",
  "symlognormal",
  "unit_lindley"
)
builtin_families <- c(
  "beta",
  "frechet",
  "gamma",
  "gaussian",
  "inverse.gaussian",
  "lognormal",
  "weibull"
)

test_that("every family resolves in every lookup", {
  for (f in c(custom_families, builtin_families)) {
    family <- brms_family_lookup(f)
    expect_s3_class(family, "brmsfamily")
    expect_true(is.function(rng_lookup(f)))
    aux <- aux_family_parameters_lookup(f)
    expect_type(aux, "character")
    limits <- aux_limits_lookup(f)
    expect_length(limits$lb, length(aux))
    expect_length(limits$ub, length(aux))
    expect_true(is.function(link_lookup("identity", f)))
  }
})

test_that("lookups match direct construction", {
  for (f in custom_families) {
    expect_equal(brms_family_lookup(f), do.call(f, list()))
  }
  expect_equal(brms_family_lookup("gaussian"), brms::brmsfamily("gaussian"))
  expect_equal(
    brms_family_lookup("weibull", "softplus"),
    brms::brmsfamily("weibull", "softplus")
  )
  expect_equal(
    brms_family_lookup("kumaraswamy", "probit"),
    kumaraswamy(link = "probit")
  )
})

test_that("custom family bounds match the constructor", {
  to_bound <- function(x, default) {
    x <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)[-1]))
    replace(x, is.na(x), default)
  }
  for (f in custom_families) {
    family <- do.call(f, list())
    expect_equal(
      aux_limits_lookup(f),
      list(lb = to_bound(family$lb, -Inf), ub = to_bound(family$ub, Inf))
    )
  }
  expect_equal(aux_limits_lookup("lomax"), list(lb = 1, ub = Inf))
  expect_equal(
    aux_limits_lookup("shifted_lognormal_uniform"),
    list(lb = c(0, 0, 0), ub = c(Inf, 1, 1))
  )
  expect_equal(aux_limits_lookup("frechet"), list(lb = 1, ub = Inf))
})

test_that("aux parameters exclude mu and are empty without aux", {
  expect_identical(aux_family_parameters_lookup("gaussian"), "sigma")
  expect_identical(aux_family_parameters_lookup("unit_lindley"), character(0))
  expect_equal(
    aux_limits_lookup("unit_lindley"),
    list(lb = numeric(0), ub = numeric(0))
  )
})

test_that("unknown identifiers error clearly", {
  expect_error(brms_family_lookup("weibul"), "Unknown family 'weibul'")
  expect_error(rng_lookup("gamm"), "Unknown family 'gamm'")
  expect_error(aux_family_parameters_lookup("gaussion"), "Unknown family")
  expect_error(aux_limits_lookup(1), "Unknown family")
  expect_error(link_lookup("identiy"), "Unknown link 'identiy'")
})

test_that("link_lookup uses the built-in link of transformed families", {
  expect_identical(link_lookup("identity", "logitnormal"), logit)
  expect_identical(link_lookup("identity", "symlognormal", TRUE), inv_symlog)
  expect_identical(link_lookup("identity", "lognormal", TRUE), exp)
  # other or unregistered families leave the link unchanged
  expect_identical(link_lookup("log", "gaussian"), log)
  expect_identical(link_lookup("log", "student", TRUE), exp)
  expect_identical(link_lookup("symlog"), symlog)
})
