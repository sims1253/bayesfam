registry_families <- names(bayesfam:::family_registry)

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

test_that("registry covers all exported bayesfam families and common brms families", {
  for (f in custom_families) {
    expect_contains(registry_families, f)
  }
  for (f in builtin_families) {
    expect_contains(registry_families, f)
  }
})

test_that("every registry entry provides family, rng, aux names and bounds", {
  for (f in registry_families) {
    expect_error(brms_family_lookup(f), NA)
    expect_s3_class(brms_family_lookup(f), "brmsfamily")
    expect_true(is.function(rng_lookup(f)))
    aux <- aux_family_parameters_lookup(f)
    expect_type(aux, "character")
    limits <- aux_limits_lookup(f)
    expect_type(limits$lb, "double")
    expect_type(limits$ub, "double")
    expect_equal(length(limits$lb), length(aux))
    expect_equal(length(limits$ub), length(aux))
    expect_true(all(limits$lb <= limits$ub))
  }
})

test_that("lookup construction matches direct construction", {
  for (f in custom_families) {
    expect_equal(brms_family_lookup(f), do.call(f, list()))
  }
  expect_equal(
    brms_family_lookup("beta"),
    brms::brmsfamily("beta")
  )
  expect_equal(
    brms_family_lookup("weibull", "softplus"),
    brms::brmsfamily("weibull", "softplus")
  )
  expect_equal(
    brms_family_lookup("gaussian", "log"),
    brms::brmsfamily("gaussian", "log")
  )
  # aux links keep constructor defaults when only the mu link is given
  expect_equal(
    brms_family_lookup("kumaraswamy", "probit"),
    kumaraswamy(link = "probit")
  )
})

test_that("constructor defaults are preserved when link is NULL", {
  expect_equal(brms_family_lookup("gaussian")$link, "identity")
  expect_equal(brms_family_lookup("beta")$link, "logit")
  expect_equal(brms_family_lookup("kumaraswamy")$links, kumaraswamy()$links)
  expect_equal(brms_family_lookup("lomax")$links, lomax()$links)
})

test_that("newly exported families are covered", {
  expect_s3_class(brms_family_lookup("generalized_gamma"), "brmsfamily")
  expect_s3_class(brms_family_lookup("symlognormal"), "brmsfamily")
  expect_identical(rng_lookup("unit_lindley"), runit_lindley)
  expect_identical(rng_lookup("symlognormal"), rsymlognormal)
  expect_identical(rng_lookup("generalized_gamma"), rgeneralized_gamma)
})

test_that("zero-auxiliary families return empty aux information", {
  expect_equal(aux_family_parameters_lookup("unit_lindley"), character(0))
  expect_equal(
    aux_limits_lookup("unit_lindley"),
    list(lb = numeric(0), ub = numeric(0))
  )
})

test_that("single-dpar families return their aux, not c(NA, 'mu')", {
  expect_equal(aux_family_parameters_lookup("gaussian"), "sigma")
  expect_equal(aux_family_parameters_lookup("beta"), "phi")
  expect_equal(aux_family_parameters_lookup("weibull"), "shape")
  expect_false(any(is.na(aux_family_parameters_lookup("gaussian"))))
})

test_that("registry metadata for selected families", {
  expect_identical(rng_lookup("weibull"), rweibull_median)
  expect_identical(rng_lookup("frechet"), rfrechet_median)
  expect_identical(rng_lookup("beta"), rbeta_mean)
  expect_identical(rng_lookup("gamma"), rgamma_mean)
  expect_identical(rng_lookup("lognormal"), rlognormal)
  expect_identical(rng_lookup("inverse.gaussian"), brms::rinv_gaussian)
  expect_identical(rng_lookup("gaussian"), rnorm)
  expect_equal(aux_family_parameters_lookup("generalized_gamma"), c("sigma", "Q"))
  expect_equal(
    aux_family_parameters_lookup("shifted_lognormal_uniform"),
    c("sigma", "mix", "shiftprop")
  )
  expect_equal(aux_limits_lookup("lomax"), list(lb = 1, ub = Inf))
  expect_equal(
    aux_limits_lookup("shifted_lognormal_uniform"),
    list(lb = c(0, 0, 0), ub = c(Inf, 1, 1))
  )
})

test_that("misspelled identifiers error clearly", {
  expect_error(brms_family_lookup("weibul"), "weibul")
  expect_error(brms_family_lookup("generalized gamma"), "generalized gamma")
  expect_error(rng_lookup("gamm"), "gamm")
  expect_error(rng_lookup("unit_lindely"), "unit_lindely")
  expect_error(aux_family_parameters_lookup("gaussion"), "gaussion")
  expect_error(aux_limits_lookup("betaa"), "betaa")
  expect_error(brms_family_lookup(1), "single string")
})

test_that("link_lookup knows all links incl. symlog", {
  expect_equal(link_lookup("logit")(0.5), logit(0.5))
  expect_equal(link_lookup("cauchit")(0.5), cauchit(0.5))
  expect_equal(link_lookup("cloglog")(0.5), cloglog(0.5))
  expect_equal(link_lookup("identity")(0.5), 0.5)
  expect_equal(link_lookup("log")(1), 0)
  expect_equal(link_lookup("softplus")(2), softplus(2))
  expect_equal(link_lookup("symlog")(1.5), symlog(1.5))
})

test_that("link_lookup response functions are the link inverses", {
  test_points <- list(
    logit = c(0.1, 0.3, 0.5, 0.7, 0.9),
    cauchit = c(0.1, 0.3, 0.5, 0.7, 0.9),
    cloglog = c(0.1, 0.3, 0.5, 0.7, 0.9),
    identity = c(-3, -1, 1, 3),
    log = c(0.5, 1, 2, 5),
    softplus = c(0.5, 1, 2, 5),
    symlog = c(-3, -1, 1, 3)
  )
  for (name in names(bayesfam:::family_registry_links)) {
    x <- test_points[[name]]
    expect_eps(
      link_lookup(name, inv = TRUE)(link_lookup(name)(x)),
      x,
      eps = 1e-10
    )
  }
})

test_that("link_lookup family overrides", {
  expect_equal(link_lookup("identity", "logitnormal", FALSE)(0.5), logit(0.5))
  expect_equal(link_lookup("identity", "cauchitnormal", TRUE)(1), inv_cauchit(1))
  expect_equal(link_lookup("identity", "cloglognormal", FALSE)(0.5), cloglog(0.5))
  expect_equal(link_lookup("identity", "softplusnormal", TRUE)(1), inv_softplus(1))
  expect_equal(link_lookup("identity", "symlognormal", FALSE)(1), symlog(1))
  expect_equal(link_lookup("identity", "symlognormal", TRUE)(1), inv_symlog(1))
  expect_equal(link_lookup("identity", "lognormal", FALSE)(2), log(2))
  # families without an override use the given link
  expect_equal(link_lookup("log", "gaussian", FALSE)(2), log(2))
  expect_equal(link_lookup("identity", "gaussian", FALSE)(2), 2)
  expect_equal(link_lookup("log", "gaussian", TRUE)(2), exp(2))
})

test_that("link_lookup errors on unknown link or family identifiers", {
  expect_error(link_lookup("logit2"), "logit2")
  expect_error(link_lookup("identity", family = "weibul"), "weibul")
  expect_error(link_lookup("identiy", family = "gaussian"), "identiy")
})
