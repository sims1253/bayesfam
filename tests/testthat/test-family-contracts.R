# Contracts every custom family must meet (#28); fixtures in helper-families.R.
families <- sub(
  "^log_lik_",
  "",
  ls(asNamespace("bayesfam"), pattern = "^log_lik_")
)

# posterior_epred callbacks that do not return the mean, and the condition
# they signal instead
epred_exceptions <- c(
  cauchitnormal = "error", # no closed-form mean
  cloglognormal = "error", # no closed-form mean
  generalized_gamma = "error", # not implemented
  gompertz = "error", # not implemented
  softplusnormal = "error", # no closed-form mean
  logitnormal = "warning", # returns the median, the mean is added in #46
  symlognormal = "warning" # returns nothing, the mean is added in #46
)

test_that("every custom family has a fixture", {
  expect_setequal(names(family_fixtures), families)
  for (f in families) {
    fx <- family_fixture(f)
    expect_named(fx$pars, fx$fam$dpars)
  }
})

# A brmsprep as the callbacks receive it: S draws of an S x N mu matrix and
# length-S vectors for the other dpars, cycling through the fixture values.
S <- 7
mock_prep <- function(fx) {
  N <- length(fx$y)
  dpars <- lapply(fx$pars, rep_len, S)
  dpars$mu <- matrix(rep_len(fx$pars$mu, S * N), S, N)
  data <- c(list(Y = fx$y), lapply(fx$vreal, rep, N))
  structure(
    list(ndraws = S, nobs = N, dpars = dpars, data = data),
    class = "brmsprep"
  )
}

# Arguments of the reference lpdf for observation i of a prep
lpdf_args <- function(fx, prep, i) {
  dpars <- lapply(fx$fam$dpars, function(dp) brms::get_dpar(prep, dp, i = i))
  c(dpars, lapply(prep$data[names(fx$vreal)], `[`, i))
}

# Integral of g(y) times the density over the support
expectation <- function(fx, args, g = function(y) 1) {
  integrand <- function(y) {
    g(y) * exp(do.call(fx$lpdf, c(list(y), unname(args))))
  }
  integrate(integrand, fx$support[1], fx$support[2], rel.tol = 1e-6)$value
}

for (f in families) {
  test_that(paste(f, "density integrates to 1"), {
    fx <- family_fixture(f)
    for (j in 1:2) {
      args <- c(lapply(unname(fx$pars), `[`, j), fx$vreal)
      expect_equal(expectation(fx, args), 1, tolerance = 1e-4)
    }
  })

  test_that(paste(f, "callbacks agree with the density and brms shapes"), {
    fx <- family_fixture(f)
    prep <- mock_prep(fx)
    N <- prep$nobs
    for (i in seq_len(N)) {
      args <- lpdf_args(fx, prep, i)
      expect_equal(
        fx$fam$log_lik(i, prep),
        do.call(fx$lpdf, c(list(fx$y[i]), unname(args)))
      )
      # posterior_predict draws S values from the family's RNG
      set.seed(i)
      yrep <- fx$fam$posterior_predict(i, prep)
      set.seed(i)
      expect_equal(yrep, do.call(fx$rng, c(list(S), unname(args))))
    }

    if (f %in% names(epred_exceptions)) {
      expect_condition(
        fx$fam$posterior_epred(prep),
        class = epred_exceptions[[f]]
      )
      return()
    }
    epred <- fx$fam$posterior_epred(prep)
    expect_equal(dim(epred), c(S, N))
    expect_true(all(is.finite(epred)))
    # epred is the mean: compare the first draw of observation 1 to E[y]
    args <- lapply(lpdf_args(fx, prep, 1), `[`, 1)
    expect_equal(epred[1, 1], expectation(fx, args, identity), tolerance = 1e-4)
  })
}

test_that("Stan lpdfs match the R log densities", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  # generalized_normal redefines symlognormal's `int sign(real)` (removed in #49)
  fxs <- lapply(setdiff(families, "generalized_normal"), family_fixture)
  scode <- unlist(lapply(fxs, function(fx) {
    lapply(fx$fam$stanvars, `[[`, "scode")
  }))
  model <- rstan::stanc(
    model_code = paste("functions {", paste(scode, collapse = "\n"), "}")
  )
  stan <- new.env()
  rstan::expose_stan_functions(model, env = stan)
  for (fx in fxs) {
    for (j in 1:2) {
      args <- unname(c(lapply(fx$pars, `[`, j), fx$vreal))
      for (y in fx$y) {
        expect_equal(
          do.call(stan[[paste0(fx$fam$name, "_lpdf")]], c(list(y), args)),
          do.call(fx$lpdf, c(list(y), args)),
          label = paste0(fx$fam$name, "_lpdf(", y, ")")
        )
      }
    }
  }
})
