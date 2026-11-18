# Contracts every custom family must meet (#28); fixtures in helper-families.R.
families <- sub(
  "^log_lik_",
  "",
  ls(asNamespace("bayesfam"), pattern = "^log_lik_")
)

test_that("every custom family has a fixture", {
  expect_identical(sort(names(family_fixtures)), sort(families))
  for (f in families) {
    fx <- family_fixture(f)
    expect_named(fx$pars, fx$fam$dpars)
  }
})

# The range of each default inverse link must lie inside the declared bounds
# of its distributional parameter: otherwise the sampler can reach parameter
# values the family's own density rejects.
test_that("default inverse links map into each dpar's declared bounds", {
  # data-dependent bounds ("min_Y") count as unbounded here
  bound <- function(b, default) {
    b <- suppressWarnings(as.numeric(unlist(b)))
    replace(b, is.na(b), default)
  }
  # the linear predictor is unbounded, so probe far into both tails
  eta <- seq(-37, 37, by = 0.1)
  for (f in families) {
    fam <- family_fixture(f)$fam
    lb <- bound(fam$lb, -Inf)
    ub <- bound(fam$ub, Inf)
    for (i in seq_along(fam$dpars)) {
      dp <- fam$dpars[i]
      link <- if (dp == "mu") fam$link else fam[[paste0("link_", dp)]]
      # brms' own response functions, as used in the fitted model
      range <- brms:::inv_link(eta, link)
      expect_true(
        all(range >= lb[i] & range <= ub[i]),
        label = sprintf(
          "%s: inverse '%s' link of '%s' reaches [%s, %s] outside bounds [%s, %s]",
          f, link, dp, format(min(range)), format(max(range)),
          format(lb[i]), format(ub[i])
        )
      )
    }
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

    epred <- fx$fam$posterior_epred(prep)
    expect_equal(dim(epred), c(S, N))
    expect_true(all(is.finite(epred)))
    # epred is the mean: compare the first draw of observation 1 to E[y]
    args <- lapply(lpdf_args(fx, prep, 1), `[`, 1)
    expect_equal(epred[1, 1], expectation(fx, args, identity), tolerance = 1e-4)
  })
}

# Means without a textbook closed form, checked tighter and over more
# parameter sets than the fixtures. The unit-interval sets keep the mass
# that rounds to 1 in double precision negligible, so the reference
# integral over y can resolve them.
epred_cases <- list(
  gompertz = list(
    c(mu = 1, beta = 0.5), c(mu = 2, beta = 1), c(mu = 0.3, beta = 4),
    c(mu = 5, beta = 2)
  ),
  generalized_gamma = list(
    c(mu = 0, sigma = 0.5, Q = -0.5), c(mu = 0.5, sigma = 1, Q = 1.5),
    c(mu = 1, sigma = 0.3, Q = 0), c(mu = -1, sigma = 0.8, Q = 0.01),
    c(mu = 0, sigma = 2, Q = -0.4)
  ),
  cauchitnormal = list(
    c(mu = -0.5, sigma = 0.5), c(mu = 0.5, sigma = 1), c(mu = 2, sigma = 3),
    c(mu = 0, sigma = 10)
  ),
  cloglognormal = list(
    c(mu = -1, sigma = 0.3), c(mu = -0.5, sigma = 0.5),
    c(mu = -4, sigma = 1.2), c(mu = -6, sigma = 1.5)
  ),
  softplusnormal = list(
    c(mu = 0.5, sigma = 0.5), c(mu = 1, sigma = 1), c(mu = -2, sigma = 3),
    c(mu = 5, sigma = 0.1)
  )
)

for (f in names(epred_cases)) {
  test_that(paste(f, "posterior_epred is the mean of the density"), {
    fx <- family_fixture(f)
    support <- fx$support
    for (p in epred_cases[[f]]) {
      N <- 3
      dpars <- lapply(as.list(p), rep, S)
      dpars$mu <- matrix(p[["mu"]], S, N)
      prep <- structure(
        list(ndraws = S, nobs = N, dpars = dpars),
        class = "brmsprep"
      )
      epred <- fx$fam$posterior_epred(prep)
      expect_equal(dim(epred), c(S, N))
      integrand <- function(y) {
        inside <- y > support[1] & y < support[2]
        out <- 0 * y
        out[inside] <- y[inside] *
          exp(do.call(fx$lpdf, c(list(y[inside]), unname(as.list(p)))))
        out
      }
      # split near the median, where the mass sits
      set.seed(1)
      split <- median(do.call(fx$rng, c(list(1e4), unname(as.list(p)))))
      mean <- integrate(integrand, support[1], split,
        rel.tol = 1e-10, subdivisions = 2000L
      )$value + integrate(integrand, split, support[2],
        rel.tol = 1e-10, subdivisions = 2000L
      )$value
      expect_equal(
        epred, matrix(mean, S, N),
        tolerance = 1e-6,
        label = paste(f, paste(names(p), p, sep = "=", collapse = ", "))
      )
    }
  })
}

test_that("generalized_gamma posterior_epred is Inf where the mean diverges", {
  # E[y] exists only for Q > -1/sigma; the third draw is just inside
  prep <- structure(list(ndraws = 3, nobs = 2, dpars = list(
    mu = matrix(0, 3, 2), sigma = c(1, 2, 0.5), Q = c(-1, -0.6, -0.5)
  )), class = "brmsprep")
  expect_equal(
    generalized_gamma()$posterior_epred(prep),
    matrix(c(Inf, Inf, 4 / 3), 3, 2)
  )
})

# The Stan code of family `f`, without includes of brms' own function chunks:
# they declare overloads, which expose_stan_functions() cannot export to R.
# The only one is fun_logm1.stan, which lomax includes for expp1, the inverse
# of its default logm1 link. The model block calls expp1; the lpdf and RNG do
# not, so they expose complete without it. Stop if that stops being true.
family_scode <- function(f) {
  scode <- unlist(lapply(family_fixture(f)$fam$stanvars, `[[`, "scode"))
  include <- grepl("^#include", scode)
  own <- paste(scode[!include], collapse = "\n")
  stopifnot(
    all(scode[include] == "#include 'fun_logm1.stan'"),
    !any(include) || !grepl("\\b(expp1|logm1)\\s*\\(", own)
  )
  own
}

# The name of the one RNG in the Stan code of family `f`
stan_rng_name <- function(f) {
  names <- regmatches(
    family_scode(f),
    gregexpr("real\\s+\\w+_rng(?=\\s*\\()", family_scode(f), perl = TRUE)
  )[[1]]
  sub("^real\\s+", "", names)
}

test_that("every family has an RNG in its Stan code", {
  for (f in families) {
    expect_equal(length(stan_rng_name(f)), 1, label = paste(f, "Stan RNG count"))
  }
})

# All family functions in one compiled model, exposed to R once. Each RNG
# gets a wrapper that returns n draws, so one call with one seed gives an
# independent sample.
stan_functions <- local({
  stan <- NULL
  function() {
    if (is.null(stan)) {
      wrappers <- vapply(families, function(f) {
        rng <- stan_rng_name(f)
        k <- length(family_fixture(f)$pars) + length(family_fixture(f)$vreal)
        args <- paste0("a", seq_len(k))
        sprintf(
          "vector draws_%s(int n, %s) {
            vector[n] y;
            for (i in 1:n) y[i] = %s(%s);
            return y;
          }",
          rng, paste("real", args, collapse = ", "), rng,
          paste(args, collapse = ", ")
        )
      }, character(1))
      model <- rstan::stanc(model_code = paste(
        "functions {",
        paste(vapply(families, family_scode, character(1)), collapse = "\n"),
        paste(wrappers, collapse = "\n"),
        "}"
      ))
      stan <<- new.env()
      rstan::expose_stan_functions(model, env = stan)
    }
    stan
  }
})

test_that("Stan lpdfs match the R log densities", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  fxs <- lapply(families, family_fixture)
  stan <- stan_functions()
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

test_that("Stan RNGs draw from the same distribution as the R RNGs", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  stan <- stan_functions()
  n <- 20000
  ns <- asNamespace("bayesfam")
  for (f in families) {
    fx <- family_fixture(f)
    cdf <- paste0("p", fx$fam$name)
    for (j in 1:2) {
      args <- unname(c(lapply(fx$pars, `[`, j), fx$vreal))
      draws <- do.call(
        stan[[paste0("draws_", stan_rng_name(f))]],
        c(list(n), args, list(base_rng__ = rstan::get_rng(j)))
      )
      expect_length(draws, n)
      # against the R CDF where there is one, else against the R RNG
      test <- if (exists(cdf, ns, inherits = FALSE)) {
        do.call(ks.test, c(list(draws, get(cdf, ns)), args))
      } else {
        set.seed(j)
        ks.test(draws, do.call(fx$rng, c(list(n), args)))
      }
      expect_gt(
        test$p.value, 0.001,
        label = sprintf("%s KS p-value (parameter set %d)", f, j)
      )
    }
  }
})

test_that("the unit Lindley Stan RNG draws from the density for small mu", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  # exp(-1 / mu) underflows to zero for mu below about 0.00134, where the
  # inverse CDF gave NaN
  stan <- stan_functions()
  mu <- 0.001
  n <- 1e5
  y <- stan$draws_unit_lindley_rng(n, mu, base_rng__ = rstan::get_rng(1L))
  expect_true(all(is.finite(y) & y > 0 & y < 1))
  # the mean of this parameterisation is mu; allow five standard errors
  expect_lt(abs(mean(y) - mu), 5 * sd(y) / sqrt(n))
})
