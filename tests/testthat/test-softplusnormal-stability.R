# Issue #31: softplus/inv_softplus and the softplusnormal density/RNG used
# log(exp(x) - 1) / log(exp(x) + 1) which overflow or cancel on valid inputs.

test_that("softplus link is stable for large, tiny and boundary inputs (#31)", {
  # large positive locations: naive log(exp(x) - 1) overflows to NaN/Inf
  expect_equal(softplus(1000), 1000)
  expect_equal(softplus(700), 700)
  # tiny positive values: naive form gives log(0) = -Inf
  expect_equal(softplus(1e-20), log(1e-20))
  expect_equal(softplus(1e-300), log(1e-300))
  # boundary behavior preserved
  expect_equal(softplus(0), -Inf)
  expect_equal(softplus(Inf), Inf)
  # outside the positive domain: NaN (with warning), as before
  expect_warning(expect_true(is.na(softplus(-1))))
})

test_that("inv_softplus is stable for large magnitudes and tiny values (#31)", {
  # naive log(exp(x) + 1) overflows to Inf for large x
  expect_equal(inv_softplus(1000), 1000)
  expect_equal(inv_softplus(700), 700)
  # tiny values: naive form cancels to 0
  expect_equal(inv_softplus(-40), exp(-40))
  expect_equal(inv_softplus(-745), exp(-745))
  # boundary behavior preserved
  expect_equal(inv_softplus(-Inf), 0)
  expect_equal(inv_softplus(Inf), Inf)
  # accuracy where the naive form loses the + 1 entirely
  expect_equal(inv_softplus(20), 20 + exp(-20))
})

test_that("softplus / inv_softplus round trips (#31)", {
  # softplus(inv_softplus(u)) = u for real u
  u <- c(-700, -40, -10, -1, -0.01, 0, 0.01, 1, 10, 40, 700, 1000)
  expect_equal(softplus(inv_softplus(u)), u)
  # inv_softplus(softplus(y)) = y for positive y
  y <- c(1e-10, 1e-4, 0.01, 1, 5, 20, 700, 1e10)
  expect_equal(inv_softplus(softplus(y)), y)
})

test_that("softplus and inv_softplus agree with the naive formulas in the safe range (#31)", {
  x_naive_safe <- c(0.05, 0.1, 0.5, 1, 2, 5, 10, 30)
  expect_equal(softplus(x_naive_safe), log(exp(x_naive_safe) - 1))
  u_naive_safe <- c(-30, -10, -2, -0.5, 0, 0.5, 2, 10, 30)
  expect_equal(inv_softplus(u_naive_safe), log(exp(u_naive_safe) + 1))
})

test_that("dsoftplusnormal finite centered log densities at large locations (#31)", {
  lpdf <- dsoftplusnormal(1000, mu = 1000, sigma = 1, log = TRUE)
  expect_true(is.finite(lpdf))
  expect_equal(lpdf, -log(sqrt(2 * pi)), tolerance = 1e-12)
  # ordinary scale stays finite too
  expect_true(is.finite(dsoftplusnormal(1000, mu = 1000, sigma = 1)))
  # observation far from mu keeps a finite (very negative) log density
  expect_true(is.finite(dsoftplusnormal(1, mu = 1000, sigma = 1, log = TRUE)))
})

test_that("dsoftplusnormal matches the naive formula in the safe range (#31)", {
  x <- c(0.05, 0.1, 0.5, 1, 2, 5, 10, 20, 30)
  mu <- c(-5, 0, 2, 10, 30)
  sigma <- c(0.3, 1, 5)
  for (m in mu) {
    for (s in sigma) {
      naive <- suppressWarnings(
        -(log(s) + 0.5 * log(2 * pi)) +
          x -
          log(exp(x) - 1) +
          -0.5 * ((log(exp(x) - 1) - m) / s)^2
      )
      expect_equal(dsoftplusnormal(x, mu = m, sigma = s, log = TRUE), naive)
    }
  }
})

test_that("dsoftplusnormal integrates to one (#31)", {
  for (par in list(c(0.5, 1), c(2, 2), c(10, 0.5), c(-3, 2))) {
    total <- integrate(
      function(x) dsoftplusnormal(x, mu = par[1], sigma = par[2]),
      lower = 0,
      upper = Inf,
      rel.tol = 1e-10
    )
    expect_equal(total$value, 1, tolerance = 1e-8)
  }
})

test_that("rsoftplusnormal supports large locations (#31)", {
  set.seed(20240930)
  draws <- rsoftplusnormal(2000, mu = 1000, sigma = 1)
  expect_true(all(is.finite(draws)))
  expect_equal(length(draws), 2000)
  # median of inv_softplus(Z) is inv_softplus(median(Z)) = inv_softplus(1000)
  expect_equal(median(draws), inv_softplus(1000), tolerance = 0.1)
})

test_that("injected Stan functions agree with the R density and RNG (#31)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("pkgload")
  scode <- softplusnormal()$stanvars[[1]]$scode
  # guard against regressions to the unstable formulas
  expect_false(grepl("log\\(exp\\(", scode))
  expect_true(grepl("log1m_exp", scode, fixed = TRUE))
  expect_true(grepl("log1p_exp", scode, fixed = TRUE))

  model_code <- paste("functions {", scode, "}", "model { target += 0; }")
  compiled <- tryCatch(
    {
      sm <- rstan::stan_model(model_code = model_code, verbose = FALSE)
      rstan::expose_stan_functions(sm)
    },
    error = function(e) NULL
  )
  skip_if(is.null(compiled), "Stan functions could not be compiled")
  # expose_stan_functions places the R wrappers into the global environment
  stan_lpdf <- get("softplusnormal_lpdf", envir = globalenv())
  stan_rng <- get("softplusnormal_rng", envir = globalenv())

  y <- c(1e-6, 0.001, 0.1, 1, 10, 100, 1000)
  mu <- c(-5, 0, 5, 1000)
  sigma <- c(0.5, 2)
  for (m in mu) {
    for (s in sigma) {
      r_lpdf <- dsoftplusnormal(y, mu = m, sigma = s, log = TRUE)
      stan_values <- vapply(y, stan_lpdf, numeric(1), mu = m, sigma = s)
      expect_equal(stan_values, r_lpdf, tolerance = 1e-10)
    }
  }
  # Stan RNG stays finite at large locations
  expect_true(all(is.finite(
    vapply(1:20, function(i) stan_rng(1000, 1), numeric(1))
  )))
})
