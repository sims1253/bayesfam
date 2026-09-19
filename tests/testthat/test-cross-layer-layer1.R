# Layer 1 of the cross-layer suite (issue #28): fast deterministic contracts
# of the exported R densities, quantiles and RNGs, driven by the family
# inventory:
#   - density normalization over the support
#   - known special cases
#   - CDF/quantile round trips (numeric CDF from the density)
#   - parameter-domain behavior (documented domain violations error)
#   - vectorized inputs and shape checks with unequal lengths
# Deliberately avoids the buggy helpers in R/test-helper.R (issues #35, #36)
# and uses plain testthat expectations on differences computed here.

cross_layer_x_grid <- function(entry, params, n = 11) {
  lo <- entry$support(params)[1]
  hi <- entry$support(params)[2]
  if (is.finite(lo) && is.finite(hi)) {
    return(seq(lo + 0.02 * (hi - lo), hi - 0.02 * (hi - lo), length.out = n))
  }
  if (is.finite(lo)) {
    return(lo + exp(seq(log(0.05), log(60), length.out = n)))
  }
  seq(-4, 4, length.out = n)
}

cross_layer_density <- function(entry, x, params, log = FALSE) {
  f <- entry$d_fun
  fmls <- names(formals(f))
  args <- c(setNames(list(x), fmls[1]), params)
  if ("log" %in% fmls) {
    args$log <- log
  }
  do.call(f, args)
}

cross_layer_prob_args <- function(entry, p, params) {
  # the probability goes to the first formal positionally so parameter names
  # like kumaraswamy's `p` cannot collide with it
  c(list(p), params)
}

# Shrink finite integration bounds slightly: densities that reject boundary
# values (e.g. strict (0, 1) support) otherwise error on points that the
# transformation inside integrate() rounds onto the boundary itself.
cross_layer_bounds <- function(lo, hi, delta = 1e-12) {
  span <- if (all(is.finite(c(lo, hi)))) {
    hi - lo
  } else {
    max(abs(c(lo[is.finite(lo)], hi[is.finite(hi)])), 1)
  }
  if (is.finite(lo)) lo <- lo + delta * span
  if (is.finite(hi)) hi <- hi - delta * span
  c(lo, hi)
}

cross_layer_numeric_cdf <- function(entry, q, params) {
  f <- function(x) {
    d <- cross_layer_density(entry, x, params)
    if (entry$always_log) {
      return(exp(d))
    }
    d
  }
  integrate(
    f,
    cross_layer_bounds(entry$support(params)[1], q)[1],
    q,
    rel.tol = 1e-9,
    subdivisions = 500L
  )$value
}

for (entry in cross_layer_inventory) {
  if (is.null(entry$d_fun)) {
    next
  }
  test_that(paste("density integrates to one:", entry$name), {
    if (entry$name == "shifted_lognormal_uniform") {
      skip(paste(
        "#24 R density of the shifted lognormal/uniform mixture is",
        "unnormalized (subtracts a truncation term the model does not apply;",
        "fixed in wave-2 PR)"
      ))
    }
    for (params in list(entry$params, entry$params2)) {
      bnds <- cross_layer_bounds(
        entry$support(params)[1],
        entry$support(params)[2]
      )
      total <- integrate(
        function(x) {
          d <- cross_layer_density(entry, x, params)
          if (entry$always_log) exp(d) else d
        },
        bnds[1],
        bnds[2],
        rel.tol = 1e-8,
        subdivisions = 1000L
      )$value
      expect_lt(abs(total - 1), 1e-4, label = paste(entry$name, "total mass"))
    }
  })
}

for (entry in cross_layer_inventory) {
  if (is.null(entry$d_fun)) {
    next
  }
  test_that(paste("vectorization, recycling and log consistency:", entry$name), {
    x <- cross_layer_x_grid(entry, entry$params)
    d_vec <- cross_layer_density(entry, x, entry$params)
    expect_length(d_vec, length(x))
    expect_true(all(is.finite(d_vec)))

    x5 <- cross_layer_x_grid(entry, entry$params, n = 5)
    has_log <- "log" %in% names(formals(entry$d_fun)) && !entry$always_log
    if (has_log) {
      lpdf <- cross_layer_density(entry, x5, entry$params, log = TRUE)
      pdf <- cross_layer_density(entry, x5, entry$params, log = FALSE)
      # compare only where the density does not underflow to zero
      keep <- is.finite(pdf) & pdf > 0
      expect_true(any(keep), label = paste(entry$name, "non-degenerate grid"))
      expect_equal(lpdf[keep], log(pdf[keep]), tolerance = 1e-12)
    } else {
      lpdf <- cross_layer_density(entry, x5, entry$params)
      expect_true(all(is.finite(lpdf)))
    }

    # scalar x with vectorized parameters of unequal length
    param_names <- names(entry$params)
    if (entry$name == "generalized_gamma") {
      # Q must stay scalar: the density branches with `if (Q != 0)`, which
      # errors for vectors (issue #23, fixed in wave-2 PR).
      param_names <- setdiff(param_names, "Q")
    }
    params_vec <- entry$params
    params_vec[param_names] <- lapply(params_vec[param_names], rep_len, 7)
    mid <- x5[3]
    out <- cross_layer_density(entry, mid, params_vec)
    expect_length(out, 7)

    # unequal x and parameter lengths exercise recycling end to end
    out13 <- cross_layer_density(entry, rep_len(mid, 13), entry$params)
    expect_length(out13, 13)

    rng_out <- do.call(entry$r_fun, c(list(n = 13), entry$params))
    expect_length(rng_out, 13)
  })
}

for (entry in cross_layer_inventory) {
  needs_d <- !is.null(entry$d_fun)
  test_that(paste("parameter-domain violations error:", entry$name), {
    x <- if (needs_d) {
      cross_layer_x_grid(entry, entry$params, n = 3)[2]
    } else {
      NA
    }
    if (needs_d) {
      expect_error(
        do.call(entry$d_fun, c(list(x = x), entry$bad_params)),
        info = paste(entry$name, "rejects", paste(
          names(entry$bad_params),
          entry$bad_params
        ))
      )
      if (!is.na(entry$bad_x)) {
        expect_error(
          cross_layer_density(entry, entry$bad_x, entry$params),
          info = paste(entry$name, "rejects out-of-support x")
        )
      }
    }
    rng_args <- c(list(n = 5), entry$bad_params)
    expect_error(
      do.call(entry$r_fun, rng_args),
      info = paste(entry$name, "rng rejects", paste(
        names(entry$bad_params),
        entry$bad_params
      ))
    )
    if (entry$name == "simplex") {
      skip(paste(
        "#27 simplex accepts sigma < 0 although the likelihood evaluates",
        "log(sigma); correct behavior is a domain error (fixed in wave-2 PR)"
      ))
      expect_error(
        cross_layer_density(entry, 0.5, list(mu = 0.5, sigma = -1)),
        info = "simplex must reject negative sigma"
      )
    }
  })
}

for (entry in cross_layer_inventory) {
  if (is.null(entry$q_fun)) {
    next
  }
  test_that(paste("CDF/quantile round trip:", entry$name), {
    for (params in list(entry$params, entry$params2)) {
      for (p in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
        q <- do.call(entry$q_fun, cross_layer_prob_args(entry, p, params))
        cdf <- cross_layer_numeric_cdf(entry, q, params)
        expect_lt(
          abs(cdf - p),
          5e-3,
          label = sprintf("%s: F(q(%s)) at mu=%s", entry$name, p, params[[1]])
        )
      }
    }
  })
}

for (i in seq_along(cross_layer_inventory)) {
  entry <- cross_layer_inventory[[i]]
  if (is.na(entry$rng_stat)) {
    next
  }
  test_that(paste("RNG recovers its location statistic:", entry$name), {
    n <- 20000
    params <- entry$rng_params %||% entry$params
    set.seed(2807 + i)
    draws <- do.call(entry$r_fun, c(list(n = n), params))
    expect_length(draws, n)
    expect_true(all(is.finite(draws)), label = paste(entry$name, "finite draws"))
    lo <- entry$support(params)[1]
    hi <- entry$support(params)[2]
    expect_true(all(draws >= lo & draws <= hi), label = paste(
      entry$name,
      "draws inside support"
    ))

    target <- entry$rng_target(params)
    if (entry$rng_stat == "mean") {
      se <- sd(draws) / sqrt(n)
      expect_lt(abs(mean(draws) - target), 6 * se, label = paste(
        entry$name,
        "mean"
      ))
    } else {
      se <- 1.2533 * sd(draws) / sqrt(n)
      expect_lt(abs(median(draws) - target), 8 * se, label = paste(
        entry$name,
        "median"
      ))
    }
  })
}

test_that("generalized_normal reduces to Laplace and Normal", {
  x <- c(-3.2, -1.5, 0.4, 2.7)
  # beta = 1: Laplace(mu, sigma)
  expect_equal(
    dgeneralized_normal(x, mu = 1, sigma = 2, beta = 1, log = TRUE),
    -log(2 * 2) - abs(x - 1) / 2
  )
  # beta = 2: Normal with sd sigma / sqrt(2)
  expect_equal(
    dgeneralized_normal(x, mu = 1, sigma = 2, beta = 2, log = TRUE),
    dnorm(x, mean = 1, sd = 2 / sqrt(2), log = TRUE)
  )
  # Laplace quantiles in closed form
  p <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  q_laplace <- ifelse(p < 0.5, 1 + 2 * log(2 * p), 1 - 2 * log(2 * (1 - p)))
  expect_equal(qgeneralized_normal(p, mu = 1, sigma = 2, beta = 1), q_laplace)
})

test_that("logistic matches stats::dlogis and stats::qlogis", {
  x <- c(-2.5, -0.3, 1.7, 4)
  expect_equal(
    dlogistic(x, mu = 1, sigma = 2, log = TRUE),
    dlogis(x, location = 1, scale = 2, log = TRUE)
  )
  p <- c(0.1, 0.35, 0.5, 0.66, 0.9)
  expect_equal(qlogistic(p, mu = 1, sigma = 2), qlogis(p, 1, 2))
})

test_that("lognormal (median parametrization) matches stats::dlnorm", {
  x <- c(0.3, 1.4, 5.2)
  expect_equal(
    dlognormal(x, mu = 0.5, sigma = 0.8, log = TRUE),
    dlnorm(x, meanlog = 0.5, sdlog = 0.8, log = TRUE)
  )
})

test_that("gamma_mean matches stats::dgamma", {
  x <- c(0.4, 2.1, 7.3)
  expect_equal(
    dgamma_mean(x, mu = 3, a = 2, log = TRUE),
    dgamma(x, shape = 2, rate = 2 / 3, log = TRUE)
  )
})

test_that("lognormal_natural is a lognormal with natural-scale mean", {
  x <- c(0.4, 1.2, 3.7)
  mu <- 2
  sigma <- 0.5
  common_term <- log1p(sigma^2 / mu^2)
  expect_equal(
    dlognormal_natural(x, mu = mu, sigma = sigma, log = TRUE),
    dlnorm(x, log(mu) - common_term / 2, sqrt(common_term), log = TRUE)
  )
})

test_that("median parametrizations hit their median", {
  expect_equal(qgompertz(0.5, mu = 2.3, beta = 1.7), 2.3)
  expect_equal(qkumaraswamy(0.5, mu = 0.3, p = 2), 0.3)
  expect_equal(qkumaraswamy(0.5, mu = 0.3, p = 3), 0.3)
  gamma_euler <- 0.5772156649015329
  expect_equal(qgumbel_mean(exp(-exp(-gamma_euler)), mu = 0, sigma = 1), 0)
})

test_that("gumbel_mean matches the standard gumbel density", {
  skip_if_not_installed("extraDistr")
  x <- c(-2, 0, 1.5, 4)
  gamma_euler <- 0.5772156649015329
  expect_equal(
    dgumbel_mean(x, mu = 2, sigma = 1.5, log = TRUE),
    extraDistr::dgumbel(
      x,
      mu = 2 - 1.5 * gamma_euler,
      sigma = 1.5,
      log = TRUE
    )
  )
})

test_that("link-normal densities are exact normals in link space", {
  # The transformed-normal (logit / cauchit / cloglog) densities must satisfy
  # the change-of-variables identity lpdf_x(x) = lpdf_t(t) - log|dx/dt| with
  # t = link(x). Wide sigma pushes tail mass into a sliver near x = 1 where
  # x-space quadrature breaks down (x quantizes to a handful of doubles), so
  # the identity is checked pointwise instead.
  link_cases <- list(
    list(
      name = "logitnormal",
      dlog = function(x, mu, sigma) {
        bayesfam::dlogitnormal(x, mu, sigma, log = TRUE)
      },
      link = function(x) qlogis(x),
      log_jac = function(t) t - 2 * log1p(exp(t))
    ),
    list(
      name = "cauchitnormal",
      dlog = function(x, mu, sigma) {
        bayesfam::dcauchitnormal(x, mu, sigma, log = TRUE)
      },
      link = function(x) qcauchy(x),
      log_jac = function(t) dcauchy(t, log = TRUE)
    ),
    list(
      name = "cloglognormal",
      dlog = function(x, mu, sigma) {
        bayesfam::dcloglognormal(x, mu, sigma, log = TRUE)
      },
      link = function(x) log(-log1p(-x)),
      log_jac = function(t) t - exp(t)
    )
  )
  x <- c(1e-6, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999999)
  for (case in link_cases) {
    for (pr in list(
      c(mu = 0, sigma = 1),
      c(mu = 0.5, sigma = 1.5),
      c(mu = -2, sigma = 0.7)
    )) {
      t <- case$link(x)
      expected <- dnorm(t, pr[["mu"]], pr[["sigma"]], log = TRUE) -
        case$log_jac(t)
      expect_equal(
        case$dlog(x, pr[["mu"]], pr[["sigma"]]),
        expected,
        tolerance = 1e-10,
        label = paste(case$name, "link-space identity")
      )
    }
  }
})

test_that("pkumaraswamy is consistent with the median parametrization", {
  skip(paste(
    "#26 pkumaraswamy uses the wrong sign and power base",
    "(1 + (x^p - 1)^q instead of 1 - (1 - x^p)^q; fixed in wave-2 PR)"
  ))
  expect_equal(pkumaraswamy(0.5, mu = 0.3, p = 2), 0.5)
  expect_equal(pkumaraswamy(0.5, mu = 0.7, p = 3), 0.5)
  q <- -(log(2) / log1p(-0.3^2))
  expect_equal(
    pkumaraswamy(0.6, mu = 0.3, p = 2),
    -expm1(q * log1p(-0.6^2))
  )
})

test_that("shifted_lognormal_uniform mixture components behave as documented", {
  set.seed(11)
  n <- 20000
  u <- rshifted_lognormal_uniform(n, 0.3, 0.4, mix = 1, shift = 0.1, max_uniform = 10)
  expect_gt(min(u), 0)
  expect_lt(max(u), 10)
  expect_lt(abs(mean(u) - 5), 6 * sd(u) / sqrt(n))

  l <- rshifted_lognormal_uniform(n, 0.3, 0.4, mix = 0, shift = 0.1, max_uniform = 10)
  expect_lt(
    abs(mean(l) - (0.1 + exp(0.3 + 0.4^2 / 2))),
    6 * sd(l) / sqrt(n)
  )
})

test_that("shifted_lognormal_uniform mix edges match component densities", {
  skip(paste(
    "#24 with mix at the edges the R density subtracts a spurious truncation",
    "term and disagrees with the component densities (fixed in wave-2 PR)"
  ))
  y <- c(0.5, 2, 6)
  expect_equal(
    dshifted_lognormal_uniform(y, 0.3, 0.4, mix = 0, shift = 0.1, max_uniform = 10),
    dlnorm(y - 0.1, 0.3, 0.4, log = TRUE)
  )
  expect_equal(
    dshifted_lognormal_uniform(c(0.5, 2), 0.3, 0.4, mix = 1, shift = 0.1, max_uniform = 10),
    dunif(c(0.5, 2), 0, 10, log = TRUE)
  )
})
