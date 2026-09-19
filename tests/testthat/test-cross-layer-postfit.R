# Layer 3 of the cross-layer suite (issue #28): post-fit callback smoke tests
# for log_lik / posterior_predict / posterior_epred of every custom family.
#
# Full brms MCMC fits are too slow for this suite, so the callbacks are driven
# with mock prep objects (make_cross-layer style, class "brmsprep") built from
# the inventory fixtures: S = 7 draws and N = 3 observations (unequal on
# purpose), per-observation mu vectors (distributional predictors), and
# observation-level vreal data for shifted_lognormal_uniform.
#
# Callbacks whose math is broken by bugs fixed on sibling branches (#22, #23,
# #24, #25) are exercised against the correct references but skipped with the
# issue number, so they light up when the wave-2 fixes land here. Callbacks
# documented as unsupported must stop / warn as documented and then skip.

cross_layer_n_draws <- 7

for (name in cross_layer_custom_families) {
  entry <- cross_layer_inventory[[name]]
  fam <- entry$family_fun()
  pf <- entry$postfit
  nobs <- length(pf$y)
  prep <- postfit_prep(entry, ndraws = cross_layer_n_draws)

  test_that(paste("log_lik callback:", name), {
    if (entry$ll3$status == "skip") {
      skip(paste0(entry$ll3$reason, "; reference kept below for the fixed code"))
    }
    for (i in seq_len(nobs)) {
      out <- fam$log_lik(i, prep)
      expect_equal(length(out), cross_layer_n_draws, label = paste(name, "length i =", i))
      expect_true(all(is.finite(out)), label = paste(name, "finite"))
      args <- postfit_draw_args(entry, prep, i)
      ref <- do.call(entry$d_fun, c(list(x = pf$y[i], log = TRUE), args))
      expect_equal(
        out,
        ref,
        tolerance = 1e-12,
        label = paste(name, "log_lik equals R density at draws, i =", i)
      )
    }
  })

  test_that(paste("posterior_predict callback:", name), {
    if (entry$pp3$status == "skip") {
      skip(paste0(entry$pp3$reason, "; reference kept below for the fixed code"))
    }
    for (i in seq_len(nobs)) {
      set.seed(4200 + i)
      ref <- do.call(
        entry$r_fun,
        c(list(n = cross_layer_n_draws), postfit_draw_args(entry, prep, i))
      )
      set.seed(4200 + i)
      out <- fam$posterior_predict(i, prep)
      expect_equal(length(out), cross_layer_n_draws, label = paste(name, "length i =", i))
      expect_equal(
        out,
        ref,
        tolerance = 1e-12,
        label = paste(name, "posterior_predict equals family RNG, i =", i)
      )
      lo_hi_params <- pf$params
      if (name == "shifted_inv_gaussian") {
        names(lo_hi_params)[names(lo_hi_params) == "ndt"] <- "shift"
      }
      lo <- entry$support(lo_hi_params)[1]
      hi <- entry$support(lo_hi_params)[2]
      expect_true(all(out >= lo & out <= hi), label = paste(
        name,
        "draws in support, i =",
        i
      ))
    }
  })

  test_that(paste("posterior_epred callback:", name), {
    if (entry$ep3$status == "skip") {
      skip(paste0(entry$ep3$reason, "; reference kept below for the fixed code"))
    }
    if (entry$ep3$status == "unsupported") {
      if (name == "symlognormal") {
        warns <- capture_warnings(epred <- fam$posterior_epred(prep))
        expect_gt(
          length(warns),
          0,
          label = paste(name, "posterior_epred must warn (documented)")
        )
      } else {
        expect_error(
          fam$posterior_epred(prep),
          info = paste(name, ":", entry$ep3$reason)
        )
      }
      skip(paste(
        "posterior_epred documented unsupported for",
        name,
        "-",
        entry$ep3$reason
      ))
    }
    if (entry$ep3$status == "warning") {
      warns <- capture_warnings(epred <- fam$posterior_epred(prep))
      expect_gt(length(warns), 0, label = paste(name, "documented warning"))
    } else {
      epred <- fam$posterior_epred(prep)
    }
    expect_equal(
      dim(epred),
      c(cross_layer_n_draws, nobs),
      label = paste(name, "epred is S x N")
    )
    ref <- if (is.null(entry$epred_ref)) prep$dpars$mu else entry$epred_ref(prep)
    expect_equal(epred, ref, label = paste(name, "epred values"))
  })
}

test_that("log_lik callback: shifted_lognormal_uniform plumbing", {
  # The values the callback produces are wrong until #24 (wave-2 PR) lands,
  # but the plumbing must return one finite log-likelihood per posterior draw.
  entry <- cross_layer_inventory[["shifted_lognormal_uniform"]]
  fam <- entry$family_fun()
  prep <- postfit_prep(entry, ndraws = cross_layer_n_draws)
  for (i in seq_along(entry$postfit$y)) {
    out <- fam$log_lik(i, prep)
    expect_length(out, cross_layer_n_draws)
    expect_true(all(is.finite(out)))
  }
})

test_that("posterior_epred callback: shifted_lognormal_uniform N = 1", {
  # With a single observation the #25 recycling bug cannot trigger, so the
  # closed-form mean can be checked directly.
  fam <- shifted_lognormal_uniform()
  S <- 5
  mu <- matrix(c(0.2, 0.25, 0.3, 0.35, 0.4), nrow = S, ncol = 1)
  mix <- matrix(c(0.1, 0.15, 0.2, 0.25, 0.3), nrow = S, ncol = 1)
  prep <- make_cross_layer_prep(
    fam,
    y = 2,
    params = list(mu = mu, sigma = 0.3, mix = mix, shiftprop = 0.5),
    ndraws = S,
    vreal1 = 0.8,
    vreal2 = 10
  )
  epred <- fam$posterior_epred(prep)
  expect_equal(dim(epred), c(S, 1))
  shift <- 0.5 * 0.8
  expected <- mix * 0.5 * 10 +
    (1 - mix) * (shift + exp(mu + 0.3^2 / 2))
  expect_equal(as.vector(epred), as.vector(expected))
})

test_that("posterior_epred callback: shifted_lognormal_uniform recycling", {
  skip(paste(
    "#25 posterior_epred multiplies S x N parameter matrices with length-N",
    "vreal data, recycling observation bounds across draws whenever",
    "S != N; rows must be identical when parameters do not vary by draw",
    "(fixed in wave-2 PR)"
  ))
  fam <- shifted_lognormal_uniform()
  S <- 4
  N <- 2
  prep <- make_cross_layer_prep(
    fam,
    y = c(1, 3),
    params = list(mu = 0.3, sigma = 0.3, mix = 0.2, shiftprop = 0.5),
    ndraws = S,
    vreal1 = c(1, 10),
    vreal2 = c(10, 20)
  )
  epred <- fam$posterior_epred(prep)
  expect_equal(dim(epred), c(S, N))
  expected <- matrix(
    c(0.2 * 5 + 0.8 * (0.5 + exp(0.3 + 0.045)), 0.2 * 10 + 0.8 * (5 + exp(0.3 + 0.045))),
    nrow = S,
    ncol = N,
    byrow = TRUE
  )
  expect_equal(epred, expected)
})
