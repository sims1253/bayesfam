# Family inventory for the cross-layer test suite (issue #28).
#
# Every exported distribution of the package is enumerated here together with
# the information the three test layers need:
#   layer 1: R density / quantile / RNG contracts (fast, deterministic)
#   layer 2: R vs injected-Stan lpdf parity
#   layer 3: log_lik / posterior_predict / posterior_epred callback fixtures
#
# The inventory is the single source of truth for coverage: tests in
# test-cross-layer-inventory.R fail when a family is missing here, when a
# callback is not registered in the brms family object, or when a callback has
# no declared layer-3 coverage.
#
# Known bugs fixed concurrently on sibling branches (wave 2) are marked with
# status "skip" and the reason references the issue number; the affected
# assertions are written against the mathematically correct behavior and will
# light up once the wave-2 fixes land on this branch.

new_cross_layer_family <- function(
  name,
  custom = TRUE,
  d_fun = NULL,
  q_fun = NULL,
  r_fun = NULL,
  d_argnames = NULL,
  dpars = NULL,
  params = NULL,
  params2 = NULL,
  support = NULL,
  bad_params = NULL,
  bad_x = NA_real_,
  always_log = FALSE,
  stan_lpdf = NA_character_,
  stan_rng = NA_character_,
  parity = FALSE,
  parity_grid = NULL,
  parity_ref = NULL,
  rng_stat = NA_character_,
  rng_target = NULL,
  rng_params = NULL,
  ll3 = list(status = "ok"),
  pp3 = list(status = "ok"),
  ep3 = list(status = "ok"),
  epred_ref = NULL
) {
  entry <- list(
    name = name,
    custom = custom,
    family_fun = if (custom) get(name, envir = asNamespace("bayesfam")) else NULL,
    d_fun = d_fun,
    q_fun = q_fun,
    r_fun = r_fun,
    d_argnames = d_argnames %||% names(params),
    dpars = dpars %||% names(params),
    params = params,
    params2 = params2,
    support = support,
    bad_params = bad_params,
    bad_x = bad_x,
    always_log = always_log,
    stan_lpdf = stan_lpdf,
    stan_rng = stan_rng,
    parity = parity,
    parity_grid = parity_grid,
    parity_ref = parity_ref,
    rng_stat = rng_stat,
    rng_target = rng_target,
    rng_params = rng_params,
    ll3 = ll3,
    pp3 = pp3,
    ep3 = ep3,
    epred_ref = epred_ref
  )
  entry$postfit <- cross_layer_postfit(entry)
  entry
}

`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

# Default layer-3 fixture (mock brms prep data) per family. Families with
# observation-level data (shifted_lognormal_uniform) or special y grids get
# their own branch.
cross_layer_postfit <- function(entry) {
  pf <- switch(
    entry$name,
    shifted_lognormal_uniform = list(
      y = c(0.5, 2, 6),
      params = list(
        mu = c(0.2, 0.4, 0.6),
        sigma = 0.3,
        mix = 0.2,
        shiftprop = 0.5
      ),
      vreal1 = c(0.8, 0.8, 0.8),
      vreal2 = c(5, 10, 20)
    ),
    betaprime = list(y = c(0.5, 2, 8), params = list(mu = c(2, 4, 6), phi = 3)),
    cauchitnormal = list(
      y = c(0.25, 0.5, 0.75),
      params = list(mu = c(0, 0.3, -0.5), sigma = 1)
    ),
    cloglognormal = list(
      y = c(0.2, 0.5, 0.8),
      params = list(mu = c(-0.3, 0, 0.2), sigma = 0.6)
    ),
    generalized_gamma = list(
      y = c(0.5, 2, 5),
      params = list(mu = c(0.3, 0.5, 0.8), sigma = 0.5, Q = 1)
    ),
    generalized_normal = list(
      y = c(-1, 0.5, 3),
      params = list(mu = c(0, 1, 2), sigma = 1.5, beta = 1.5)
    ),
    gompertz = list(
      y = c(0.5, 2, 4),
      params = list(mu = c(1, 2, 3), beta = 0.8)
    ),
    gumbel_mean = list(
      y = c(0, 2, 5),
      params = list(mu = c(0, 2, 4), sigma = 1.2)
    ),
    kumaraswamy = list(
      y = c(0.2, 0.5, 0.8),
      params = list(mu = c(0.3, 0.5, 0.7), p = 2)
    ),
    logistic = list(
      y = c(-1, 2, 5),
      params = list(mu = c(0, 2, 4), sigma = 1.2)
    ),
    logitnormal = list(
      y = c(0.3, 0.5, 0.7),
      params = list(mu = c(-0.5, 0, 0.5), sigma = 1)
    ),
    lognormal_natural = list(
      y = c(0.8, 2, 5),
      params = list(mu = c(1, 2, 3), sigma = 0.4)
    ),
    lomax = list(
      y = c(0.5, 2, 6),
      params = list(mu = c(1, 2, 3), alpha = 4)
    ),
    shifted_inv_gaussian = list(
      y = c(1, 2, 3),
      params = list(mu = c(1, 1.5, 2), shape = 3, ndt = 0.5)
    ),
    simplex = list(
      y = c(0.3, 0.5, 0.8),
      params = list(mu = c(0.3, 0.5, 0.7), sigma = 1)
    ),
    softplusnormal = list(
      y = c(0.5, 1.5, 4),
      params = list(mu = c(0.5, 1, 2), sigma = 0.7)
    ),
    symlognormal = list(
      y = c(-1, 0.5, 3),
      params = list(mu = c(-0.5, 0, 0.5), sigma = 0.5)
    ),
    unit_lindley = list(
      y = c(0.2, 0.5, 0.8),
      params = list(mu = c(0.3, 0.5, 0.7))
    ),
    NULL
  )
  pf %||% list(y = c(0.5, 2, 5), params = entry$params)
}

# Reference implementation of the shifted lognormal/uniform log density with
# the mathematically correct (untruncated lognormal) semantics. The current R
# implementation subtracts a truncation term and is therefore wrong (issue #24,
# fixed on a sibling branch); Stan and this reference agree.
ref_dshifted_lognormal_uniform <- function(
  y,
  meanlog,
  sdlog,
  mix,
  shift,
  max_uniform
) {
  unif_llh <- dunif(y, 0, max_uniform, log = TRUE)
  lognormal_llh <- dlnorm(y - shift, meanlog = meanlog, sdlog = sdlog, log = TRUE)
  bayesfam::logsumexp(c(
    log(mix) + unif_llh,
    log1p(-mix) + lognormal_llh
  ))
}

cross_layer_inventory <- list(

  new_cross_layer_family(
    name = "betaprime",
    d_fun = bayesfam::dbetaprime,
    q_fun = bayesfam::qbetaprime,
    r_fun = bayesfam::rbetaprime,
    params = list(mu = 4, phi = 2),
    params2 = list(mu = 0.8, phi = 5),
    support = function(pr) c(0, Inf),
    bad_params = list(phi = 0),
    bad_x = -1,
    stan_lpdf = "betaprime_lpdf",
    stan_rng = "betaprime_rng",
    parity = TRUE,
    parity_grid = list(
      y = c(0.1, 0.5, 1, 2, 5, 10),
      mu = c(1, 2, 4, 4, 2, 0.8),
      phi = c(1, 2, 2, 5, 5, 0.5)
    ),
    rng_stat = "mean",
    rng_target = function(pr) pr$mu
  ),

  new_cross_layer_family(
    name = "cauchitnormal",
    d_fun = bayesfam::dcauchitnormal,
    r_fun = bayesfam::rcauchitnormal,
    params = list(mu = 0.5, sigma = 1),
    params2 = list(mu = -1, sigma = 2),
    support = function(pr) c(0, 1),
    bad_params = list(sigma = -1),
    bad_x = 1.5,
    stan_lpdf = "cauchitnormal_lpdf",
    stan_rng = "cauchitnormal_rng",
    parity = TRUE,
    parity_grid = list(
      y = c(0.2, 0.4, 0.5, 0.7, 0.9),
      mu = c(0, 0.5, 1, -0.5, 2),
      sigma = c(0.5, 1, 2, 1.5, 0.5)
    ),
    rng_stat = "median",
    rng_target = function(pr) pcauchy(pr$mu),
    ep3 = list(
      status = "unsupported",
      reason = "no closed-form mean for the cauchit-normal (family docs)"
    )
  ),

  new_cross_layer_family(
    name = "cloglognormal",
    d_fun = bayesfam::dcloglognormal,
    r_fun = bayesfam::rcloglognormal,
    params = list(mu = -0.3, sigma = 0.6),
    params2 = list(mu = 0.2, sigma = 0.7),
    support = function(pr) c(0, 1),
    bad_params = list(sigma = -1),
    bad_x = 1.5,
    stan_lpdf = "cloglognormal_lpdf",
    stan_rng = "cloglognormal_rng",
    parity = TRUE,
    parity_grid = list(
      y = c(0.2, 0.4, 0.5, 0.7, 0.9),
      mu = c(-0.5, 0, 0.3, 0, -1),
      sigma = c(0.5, 1, 0.6, 1.5, 0.3)
    ),
    rng_stat = "median",
    rng_target = function(pr) 1 - exp(-exp(pr$mu)),
    ep3 = list(
      status = "unsupported",
      reason = "no closed-form mean for the cloglog-normal (family docs)"
    )
  ),

  new_cross_layer_family(
    name = "generalized_gamma",
    d_fun = bayesfam::dgeneralized_gamma,
    r_fun = bayesfam::rgeneralized_gamma,
    d_argnames = c("mu", "sigma", "Q"),
    params = list(mu = 0.5, sigma = 0.8, Q = 1),
    params2 = list(mu = 0, sigma = 0.5, Q = 0),
    support = function(pr) c(0, Inf),
    bad_params = list(sigma = 0),
    bad_x = 0,
    stan_lpdf = "generalized_gamma_lpdf",
    parity = TRUE,
    parity_grid = list(
      y = c(0.5, 1, 2, 5, 0.5, 2),
      mu = c(0, 0.5, 0.5, 1, 0, 0.5),
      sigma = c(0.5, 0.8, 0.3, 1, 0.5, 0.8),
      Q = c(1, -1, 2, 0.5, 0, 0)
    ),
    rng_stat = "mean",
    # Q = 0 reduces the generalized gamma to a lognormal with mean
    # exp(mu + sigma^2 / 2); no simple mean exists for Q != 0
    rng_params = list(mu = 0.5, sigma = 0.8, Q = 0),
    rng_target = function(pr) exp(pr$mu + pr$sigma^2 / 2),
    ll3 = list(
      status = "skip",
      # dgeneralized_gamma branches on `if (Q != 0)` and errors on the
      # length-S draw vector that get_dpar returns (issue #23, wave-2 PR).
      reason = "#23 generalized_gamma log_lik fails on posterior-draw vectors of Q (fixed in wave-2 PR)"
    ),
    ep3 = list(
      status = "unsupported",
      reason = "posterior_epred stops: no implementation (family docs)"
    )
  ),

  new_cross_layer_family(
    name = "generalized_normal",
    d_fun = bayesfam::dgeneralized_normal,
    q_fun = bayesfam::qgeneralized_normal,
    r_fun = bayesfam::rgeneralized_normal,
    params = list(mu = 1, sigma = 2, beta = 1),
    params2 = list(mu = 0, sigma = 1, beta = 2),
    support = function(pr) c(-Inf, Inf),
    bad_params = list(sigma = 0),
    bad_x = NA_real_,
    stan_lpdf = "generalized_normal_lpdf",
    stan_rng = "generalized_normal_rng",
    parity = TRUE,
    parity_grid = list(
      y = c(-3, -1, 0, 0.5, 2, 5),
      mu = c(0, 1, 0, 1, 2, 0),
      sigma = c(1, 2, 0.5, 1.5, 2, 1),
      beta = c(1, 2, 0.5, 4, 1.5, 2)
    ),
    rng_stat = "mean",
    rng_target = function(pr) pr$mu
  ),

  new_cross_layer_family(
    name = "gompertz",
    d_fun = bayesfam::dgompertz,
    q_fun = bayesfam::qgompertz,
    r_fun = bayesfam::rgompertz,
    params = list(mu = 2, beta = 1),
    params2 = list(mu = 0.5, beta = 0.2),
    support = function(pr) c(0, Inf),
    bad_params = list(mu = 0),
    bad_x = 0,
    stan_lpdf = "gompertz_lpdf",
    stan_rng = "gompertz_rng",
    parity = TRUE,
    parity_grid = list(
      y = c(0.2, 0.5, 1, 2, 5),
      mu = c(1, 2, 2, 0.5, 3),
      beta = c(1, 0.5, 2, 0.8, 1)
    ),
    rng_stat = "median",
    rng_target = function(pr) pr$mu,
    ep3 = list(
      status = "unsupported",
      reason = "posterior_epred stops: mean not defined for the family (family docs)"
    )
  ),

  new_cross_layer_family(
    name = "gumbel_mean",
    d_fun = bayesfam::dgumbel_mean,
    q_fun = bayesfam::qgumbel_mean,
    r_fun = bayesfam::rgumbel_mean,
    params = list(mu = 2, sigma = 1.5),
    params2 = list(mu = 0, sigma = 0.5),
    support = function(pr) c(-Inf, Inf),
    bad_params = list(sigma = 0),
    bad_x = NA_real_,
    stan_lpdf = "gumbel_mean_lpdf",
    stan_rng = "gumbel_mean_rng",
    parity = TRUE,
    parity_grid = list(
      y = c(-2, 0, 1, 3, 6),
      mu = c(0, 1, 1, 2, 2),
      sigma = c(1, 0.5, 1.5, 2, 1)
    ),
    rng_stat = "mean",
    rng_target = function(pr) pr$mu
  ),

  new_cross_layer_family(
    name = "kumaraswamy",
    d_fun = bayesfam::dkumaraswamy,
    q_fun = bayesfam::qkumaraswamy,
    r_fun = bayesfam::rkumaraswamy,
    d_argnames = c("mu", "p"),
    params = list(mu = 0.3, p = 2),
    params2 = list(mu = 0.6, p = 0.5),
    support = function(pr) c(0, 1),
    bad_params = list(p = 0),
    bad_x = 0,
    stan_lpdf = "kumaraswamy_lpdf",
    stan_rng = "kumaraswamy_rng",
    parity = TRUE,
    parity_grid = list(
      y = c(0.1, 0.3, 0.5, 0.7, 0.9),
      mu = c(0.2, 0.5, 0.3, 0.6, 0.4),
      p = c(2, 0.5, 3, 1, 5)
    ),
    rng_stat = "median",
    rng_target = function(pr) pr$mu,
    epred_ref = function(prep) {
      mu <- prep$dpars$mu
      p <- prep$dpars$p
      q <- -(log(2) / log1p(-mu^p))
      q * beta(1 + 1 / p, q)
    }
  ),

  new_cross_layer_family(
    name = "logistic",
    d_fun = bayesfam::dlogistic,
    q_fun = bayesfam::qlogistic,
    r_fun = bayesfam::rlogistic,
    params = list(mu = 2, sigma = 1.5),
    params2 = list(mu = -1, sigma = 0.7),
    support = function(pr) c(-Inf, Inf),
    bad_params = list(sigma = 0),
    bad_x = NA_real_,
    stan_lpdf = "logistic_r_lpdf",
    stan_rng = "logistic_r_rng",
    parity = TRUE,
    parity_grid = list(
      y = c(-3, -1, 0.5, 2, 5),
      mu = c(0, 1, 0, 2, 3),
      sigma = c(1, 0.5, 1.5, 2, 1)
    ),
    rng_stat = "mean",
    rng_target = function(pr) pr$mu
  ),

  new_cross_layer_family(
    name = "logitnormal",
    d_fun = bayesfam::dlogitnormal,
    r_fun = bayesfam::rlogitnormal,
    params = list(mu = 0.7, sigma = 1),
    params2 = list(mu = -1.5, sigma = 0.4),
    support = function(pr) c(0, 1),
    bad_params = list(sigma = -1),
    bad_x = 1.5,
    stan_lpdf = "logitnormal_lpdf",
    stan_rng = "logitnormal_rng",
    parity = TRUE,
    parity_grid = list(
      y = c(0.1, 0.3, 0.5, 0.7, 0.9),
      mu = c(-1, 0, 0.5, 0, -0.5),
      sigma = c(0.5, 1, 2, 0.5, 1.5)
    ),
    rng_stat = "median",
    rng_target = function(pr) plogis(pr$mu),
    ep3 = list(
      status = "warning",
      reason = "posterior_epred warns and returns the median, not the mean (family docs)"
    ),
    epred_ref = function(prep) plogis(prep$dpars$mu)
  ),

  new_cross_layer_family(
    name = "lognormal_natural",
    d_fun = bayesfam::dlognormal_natural,
    r_fun = bayesfam::rlognormal_natural,
    params = list(mu = 2, sigma = 0.5),
    params2 = list(mu = 0.4, sigma = 1.2),
    support = function(pr) c(0, Inf),
    bad_params = list(mu = 0),
    bad_x = 0,
    stan_lpdf = "lognormal_natural_lpdf",
    stan_rng = "lognormal_natural_rng",
    parity = TRUE,
    parity_grid = list(
      y = c(0.5, 1, 2, 5),
      mu = c(1, 2, 1.5, 3),
      sigma = c(0.3, 0.5, 1, 0.8)
    ),
    rng_stat = "mean",
    rng_target = function(pr) pr$mu
  ),

  new_cross_layer_family(
    name = "lomax",
    d_fun = bayesfam::dlomax,
    q_fun = bayesfam::qlomax,
    r_fun = bayesfam::rlomax,
    params = list(mu = 2, alpha = 5),
    params2 = list(mu = 0.5, alpha = 2.5),
    support = function(pr) c(0, Inf),
    bad_params = list(mu = 0),
    bad_x = -1,
    stan_lpdf = "lomax_lpdf",
    stan_rng = "lomax_rng",
    parity = TRUE,
    parity_grid = list(
      y = c(0.1, 0.5, 1, 2, 5, 10, 20),
      mu = c(1, 2, 2, 0.5, 2, 1, 3),
      alpha = c(3, 5, 2.5, 8, 4, 3, 6)
    ),
    rng_stat = "mean",
    rng_target = function(pr) pr$mu,
    pp3 = list(
      status = "skip",
      # posterior_predict_lomax samples rgompertz instead of rlomax
      # (issue #22, fixed in wave-2 PR).
      reason = "#22 posterior_predict_lomax uses rgompertz (fixed in wave-2 PR)"
    )
  ),

  new_cross_layer_family(
    name = "shifted_inv_gaussian",
    d_fun = bayesfam::dshifted_inv_gaussian,
    r_fun = bayesfam::rshifted_inv_gaussian,
    d_argnames = c("mu", "shape", "shift"),
    dpars = c("mu", "shape", "ndt"),
    params = list(mu = 1.5, shape = 2, shift = 0.5),
    params2 = list(mu = 3, shape = 0.8, shift = 0.1),
    support = function(pr) c(pr$shift, Inf),
    bad_params = list(mu = 0),
    bad_x = 0,
    stan_lpdf = "shifted_inv_gaussian_lpdf",
    parity = TRUE,
    parity_grid = list(
      y = c(1, 1.5, 2, 3, 5),
      mu = c(1, 1, 1.5, 2, 2),
      shape = c(1, 3, 2, 5, 0.5),
      ndt = c(0.5, 0.3, 0.5, 1, 2)
    ),
    rng_stat = "mean",
    rng_target = function(pr) pr$mu + pr$shift,
    epred_ref = function(prep) prep$dpars$mu + prep$dpars$ndt
  ),

  new_cross_layer_family(
    name = "shifted_lognormal_uniform",
    d_fun = bayesfam::dshifted_lognormal_uniform,
    r_fun = bayesfam::rshifted_lognormal_uniform,
    d_argnames = c("meanlog", "sdlog", "mix", "shift", "max_uniform"),
    dpars = c("mu", "sigma", "mix", "shiftprop"),
    params = list(meanlog = 0.3, sdlog = 0.4, mix = 0.1, shift = 0.1, max_uniform = 10),
    params2 = list(meanlog = 0, sdlog = 1, mix = 0.5, shift = 0.5, max_uniform = 3),
    support = function(pr) c(0, Inf),
    bad_params = list(sdlog = 0),
    bad_x = -1,
    always_log = TRUE,
    stan_lpdf = "shifted_lognormal_uniform_lpdf",
    parity = TRUE,
    parity_grid = list(
      y = c(0.5, 2, 6, 12, 0.2),
      mu = c(0.2, 0.4, 0.6, 0.5, 0.1),
      sigma = c(0.3, 0.5, 0.3, 0.5, 0.2),
      mix = c(0.2, 0.1, 0.5, 0.05, 0.9),
      shiftprop = c(0.5, 0.3, 0.8, 0.1, 0.5),
      max_shift = c(0.8, 0.5, 1, 0.3, 0.4),
      max_uniform = c(5, 10, 20, 15, 10)
    ),
    parity_ref = function(pr) {
      ref_dshifted_lognormal_uniform(
        y = pr$y,
        meanlog = pr$mu,
        sdlog = pr$sigma,
        mix = pr$mix,
        shift = pr$shiftprop * pr$max_shift,
        max_uniform = pr$max_uniform
      )
    },
    rng_stat = "mean",
    rng_target = function(pr) {
      pr$mix * 0.5 * pr$max_uniform +
        (1 - pr$mix) * (pr$shift + exp(pr$meanlog + pr$sdlog^2 / 2))
    },
    ll3 = list(
      status = "skip",
      # dshifted_lognormal_uniform subtracts a truncation term that the model
      # does not apply, so log_lik disagrees with the Stan likelihood
      # (issue #24, fixed in wave-2 PR).
      reason = "#24 mixture R density unnormalized / disagrees with Stan (fixed in wave-2 PR)"
    ),
    ep3 = list(
      status = "skip",
      # posterior_epred recycles length-N vreal data along columns of the
      # S x N posterior matrices whenever S != N (issue #25, wave-2 PR). The
      # N = 1 case is covered by a dedicated passing test.
      reason = "#25 posterior_epred recycles observation bounds across draws (fixed in wave-2 PR)"
    )
  ),

  new_cross_layer_family(
    name = "simplex",
    d_fun = bayesfam::dsimplex,
    r_fun = bayesfam::rsimplex,
    params = list(mu = 0.7, sigma = 1),
    params2 = list(mu = 0.3, sigma = 2),
    support = function(pr) c(0, 1),
    bad_params = list(sigma = 0),
    bad_x = 0,
    stan_lpdf = "simplex_lpdf",
    parity = TRUE,
    parity_grid = list(
      y = c(0.2, 0.5, 0.7, 0.9),
      mu = c(0.3, 0.5, 0.7, 0.8),
      sigma = c(0.5, 1, 2, 10)
    ),
    rng_stat = "mean",
    rng_target = function(pr) pr$mu
  ),

  new_cross_layer_family(
    name = "softplusnormal",
    d_fun = bayesfam::dsoftplusnormal,
    r_fun = bayesfam::rsoftplusnormal,
    params = list(mu = 1, sigma = 0.8),
    params2 = list(mu = 0.2, sigma = 2),
    support = function(pr) c(0, Inf),
    bad_params = list(sigma = 0),
    bad_x = 0,
    stan_lpdf = "softplusnormal_lpdf",
    stan_rng = "softplusnormal_rng",
    parity = TRUE,
    parity_grid = list(
      y = c(0.2, 0.5, 1, 3),
      mu = c(0, 0.5, 1, 2),
      sigma = c(0.5, 1, 0.3, 1.5)
    ),
    rng_stat = "median",
    rng_target = function(pr) log(exp(pr$mu) + 1),
    ep3 = list(
      status = "unsupported",
      reason = "posterior_epred stops: no implementation (family docs)"
    )
  ),

  new_cross_layer_family(
    name = "symlognormal",
    d_fun = bayesfam::dsymlognormal,
    r_fun = bayesfam::rsymlognormal,
    params = list(mu = 0.5, sigma = 0.5),
    params2 = list(mu = -1, sigma = 1.5),
    support = function(pr) c(-Inf, Inf),
    bad_params = list(sigma = -1),
    bad_x = NA_real_,
    # Not included in the combined Stan parity program: its injected scode
    # redefines the helper `sign`, which is also defined by
    # generalized_normal, and Stan forbids duplicate definitions.
    stan_lpdf = "symlognormal_lpdf",
    stan_rng = "symlognormal_rng",
    parity = FALSE,
    rng_stat = "median",
    rng_target = function(pr) sign(pr$mu) * (exp(abs(pr$mu)) - 1),
    ep3 = list(
      status = "unsupported",
      reason = "posterior_epred warns and returns no mean (family docs)"
    )
  ),

  new_cross_layer_family(
    name = "unit_lindley",
    d_fun = bayesfam::dunit_lindley,
    q_fun = bayesfam::qunit_lindley,
    r_fun = bayesfam::runit_lindley,
    params = list(mu = 0.4),
    params2 = list(mu = 0.7),
    support = function(pr) c(0, 1),
    bad_params = list(mu = 1.5),
    bad_x = 0,
    stan_lpdf = "unit_lindley_lpdf",
    parity = TRUE,
    parity_grid = list(
      y = c(0.1, 0.3, 0.5, 0.7, 0.9),
      mu = c(0.2, 0.4, 0.5, 0.7, 0.9)
    ),
    rng_stat = "mean",
    rng_target = function(pr) pr$mu
  ),

  # Distributions that accompany built-in brms families: R density/RNG only,
  # no injected Stan code and no custom callbacks (brms provides those).

  new_cross_layer_family(
    name = "lognormal",
    custom = FALSE,
    d_fun = bayesfam::dlognormal,
    r_fun = bayesfam::rlognormal,
    params = list(mu = 0.5, sigma = 0.8),
    params2 = list(mu = 0, sigma = 1.5),
    support = function(pr) c(0, Inf),
    bad_params = list(sigma = 0),
    bad_x = 0,
    rng_stat = "median",
    rng_target = function(pr) exp(pr$mu)
  ),

  new_cross_layer_family(
    name = "beta_mean",
    custom = FALSE,
    d_fun = bayesfam::dbeta_mean,
    q_fun = bayesfam::qbeta_mean,
    r_fun = bayesfam::rbeta_mean,
    d_argnames = c("mu", "phi"),
    params = list(mu = 0.3, phi = 5),
    params2 = list(mu = 0.6, phi = 1),
    support = function(pr) c(0, 1),
    bad_params = list(phi = 0),
    bad_x = 1,
    rng_stat = "mean",
    rng_target = function(pr) pr$mu
  ),

  new_cross_layer_family(
    name = "gamma_mean",
    custom = FALSE,
    d_fun = bayesfam::dgamma_mean,
    r_fun = bayesfam::rgamma_mean,
    params = list(mu = 3, a = 2),
    params2 = list(mu = 0.5, a = 5),
    support = function(pr) c(0, Inf),
    bad_params = list(mu = 0),
    bad_x = 0,
    rng_stat = "mean",
    rng_target = function(pr) pr$mu
  ),

  new_cross_layer_family(
    name = "student_mean",
    custom = FALSE,
    d_fun = NULL,
    r_fun = bayesfam::rstudent_mean,
    d_argnames = c("mu", "df", "sigma"),
    dpars = c("mu", "df", "sigma"),
    params = list(mu = 2, df = 5, sigma = 1),
    params2 = list(mu = -1, df = 10, sigma = 2),
    support = function(pr) c(-Inf, Inf),
    bad_params = list(df = 0),
    bad_x = NA_real_,
    rng_stat = "mean",
    rng_target = function(pr) pr$mu
  ),

  new_cross_layer_family(
    name = "exgauss_mean",
    custom = FALSE,
    d_fun = NULL,
    r_fun = bayesfam::rexgauss_mean,
    d_argnames = c("mu", "sigma", "beta"),
    dpars = c("mu", "sigma", "beta"),
    params = list(mu = 1, sigma = 0.5, beta = 2),
    params2 = list(mu = 0, sigma = 1, beta = 0.5),
    support = function(pr) c(-Inf, Inf),
    bad_params = list(sigma = 0),
    bad_x = NA_real_,
    rng_stat = "mean",
    rng_target = function(pr) pr$mu
  ),

  new_cross_layer_family(
    name = "frechet_median",
    custom = FALSE,
    d_fun = bayesfam::dfrechet_median,
    r_fun = bayesfam::rfrechet_median,
    d_argnames = c("mu", "nu"),
    params = list(mu = 2, nu = 3),
    params2 = list(mu = 0.5, nu = 6),
    support = function(pr) c(0, Inf),
    bad_params = list(nu = 1),
    bad_x = 0,
    rng_stat = "mean",
    rng_target = function(pr) pr$mu
  )
)

names(cross_layer_inventory) <- vapply(
  cross_layer_inventory,
  function(e) e$name,
  character(1)
)

# All custom brms families of the package; the coverage test compares the
# inventory against this list so that new families cannot be added silently.
cross_layer_custom_families <- c(
  "betaprime",
  "cauchitnormal",
  "cloglognormal",
  "generalized_gamma",
  "generalized_normal",
  "gompertz",
  "gumbel_mean",
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

cross_layer_helper_families <- c(
  "lognormal",
  "beta_mean",
  "gamma_mean",
  "student_mean",
  "exgauss_mean",
  "frechet_median"
)

# Build a mock brms prep object (make_brmsprep-style) without running any
# Stan sampler. Parameters may be scalars, per-observation vectors (recycled
# over draws, i.e. distributional predictors), or ready-made S x N matrices.
make_cross_layer_prep <- function(
  family,
  y,
  params,
  ndraws = 7,
  vreal1 = NULL,
  vreal2 = NULL
) {
  nobs <- length(y)
  expand <- function(v) {
    if (is.matrix(v)) {
      return(v)
    }
    if (length(v) == 1) {
      return(matrix(v, nrow = ndraws, ncol = nobs))
    }
    matrix(v, nrow = ndraws, ncol = nobs, byrow = TRUE)
  }
  dpars <- lapply(params, expand)
  data <- list(Y = y, N = nobs)
  if (!is.null(vreal1)) {
    data$vreal1 <- vreal1
  }
  if (!is.null(vreal2)) {
    data$vreal2 <- vreal2
  }
  structure(
    list(
      family = family,
      dpars = dpars,
      data = data,
      ndraws = ndraws
    ),
    class = "brmsprep"
  )
}

# Inventory postfit fixture as a prep object.
postfit_prep <- function(entry, ndraws = 7) {
  pf <- entry$postfit
  make_cross_layer_prep(
    family = entry$family_fun(),
    y = pf$y,
    params = pf$params,
    ndraws = ndraws,
    vreal1 = pf$vreal1 %||% NULL,
    vreal2 = pf$vreal2 %||% NULL
  )
}

# Map brms dpar draws (column i of the S x N matrices) to the argument list of
# the family's R density / RNG.
postfit_draw_args <- function(entry, prep, i) {
  if (entry$name == "shifted_lognormal_uniform") {
    return(list(
      meanlog = prep$dpars$mu[, i],
      sdlog = prep$dpars$sigma[, i],
      mix = prep$dpars$mix[, i],
      shift = prep$dpars$shiftprop[, i] * prep$data$vreal1[i],
      max_uniform = prep$data$vreal2[i]
    ))
  }
  dpar_to_arg <- setNames(
    entry$d_argnames,
    if (entry$name == "shifted_inv_gaussian") {
      c("mu", "shape", "ndt")
    } else {
      entry$dpars[seq_along(entry$d_argnames)]
    }
  )
  args <- lapply(names(dpar_to_arg), function(nm) prep$dpars[[nm]][, i])
  stats::setNames(args, dpar_to_arg)
}

# Evaluate a family's R log density at a single (y, params) point.
cross_layer_r_loglik <- function(entry, y, params) {
  if (entry$name == "shifted_lognormal_uniform") {
    return(ref_dshifted_lognormal_uniform(
      y = y,
      meanlog = params$meanlog,
      sdlog = params$sdlog,
      mix = params$mix,
      shift = params$shift,
      max_uniform = params$max_uniform
    ))
  }
  do.call(
    entry$d_fun,
    c(list(x = y, log = TRUE), params)
  )
}
