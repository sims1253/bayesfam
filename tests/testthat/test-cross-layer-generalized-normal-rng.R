# Issue #29: the injected Stan generalized_normal_rng used gamma_lccdf (a log
# survival probability in (0, 1]) as if it were a quantile, bounding all draws
# within mu +/- sigma. The fix draws the magnitude as
# sigma * gamma_rng(1/beta, 1)^(1/beta) with a uniform sign, matching the R
# side. For the standard generalized normal, |z|^beta ~ Gamma(1/beta, rate 1)
# with z = (x - mu) / sigma, so
#   P(|X - mu| > sigma) = pgamma(1, shape = 1/beta, lower.tail = FALSE)
#   Var(X) = sigma^2 * Gamma(3/beta) / Gamma(1/beta)
# The old construction puts every draw inside mu +/- sigma and must fail these
# checks for every beta.

check_generalized_normal_rng_tails <- function(draws, mu, sigma, beta, label) {
  n <- length(draws)
  p_tail <- pgamma(1, shape = 1 / beta, lower.tail = FALSE)
  frac <- mean(abs(draws - mu) > sigma)
  se <- sqrt(p_tail * (1 - p_tail) / n)
  expect_lt(
    abs(frac - p_tail),
    5 * se,
    label = sprintf(
      "%s (beta = %s): tail mass %.4f, expected %.4f",
      label,
      beta,
      frac,
      p_tail
    )
  )
  expect_gt(frac, 0, label = sprintf("%s (beta = %s): draws must exceed mu +/- sigma", label, beta))

  var_true <- sigma^2 * gamma(3 / beta) / gamma(1 / beta)
  rel_tol <- if (beta < 0.75) 0.2 else 0.12
  expect_lt(
    abs(var(draws) / var_true - 1),
    rel_tol,
    label = sprintf(
      "%s (beta = %s): variance %.2f, expected %.2f",
      label,
      beta,
      var(draws),
      var_true
    )
  )

  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  q_true <- qgeneralized_normal(p, mu = mu, sigma = sigma, beta = beta)
  q_emp <- quantile(draws, p)
  expect_true(
    all(abs(q_emp - q_true) < pmax(0.03 * sigma, 0.15 * abs(q_true))),
    label = sprintf("%s (beta = %s): empirical quantiles off", label, beta)
  )
}

test_that("R-side generalized_normal RNG has the density-implied tails", {
  n <- 40000
  mu <- 1
  sigma <- 2
  for (beta in c(0.5, 1, 2, 4)) {
    set.seed(9102 + 10 * beta)
    draws <- rgeneralized_normal(n, mu = mu, sigma = sigma, beta = beta)
    check_generalized_normal_rng_tails(draws, mu, sigma, beta, "rgeneralized_normal")
  }
})

test_that("compiled Stan generalized_normal_rng has the density-implied tails", {
  skip_if_no_stan_toolchain()

  entry <- cross_layer_inventory[["generalized_normal"]]
  program <- paste0(
    "functions {\n",
    cross_layer_stan_scode(entry),
    "\n}\ndata {\n",
    "  real mu;\n  real sigma;\n  real beta;\n  int<lower=1> N;\n",
    "}\nmodel { }\n",
    "generated quantities {\n",
    "  vector[N] draws;\n",
    "  for (m in 1:N) draws[m] = generalized_normal_rng(mu, sigma, beta);\n",
    "}\n"
  )
  model <- compile_cross_layer_stan(program, "generalized-normal-rng")

  n <- 40000
  mu <- 1
  sigma <- 2
  for (beta in c(0.5, 1, 2, 4)) {
    fit <- run_fixed_param(
      model,
      data = list(mu = mu, sigma = sigma, beta = beta, N = n),
      seed = 7300 + 10 * beta
    )
    draws <- stan_gq(fit, "draws")
    check_generalized_normal_rng_tails(
      draws,
      mu,
      sigma,
      beta,
      "Stan generalized_normal_rng"
    )
  }
})
