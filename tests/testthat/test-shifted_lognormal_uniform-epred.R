test_that("posterior_epred_shifted_lognormal_uniform expands observation data by row (#25)", {
  # S != N with unequal per-observation bounds and constant parameters
  # across draws: every row must be identical (the old recycling bug swapped
  # the middle row's bounds between observations).
  S <- 3
  N <- 2
  prep <- structure(list(
    ndraws = S,
    dpars = list(
      mu = matrix(0.5, S, N),
      sigma = matrix(0.7, S, N),
      mix = matrix(0.1, S, N),
      shiftprop = matrix(0.5, S, N)
    ),
    data = list(vreal1 = c(1, 10), vreal2 = c(5, 20))
  ), class = "brmsprep")

  epred <- posterior_epred_shifted_lognormal_uniform(prep)
  expect_true(is.matrix(epred))
  expect_equal(dim(epred), c(S, N))

  expected_row <- function(j) {
    shift <- 0.5 * c(1, 10)[j]
    0.1 * (0.5 * c(5, 20)[j]) +
      (1 - 0.1) * (shift + exp(0.5 + 0.7^2 / 2))
  }
  for (s in seq_len(S)) {
    expect_eps(
      epred[s, ],
      vapply(seq_len(N), expected_row, numeric(1)),
      eps = 1e-12
    )
  }
})

test_that("posterior_epred_shifted_lognormal_uniform handles varying parameters (#25)", {
  # predicted (per-draw and per-observation) parameters as S x N matrices
  S <- 4
  N <- 3
  mu <- matrix(seq_len(S * N) * 0.1, S, N)
  sigma <- matrix(0.5 + seq_len(S * N) * 0.01, S, N)
  mix <- matrix(0.2, S, N)
  shiftprop <- matrix(0.25, S, N)
  max_shift <- c(1, 2, 4)
  max_uniform <- c(4, 8, 16)
  prep <- structure(list(
    ndraws = S,
    dpars = list(
      mu = mu,
      sigma = sigma,
      mix = mix,
      shiftprop = shiftprop
    ),
    data = list(vreal1 = max_shift, vreal2 = max_uniform)
  ), class = "brmsprep")

  max_shift_m <- matrix(max_shift, nrow = S, ncol = N, byrow = TRUE)
  max_uniform_m <- matrix(max_uniform, nrow = S, ncol = N, byrow = TRUE)
  expected <- mix * (0.5 * max_uniform_m) +
    (1 - mix) * (shiftprop * max_shift_m + exp(mu + sigma^2 / 2))

  expect_equal(
    posterior_epred_shifted_lognormal_uniform(prep),
    expected
  )

  # draw-level parameter vectors (length S) apply row-wise
  prep_draw <- structure(list(
    ndraws = S,
    dpars = list(
      mu = seq_len(S) * 0.1,
      sigma = rep(0.7, S),
      mix = rep(0.1, S),
      shiftprop = rep(0.5, S)
    ),
    data = list(vreal1 = max_shift, vreal2 = max_uniform)
  ), class = "brmsprep")
  out_draw <- posterior_epred_shifted_lognormal_uniform(prep_draw)
  expect_equal(dim(out_draw), c(S, N))
  for (s in seq_len(S)) {
    expected_s <- 0.1 * (0.5 * max_uniform) +
      0.9 * (0.5 * max_shift + exp(s * 0.1 + 0.7^2 / 2))
    expect_eps(out_draw[s, ], expected_s, eps = 1e-12)
  }
})

test_that("posterior_epred_shifted_lognormal_uniform works for a single observation (#25)", {
  S <- 5
  prep <- structure(list(
    ndraws = S,
    dpars = list(
      mu = matrix(seq_len(S) * 0.2, S, 1),
      sigma = matrix(1, S, 1),
      mix = matrix(0.3, S, 1),
      shiftprop = matrix(0.4, S, 1)
    ),
    data = list(vreal1 = 2, vreal2 = 10)
  ), class = "brmsprep")

  epred <- posterior_epred_shifted_lognormal_uniform(prep)
  expect_equal(dim(epred), c(S, 1))
  expected <- 0.3 * 5 + 0.7 * (0.8 + exp(seq_len(S) * 0.2 + 0.5))
  expect_eps(as.vector(epred), expected, eps = 1e-12)

  # scalar parameters also recycle to S x N
  prep_scalar <- structure(list(
    ndraws = S,
    dpars = list(mu = 0, sigma = 1, mix = 0.2, shiftprop = 0.5),
    data = list(vreal1 = 1, vreal2 = 6)
  ), class = "brmsprep")
  epred_scalar <- posterior_epred_shifted_lognormal_uniform(prep_scalar)
  expect_equal(dim(epred_scalar), c(S, 1))
  expect_true(all(epred_scalar == 0.2 * 3 + 0.8 * (0.5 + exp(0.5))))
})
