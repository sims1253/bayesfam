test_that("pkumaraswamy implements the Kumaraswamy CDF (#26)", {
  mu_grid <- seq(from = 0.05, to = 0.95, length.out = 7)
  p_grid <- c(0.5, 1, 2, 3.7, 10)

  # mu parameterizes the median, so F(mu) must be 0.5 for every (mu, p)
  for (mu in mu_grid) {
    for (p in p_grid) {
      expect_eps(pkumaraswamy(mu, mu = mu, p = p), 0.5, eps = 1e-12)
    }
  }

  # CDF is monotone and within [0, 1] and agrees with the integrated density
  for (mu in c(0.2, 0.5, 0.85)) {
    for (p in p_grid) {
      x <- seq(from = 1e-6, to = 1 - 1e-6, length.out = 501)
      cdf <- pkumaraswamy(x, mu = mu, p = p)
      expect_true(all(diff(cdf) >= 0))
      expect_true(all(cdf >= 0 & cdf <= 1))
      # integrated density at selected points
      for (xq in c(0.1, 0.3, mu, 0.7, 0.9)) {
        expect_eps(
          pkumaraswamy(xq, mu = mu, p = p),
          integrate(
            dkumaraswamy,
            lower = 1e-9,
            upper = xq,
            mu = mu,
            p = p,
            subdivisions = 1000L,
            rel.tol = 1e-10
          )$value,
          eps = 1e-4
        )
      }
    }
  }

  # the previously broken case: F(0.5) at mu = 0.5 for non-integer q
  expect_eps(pkumaraswamy(0.5, mu = 0.5, p = 2), 0.5, eps = 1e-12)

  # CDF / quantile round trip
  for (mu in c(0.2, 0.5, 0.85)) {
    for (p in p_grid) {
      u <- seq(from = 0.01, to = 0.99, length.out = 25)
      qq <- qkumaraswamy(u, mu = mu, p = p)
      expect_eps(pkumaraswamy(qq, mu = mu, p = p), u, eps = 1e-8)
    }
  }

  # vectorized input works
  x <- c(0.1, 0.4, 0.8)
  expect_equal(
    pkumaraswamy(x, mu = 0.4, p = 2),
    vapply(x, function(xi) pkumaraswamy(xi, mu = 0.4, p = 2), numeric(1))
  )

  # error checks stay intact
  expect_error(pkumaraswamy(1.1, mu = 0.5, p = 1))
  expect_error(pkumaraswamy(0.5, mu = 0, p = 1))
  expect_error(pkumaraswamy(0.5, mu = 0.5, p = 0))
})
