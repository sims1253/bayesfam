# Issue #30: dbeta_mean used log(gamma(phi)) which overflows to Inf for
# phi >~ 171, producing Inf/NaN densities for ordinary precision values.

test_that("dbeta_mean matches stats::dbeta for large precision on a grid (#30)", {
  x <- c(1e-4, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 1 - 1e-4)
  mu_values <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  phi_values <- c(0.5, 1, 2, 10, 100, 171, 200, 500, 1000, 5000)
  for (mu in mu_values) {
    for (phi in phi_values) {
      ref_log <- stats::dbeta(
        x,
        shape1 = mu * phi,
        shape2 = (1 - mu) * phi,
        log = TRUE
      )
      expect_equal(dbeta_mean(x, mu = mu, phi = phi, log = TRUE), ref_log)
      ref <- stats::dbeta(x, shape1 = mu * phi, shape2 = (1 - mu) * phi)
      expect_equal(dbeta_mean(x, mu = mu, phi = phi), ref)
    }
  }
})

test_that("dbeta_mean is finite and correct at the issue's examples (#30)", {
  lpdf <- dbeta_mean(0.5, mu = 0.5, phi = 200, log = TRUE)
  expect_true(is.finite(lpdf))
  expect_equal(lpdf, stats::dbeta(0.5, shape1 = 100, shape2 = 100, log = TRUE))
  expect_equal(lpdf, 2.422117, tolerance = 1e-6)

  lpdf1000 <- dbeta_mean(0.5, mu = 0.5, phi = 1000, log = TRUE)
  expect_true(is.finite(lpdf1000))
  expect_false(is.nan(lpdf1000))
  expect_equal(
    lpdf1000,
    stats::dbeta(0.5, shape1 = 500, shape2 = 500, log = TRUE)
  )
})

test_that("dbeta_mean ordinary scale is finite for large precision (#30)", {
  dens <- dbeta_mean(c(0.25, 0.5, 0.75), mu = 0.5, phi = 1000)
  expect_true(all(is.finite(dens)))
  expect_true(all(dens > 0))
})
