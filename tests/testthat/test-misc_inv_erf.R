test_that("inverse gaussian error function round-trip", {
  # erf is odd, so the round trip has to recover x on both sides of 0
  x <- seq(from = -3, to = 3, length.out = 100)
  expect_eps(inv_erf(erf(x)), x, eps = 1e-10)
})

test_that("inv_erf known values", {
  expect_equal(inv_erf(0), 0)
  expect_eps(inv_erf(0.5), 0.4769362762, eps = 1e-9)
  expect_eps(inv_erf(-0.5), -0.4769362762, eps = 1e-9)
  # erf(1) ~= 0.8427007929497149
  expect_eps(inv_erf(erf(1)), 1, eps = 1e-12)
  expect_eps(inv_erf(erf(-1)), -1, eps = 1e-12)
  # symmetry
  expect_eps(inv_erf(-0.3), -inv_erf(0.3), eps = 1e-15)
  # boundaries
  expect_equal(inv_erf(1), Inf)
  expect_equal(inv_erf(-1), -Inf)
})
