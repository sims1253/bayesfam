
eps <- 1e-6
n_testset <- 100

# Test the link functions
test_that("logit-link", {
  # check, that the error is within reason
  expect_eps(0, logit(0.5), eps)
  expect_eps(-1.38629436112, logit(0.2), eps)
  expect_eps(0.847297860387, logit(0.7), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(0, -1.38629436112, 0.847297860387), logit(c(0.5, 0.2, 0.7)), eps)
  # check link against link from another package
  testset_logit <- seq(from = eps, to = 1 - eps, length.out = n_testset)
  expect_eps(brms:::logit(testset_logit), logit(testset_logit), eps)
  # check values at the boundary of the defined space
  expect_equal(logit(0), -Inf)
  expect_equal(logit(1), Inf)
  # check, that non-numeric arguments result in an error
  expect_error(logit("R"))
  # check, that wrong number of arguments produce an error
  expect_error(logit(0.1, 0.5))
  # check ranges
  expect_error(logit(-1))
  expect_error(logit(2))
})

test_that("inverse-logit-link", {
  # check, that the error is within reason
  expect_eps(0.5, inv_logit(0), eps)
  expect_eps(0.119202922022, inv_logit(-2), eps)
  expect_eps(0.73105857863, inv_logit(1), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(0.5, 0.119202922022, 0.73105857863), inv_logit(c(0, -2, 1)), eps)
  # check link against link from another package
  testset_invlogit <- seq(from = -10, to = 10, length.out = n_testset)
  expect_eps(brms:::inv_logit(testset_invlogit), inv_logit(testset_invlogit), eps)
  # check values, for x to Inf
  expect_equal(1, inv_logit(Inf))
  expect_equal(0, inv_logit(-Inf))
  # check, that non-numeric arguments result in an error
  expect_error(inv_logit("R"))
  # check, that wrong number of arguments produce an error
  expect_error(logit(0.1, 0.5))
})

test_that("cloglog-link", {
  # check, that the error is within reason
  expect_eps(-9.21029, cloglog(0.0001), eps)
  expect_eps(-2.250367, cloglog(0.1), eps)
  expect_eps(-0.005764308, cloglog(0.63), eps)
  expect_eps(0.475885, cloglog(0.8), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(-9.21029, -2.250367, -0.005764308, 0.475885), cloglog(c(0.0001, 0.1, 0.63, 0.8)), eps)
  # check link against link from another package
  testset_cloglog <- seq(from = eps, to = 1 - eps, length.out = n_testset)
  expect_eps(brms:::cloglog(testset_cloglog), cloglog(testset_cloglog), eps)
  # check boundary values
  expect_equal(-Inf, cloglog(0))
  expect_equal(Inf, cloglog(1))
  # check, that non-numeric arguments result in an error
  expect_error(cloglog("R"))
  # check, that wrong number of arguments produce an error
  expect_error(cloglog(0.1, 0.5))
  # check ranges
  expect_error(cloglog(-1))
  expect_error(cloglog(2))
})

test_that("inverse-cloglog-link", {
  # check, that the error is within reason
  expect_eps(0.04856801, inv_cloglog(-3), eps)
  expect_eps(0.3077994, inv_cloglog(-1), eps)
  expect_eps(0.6321206, inv_cloglog(0), eps)
  expect_eps(0.934012, inv_cloglog(1), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(0.04856801, 0.3077994, 0.6321206, 0.934012), inv_cloglog(c(-3, -1, 0, 1)), eps)
  # check against another package
  testset_invcloglog <- seq(from = -10, to = 5, length.out = n_testset)
  expect_eps(brms:::inv_cloglog(testset_invcloglog), inv_cloglog(testset_invcloglog), eps)
  # check, values of x to Inf, on the boundary of the defined space
  expect_equal(0, inv_cloglog(-Inf))
  expect_equal(1, inv_cloglog(Inf))
  # check, that non-numeric arguments result in an error
  expect_error(inv_cloglog("R"))
  # check, that wrong number of arguments produce an error
  expect_error(inv_cloglog(0.1, 0.5))
})

test_that("cauchit-link", {
  # check, that the error is within reason
  expect_eps(-3.077684, cauchit(0.1), eps)
  expect_eps(-0.3249197, cauchit(0.4), eps)
  expect_eps(0, cauchit(0.5), eps)
  expect_eps(1.376382, cauchit(0.8), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(-3.077684, -0.3249197, 0, 1.376382), cauchit(c(0.1, 0.4, 0.5, 0.8)), eps)
  # check against another package
  testset_cauchit <- seq(from = eps, to = 1 - eps, length.out = n_testset)
  expect_eps(qcauchy(testset_cauchit), cauchit(testset_cauchit), 10 * eps)
  # check values on the boundary, should logically be Inf, but just approach Inf
  expect_equal(cauchit(0), -Inf)
  expect_equal(cauchit(1), Inf)
  # check ranges
  expect_error(cauchit(-1))
  expect_error(cauchit(2))
  # check, that non-numeric arguments result in an error
  expect_error(cauchit("R"))
  # check, that wrong number of arguments produce an error
  expect_error(cauchit(0.1, 0.5))
})

test_that("inv-cauchit-link", {
  # check, that the error is within reason
  expect_eps(0.03172552, inv_cauchit(-10), eps)
  expect_eps(0.5, inv_cauchit(0), eps)
  expect_eps(0.9682745, inv_cauchit(10), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(0.03172552, 0.5, 0.9682745), inv_cauchit(c(-10, 0, 10)), eps)
  # equivalent to pcauchy
  testset_invcauchit <- seq(from = -10, to = 10, length.out = n_testset)
  expect_equal(pcauchy(testset_invcauchit), inv_cauchit(testset_invcauchit))
  # check x approaching Inf on boundry of defined space
  expect_equal(0, inv_cauchit(-Inf))
  expect_equal(1, inv_cauchit(Inf))
  # check, that non-numeric arguments result in an error
  expect_error(inv_cauchit("R"))
  # check, that wrong number of arguments produce an error
  expect_error(inv_cauchit(0.1, 0.5))
})

test_that("gaussian-error-function", {
  # check, that the error is within reason
  expect_eps(-0.8427008, erf(-1), eps)
  expect_equal(0, erf(0))
  expect_eps(0.2227026, erf(0.2), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(-0.8427008, 0, 0.2227026), erf(c(-1, 0, 0.2)), eps)
  # check x approaching Inf on boundry of defined space
  expect_equal(-1, erf(-Inf))
  expect_equal(1, erf(Inf))
  # check, that non-numeric arguments result in an error
  expect_error(erf("R"))
  # check, that wrong number of arguments produce an error
  expect_error(erf(0.1, 0.5))
})

test_that("Softplus link-function", {
  # check, that the error is within reason
  expect_eps(-2.252168, softplus(0.1), eps)
  expect_eps(0.5413249, softplus(1), eps)
  expect_eps(3.981515, softplus(4), eps)
  expect_equal(100, softplus(100)) # should be equal to machine precision (I think)
  # check vector as argument returns vector with same results
  expect_eps(c(-2.252168, 0.5413249, 3.981515, 100), softplus(c(0.1, 1, 4, 100)), eps)
  # check x approaching Inf on boundry of defined space
  expect_equal(-Inf, softplus(0))
  expect_equal(Inf, softplus(Inf))
  # check, that non-numeric arguments result in an error
  expect_error(softplus("R"))
  # check, that wrong number of arguments produce an error
  expect_error(softplus(0.1, 0.5))
  # check values outside the defined scope are NaN and throw warning
  expect_warning(expect_true(is.na(softplus(-1))))
})

test_that("Inverse Softplus link-function", {
  # check, that the error is within reason
  expect_eps(4.53989e-05, inv_softplus(-10), eps)
  expect_eps(0.3132617, inv_softplus(-1), eps)
  expect_eps(1.313262, inv_softplus(1), eps)
  expect_equal(100, inv_softplus(100)) # should be equal to machine precision (I think)
  # check vector as argument returns vector with same results
  expect_eps(c(4.53989e-05, 0.3132617, 1.313262, 100), inv_softplus(c(-10, -1, 1, 100)), eps)
  # check x approaching Inf on boundry of defined space
  expect_equal(0, inv_softplus(-Inf))
  expect_equal(Inf, inv_softplus(Inf))
  # check, that non-numeric arguments result in an error
  expect_error(inv_softplus("R"))
  # check, that wrong number of arguments produce an error
  expect_error(inv_softplus(0.1, 0.5))
})
