test_that("symlognormal", {
   # Check of brms can fit the custom family and recover the intercept and shape
  expect_brms_family(
    intercept = symlog(0),
    aux_par = 1,
    ref_intercept = 0,
    parameter_link = symlog,
    rng_link = identity,
    family = symlognormal,
    rng = rsymlognormal,
    aux_name = "sigma"
  )
})
