# Check that \|a - b\| \< eps. Works with scalars and vectors on any input. For vector input, one may defined, how many \|a - b\| \>= eps are acceptable with r argument. Given it is a "workhorse" used in almost all test-files, this function is also very well tested.

Check that \|a - b\| \< eps. Works with scalars and vectors on any
input. For vector input, one may defined, how many \|a - b\| \>= eps are
acceptable with r argument. Given it is a "workhorse" used in almost all
test-files, this function is also very well tested.

## Usage

``` r
expect_eps(a, b, eps, r = 0, relative = FALSE, note = NULL, debug = FALSE)
```

## Arguments

- a:

  numeric scalar or vector a to be compared

- b:

  numeric scalar or vector b to be compared

- eps:

  numeric scalar or vector, setting the max differences, eps \> 0

- r:

  optional numeric scalar (r = 0 in default), relative number of values,
  that may have an difference \> eps.

- relative:

  bool argument, if set will take the normale difference with euler
  metric. Default = FALSE

- note:

  optional parameter used for debugging.

- debug:

  bool argument, if set will printout the difference, in case of
  failure. Default = FALSE Used internally, too calculate acceptable
  amount of deviances. Calculated absolute value will be floored.

## Value

success or failure with message

## Details

For vector/scalar combinations, allowed are (r has to always be a
scalar!):

- all scalar

- a vector, b scalar, eps scalar -\> compare a values to be close enough
  to b

- a scalar, b vector, eps scalar -\> compare b values to be close enough
  to a

- a scalar, b scalar, eps vector -\> compare the difference a, b too all
  eps -\> case not intended, but would work

- a vector, b scalar, eps vector -\> compare each vector-difference
  entry against each eps

- all vectors -\> each entry of \|a - b\| is compared to the same entry
  in eps -\> different vector lengths != 1
  dissallowed!expect_brms_family

## Examples

``` r
print(bayesfam:::expect_eps(1, 1.1, 0.2)) # should pass
#> <expectation_success/expectation/condition>
#> As expected
# print(expect_error(bayesfam:::expect_eps(c(0, 1, 3), c(1, 1, 2), 1e-4, 1 / 3)))
# should fail (2/3 were wrong, but only 1/3 was allowed)
```
