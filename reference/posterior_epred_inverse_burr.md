# Expectation-Predict vignette for the Inverse Burr distribution, with Median parametrization.

Expectation-Predict vignette for the Inverse Burr distribution, with
Median parametrization.

## Usage

``` r
posterior_epred_inverse_burr(prep)
```

## Arguments

- prep:

  brms data

## Value

Mean of the Inverse Burr distribution given data in prep; Inf for draws
with gamma \<= 1, where the mean does not exist
