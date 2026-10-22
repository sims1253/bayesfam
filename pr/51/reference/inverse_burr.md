# Inverse Burr Stan-implementation in Median parametrization.

Inverse Burr Stan-implementation in Median parametrization.

## Usage

``` r
inverse_burr(link = "log", link_tau = "log", link_gamma = "log")
```

## Arguments

- link:

  Link function for function

- link_tau:

  Link function for tau argument

- link_gamma:

  Link function for gamma argument

## Value

brms Inverse Burr distribution family

## Examples

``` r
a <- rnorm(1000)
data <- list(a = a, y = rinverse_burr(1000, exp(0.5 * a + 1), 2, 3))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = inverse_burr(), stanvars = inverse_burr()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Start sampling
plot(fit)
```
