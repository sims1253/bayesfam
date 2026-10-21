# Simplex brms-implementation in median parametrization.

Simplex brms-implementation in median parametrization.

## Usage

``` r
simplex(link = "logit", link_sigma = "identity")
```

## Arguments

- link:

  Link function for function

- link_sigma:

  Link function for sigma argument

## Value

brms Beta-Custom distribution family

## Examples

``` r
a <- rnorm(1000)
data <- list(a = a, y = rsimplex(1000, brms::inv_logit_scaled(0.5 * a + 1), 2))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = simplex(), stanvars = simplex()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Start sampling
plot(fit)
```
