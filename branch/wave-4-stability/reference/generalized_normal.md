# Generalized Normal BRMS family

Generalized Normal BRMS family

## Usage

``` r
generalized_normal(link = "identity", link_sigma = "log", link_beta = "log")
```

## Arguments

- link:

  Link function for function

- link_sigma:

  Link function for sigma argument

- link_beta:

  Link function for beta argument

## Value

BRMS generalized_normal distribution family

## Examples

``` r
data <- list(y = rgeneralized_normal(n = 1000, mu = 2, sigma = 2, beta = 4))
fit <- brms::brm(
  formula = y ~ 1, data = data,
  family = generalized_normal(), stanvars = generalized_normal()$stanvars,
  init = 0.1
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(fit)
#> Error: object 'fit' not found
```
