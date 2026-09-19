# Lognormal Natural brms family

Lognormal Natural brms family

## Usage

``` r
lognormal_natural(link = "log", link_sigma = "log")
```

## Arguments

- link:

  link for mu, default = log

- link_sigma:

  link for sigma, default = log

## Value

lognormal natural brms family

## Examples

``` r
a <- rnorm(n = 1000)
data <- list(a = a, y = rlognormal_natural(n = 1000, mu = exp(0.5 * a + 1), sigma = exp(2)))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = lognormal_natural(), stanvars = lognormal_natural()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(fit)
#> Error: object 'fit' not found
```
