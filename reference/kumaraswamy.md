# Kumaraswamy brms-implementation in median parametrization.

Kumaraswamy brms-implementation in median parametrization.

## Usage

``` r
kumaraswamy(link = "logit", link_p = "log")
```

## Arguments

- link:

  Link function for mu

- link_p:

  Link function for p argument

## Value

brms Beta-Custom distribution family

## Examples

``` r
a <- rnorm(1000)
data <- list(a = a, y = rkumaraswamy(1000, brms::inv_logit_scaled(0.5 * a + 1), 2))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = kumaraswamy(), stanvars = kumaraswamy()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Start sampling
plot(fit)
```
