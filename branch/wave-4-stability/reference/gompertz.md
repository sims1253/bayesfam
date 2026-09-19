# Custom Gompertz brms-implementation in median parametrization.

Custom Gompertz brms-implementation in median parametrization.

## Usage

``` r
gompertz(link = "log", link_b = "log")
```

## Arguments

- link:

  Link function for function

- link_b:

  Link function for eta argument

## Value

brms gompertz distribution family

## Examples

``` r
a <- rnorm(1000)
data <- list(a = a, y = rgompertz(1000, mu = exp(0.5 * a + 1), beta = 0.1))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = gompertz(), stanvars = gompertz()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(fit)
#> Error: object 'fit' not found
```
