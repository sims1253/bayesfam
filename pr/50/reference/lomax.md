# Lomax Stan-implementation in Mean parametrization.

Lomax Stan-implementation in Mean parametrization.

## Usage

``` r
lomax(link = "log", link_alpha = "log1p")
```

## Arguments

- link:

  Link function for function

- link_alpha:

  Link function for eta argument

## Value

brms Lomax distribution family

## Examples

``` r
a <- rnorm(1000)
data <- list(a = a, y = rlomax(1000, exp(0.5 * a + 1), 2))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = lomax(), stanvars = lomax()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Start sampling
plot(fit)
```
