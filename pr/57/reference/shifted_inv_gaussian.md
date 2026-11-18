# Shifted Inverse Gauss brms family

Shifted Inverse Gauss brms family

## Usage

``` r
shifted_inv_gaussian(link = "log", link_shape = "log", link_ndt = "log")
```

## Arguments

- link:

  link for mu, default="log"

- link_shape:

  link for the shape, default="log"

- link_ndt:

  link for the shift, default="log"

## Value

Shifted Inverse Gauss brms family

## Examples

``` r
a <- rnorm(1000)
data <- list(
  a = a,
  y = rshifted_inv_gaussian(
    n = 1000, mu = exp(0.5 * a + 1),
    shape = 1, shift = 1
  )
)
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = shifted_inv_gaussian(), stanvars = shifted_inv_gaussian()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(fit)
#> Error: object 'fit' not found
```
