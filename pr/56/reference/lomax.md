# Lomax Stan-implementation in Mean parametrization.

Lomax Stan-implementation in Mean parametrization.

## Usage

``` r
lomax(link = "log", link_alpha = "logm1")
```

## Arguments

- link:

  Link function for function

- link_alpha:

  Link function for the alpha argument. The default `"logm1"` (inverse
  `exp(eta) + 1`) keeps alpha above 1.

## Value

brms Lomax distribution family

## Details

With a `"logm1"` link the Stan code includes brms' overloaded `expp1`
and `logm1` functions, which
[`rstan::expose_stan_functions()`](https://mc-stan.org/rstan/reference/expose_stan_functions.html)
cannot export, so expose `lomax_lpdf` and `lomax_rng` from a family with
other links.

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
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(fit)
#> Error: object 'fit' not found
```
