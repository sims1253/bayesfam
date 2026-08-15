# Beta prime brms custom family

Beta prime brms custom family

## Usage

``` r
betaprime(link = "log", link_phi = "log")
```

## Arguments

- link:

  Link function for function

- link_phi:

  Link function for beta argument

## Value

brms beta prime distribution family

## Examples

``` r
a <- rnorm(n = 1000)
data <- list(a = a, y = rbetaprime(n = 1000, mu = exp(0.5 * a + 1), phi = 2))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = betaprime(), stanvars = betaprime()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(fit)
#> Error: object 'fit' not found
```
