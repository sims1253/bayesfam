# Unit Lindley brms family

Unit Lindley brms family

## Usage

``` r
unit_lindley(link = "logit")
```

## Arguments

- link:

  link for mu, default = logit

## Value

Unit Lindley brms family

## Examples

``` r
a <- rnorm(n = 1000)
data <- list(a = a, y = runit_lindley(n = 1000, mu = inv_logit(0.5 * a + 1)))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = unit_lindley(), stanvars = unit_lindley()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(fit)
#> Error: object 'fit' not found
```
