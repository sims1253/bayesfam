# Custom brms family Softplus in median parametrization.

Custom brms family Softplus in median parametrization.

## Usage

``` r
softplusnormal(link = "identity", link_sigma = "log")
```

## Arguments

- link:

  Link function argument (as string) for Median argument. Left as
  identity!

- link_sigma:

  Link function argument (as string) for Shape argument

## Value

Softplus brms model-object

## Examples

``` r
a <- rnorm(1000)
data <- list(a = a, y = rsoftplusnormal(1000, 0.5 * a + 1, 2))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = softplusnormal(), stanvars = softplusnormal()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(fit)
#> Error: object 'fit' not found
```
