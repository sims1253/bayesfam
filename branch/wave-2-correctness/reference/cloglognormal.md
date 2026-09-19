# Custom brms family Cloglog-Normal in median parametrization.

Custom brms family Cloglog-Normal in median parametrization.

## Usage

``` r
cloglognormal(link = "identity", link_sigma = "log")
```

## Arguments

- link:

  Link function argument (as string) for Median argument. Left as
  identity!

- link_sigma:

  Link function argument (as string) for Shape argument

## Value

Cloglog brms model-object

## Examples

``` r
data <- rcloglognormal(1000, 0.5, 2)
# cloglognormal does not like values to close to the boundary
data <- limit_data(data, c(1e-12, 1 - 1e-12))
fit <- brms::brm(
  formula = y ~ 1, data = list(y = data),
  family = cloglognormal(), stanvars = cloglognormal()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(fit)
#> Error: object 'fit' not found
```
