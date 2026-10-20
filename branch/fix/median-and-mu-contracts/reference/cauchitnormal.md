# Custom brms family Cauchitnormal

Custom brms family Cauchitnormal

## Usage

``` r
cauchitnormal(link = "identity", link_sigma = "log")
```

## Arguments

- link:

  Link function argument (as string) for Median argument. Left as
  identity!

- link_sigma:

  Link function argument (as string) for Shape argument

## Value

Cauchitnormal brms model-object

## Examples

``` r
a <- rnorm(1000)
data <- list(a = a, y = rcauchitnormal(1000, 0.5 * a + 1, 2))
fit1 <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = cauchitnormal(), stanvars = cauchitnormal()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(fit1)
#> Error: object 'fit1' not found
```
