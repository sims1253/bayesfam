# Custom brms family Logit-Normal in median parametrization.

Custom brms family Logit-Normal in median parametrization.

## Usage

``` r
logitnormal(link = "identity", link_sigma = "log")
```

## Arguments

- link:

  Link function argument (as string) for Median argument. Left as
  identity!

- link_sigma:

  Link function argument (as string) for Shape argument

## Value

Logitnormal brms model-object

## Examples

``` r
a <- rnorm(1000)
data <- list(a = a, y = rlogitnormal(1000, 0.5 * a + 1, 2))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = logitnormal(), stanvars = logitnormal()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Start sampling
plot(fit)
```
