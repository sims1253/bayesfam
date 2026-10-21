# logistic brms custom family

logistic brms custom family

## Usage

``` r
logistic(link = "identity", link_sigma = "log")
```

## Arguments

- link:

  Link function for function

- link_sigma:

  Link function for sigma argument

## Value

BRMS logistic distribution family

## Examples

``` r
a <- rnorm(n = 1000)
data <- list(a = a, y = rlogistic(n = 1000, mu = a + 2, sigma = 2))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = logistic(), stanvars = logistic()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Start sampling
plot(fit)
```
