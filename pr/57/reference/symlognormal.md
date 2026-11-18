# Custom brms family symlog-Normal

Custom brms family symlog-Normal

## Usage

``` r
symlognormal(link = "identity", link_sigma = "log")
```

## Source

Based on Hafner, D., Pasukonis, J., Ba, J., & Lillicrap, T. (2023).
Mastering Diverse Domains through World Models.
(<https://doi.org/10.48550/arXiv.2301.04104>)

## Arguments

- link:

  Link function argument (as string) for Median argument. Left as
  identity!

- link_sigma:

  Link function argument (as string) for Shape argument

## Value

symlognormal brms model-object

## Examples

``` r
a <- rnorm(1000)
data <- list(a = a, y = rsymlognormal(1000, 0.5 * a + 1, 2))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = symlognormal(), stanvars = symlognormal()$stanvars,
  refresh = 0
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(fit)
#> Error: object 'fit' not found
```
