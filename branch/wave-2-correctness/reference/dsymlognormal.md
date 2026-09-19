# symlognormal density function

symlognormal density function

## Usage

``` r
dsymlognormal(x, mu, sigma, log = FALSE)
```

## Source

Based on Hafner, D., Pasukonis, J., Ba, J., & Lillicrap, T. (2023).
Mastering Diverse Domains through World Models.
(<https://doi.org/10.48550/arXiv.2301.04104>)

## Arguments

- x:

  Value space of the distribution, x unbound

- mu:

  Median parameter, mu unbound

- sigma:

  Sigma shape parameter, sigma \>= 0

- log:

  Bool argument, if true, returns the logarithmic density

## Value

Normal distribution density with symlog link function

## Examples

``` r
x <- seq(from = -5, to = 5, length.out = 1000)
plot(x, dsymlognormal(x, mu = -0.2, sigma = 0.4), type = "l")
```
