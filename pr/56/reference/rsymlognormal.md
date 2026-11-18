# symlognormal RNG-function

symlognormal RNG-function

## Usage

``` r
rsymlognormal(n, mu = 0, sigma = 1)
```

## Source

Based on Hafner, D., Pasukonis, J., Ba, J., & Lillicrap, T. (2023).
Mastering Diverse Domains through World Models.
(<https://doi.org/10.48550/arXiv.2301.04104>)

## Arguments

- n:

  Number of draws

- mu:

  Median parameter, mu unbound, mu already symlog transformed

- sigma:

  Sigma shape parameter, sigma \> 0

## Value

n symlognormal distributed samples

## Examples

``` r
hist(rsymlognormal(100, 0.5, 2))
```
