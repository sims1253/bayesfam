# Logarithm of the sum of exponentials.

A more numerically stable equivalent to `log(sum(exp(x)))`

## Usage

``` r
logsumexp(x)
```

## Source

https://en.wikipedia.org/wiki/LogSumExp#log-sum-exp_trick_for_log-domain_calculations

## Arguments

- x:

  a vector of values

## Value

log(sum(exp(x)))

## Details

An empty vector returns `-Inf` (the sum is empty). `NA`/`NaN` entries
propagate and return `NA_real_`. If any entry is `+Inf`, the result is
`Inf`; if all entries are `-Inf`, the result is `-Inf`.
