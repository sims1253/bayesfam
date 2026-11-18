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

An empty vector returns `-Inf` (the sum is empty). A non-finite maximum
is returned as is: `NA`/`NaN` entries propagate, otherwise any `Inf`
entry gives `Inf` and all-`-Inf` input gives `-Inf`.
