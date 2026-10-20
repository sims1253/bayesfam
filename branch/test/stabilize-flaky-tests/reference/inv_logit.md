# Logit response function

Logit response function

## Usage

``` r
inv_logit(x)
```

## Arguments

- x:

  value of x to be transformed, x unbound

## Value

inverse logit value of x, result is e (0, 1)

## Examples

``` r
x <- seq(from = -5, to = 5, length.out = 100)
plot(x, inv_logit(x), type = "l")
```
