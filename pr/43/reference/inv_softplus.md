# Softplus response function

Softplus response function

## Usage

``` r
inv_softplus(x)
```

## Arguments

- x:

  value to be transformed, x is unbound

## Value

inv_softplus of x, result is positive unbound

## Examples

``` r
x <- seq(from = 0.1, to = 5, length.out = 100)
plot(x, softplus(x), type = "l")
```
