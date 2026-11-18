# Softplus link function

Softplus link function

## Usage

``` r
softplus(x)
```

## Arguments

- x:

  value to be transformed, x positive unbound

## Value

softplus function of x, result is unbound

## Examples

``` r
x <- seq(from = 0, to = 5, length.out = 100)
plot(x, softplus(x), type = "l")
```
