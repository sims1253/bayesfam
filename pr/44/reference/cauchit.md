# Cauchit link function

Cauchit link function

## Usage

``` r
cauchit(x)
```

## Arguments

- x:

  value of x to be transformed, not defined for x of whole integer

## Value

cauchit value of x, result is unbound

## Examples

``` r
x <- seq(from = 0.1, to = 0.9, length.out = 100)
plot(x, cauchit(x), type = "l")
```
