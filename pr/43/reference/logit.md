# Logit link function

Logit link function

## Usage

``` r
logit(x)
```

## Arguments

- x:

  value of x to be transformed, x e (0, 1)

## Value

logit value of x, x unbound

## Examples

``` r
x <- seq(from = 0.1, to = 0.9, length.out = 100)
plot(x, logit(x), type = "l")
```
