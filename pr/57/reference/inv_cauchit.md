# Cauchit response function

Cauchit response function

## Usage

``` r
inv_cauchit(x)
```

## Arguments

- x:

  value of x to be transformed, any real scalar or vector allowed

## Value

inverse cauchit of x, result is e (0, 1)

## Examples

``` r
x <- seq(from = -10, to = 10, length.out = 100)
plot(x, inv_cauchit(x), type = "l")
```
