# Gaussion Error function

Gaussion Error function

## Usage

``` r
erf(x)
```

## Arguments

- x:

  value to be transformed, x unbound

## Value

erf function of x, result e (0, 1)

## Examples

``` r
x <- seq(from = -2, to = 2, length.out = 100)
plot(x, erf(x), type = "l")
```
