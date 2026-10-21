# Median parameterization of the Weibull pdf.

Median parameterization of the Weibull pdf.

## Usage

``` r
dweibull_median(x, mu, k, log = FALSE)
```

## Arguments

- x:

  Value space, x \> 0.

- mu:

  Median parameter, mu \> 0.

- k:

  Shape parameter, k \> 0.

- log:

  Optional argument. If TRUE, returns log(pdf). Normally False.

## Value

f(x \| mu, k)

## Details

Define constant sigma as \$\$\sigma(\mu, k) := \mu / \Gamma(1 + 1 /
k)\$\$

The Weibull distribution density is defined as \$\$f(y) =
\frac{k}{\sigma} \* (\frac{x}{\sigma})^{\alpha - 1} \*
exp(-(\frac{x}{\sigma})^\alpha)\$\$

## Examples

``` r
x <- seq(from = 0.01, to = 10, length.out = 1000)
plot(x, dweibull_median(x, mu = 2, k = 1), type = "l")
```
