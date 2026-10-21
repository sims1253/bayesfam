# Mean parameterization of the gamma pdf.

Mean parameterization of the gamma pdf.

## Usage

``` r
dgamma_mean(x, mu, a, log = FALSE)
```

## Arguments

- x:

  Value space, x \> 0.

- mu:

  Mean parameter of the density, mu \> 0.

- a:

  Shape parameter, a \> 0.

- log:

  optional argument. If true, returns logarithmic probability. Default =
  FALSE

## Value

f(x \| mu, k)

## Details

Define rate constant rho as: \$\$\rho(\alpha, \mu) =
\frac{\alpha}{\mu}\$\$

The Frechet distribution density is defined as \$\$f(y) =
\frac{y^{\alpha - 1} exp(-\frac{y}{\rho})} {\rho^\alpha \Gamma(\alpha)}
\$\$

## Examples

``` r
x <- seq(from = 0.01, to = 10, length.out = 1000)
plot(x, dgamma_mean(x, mu = 2, a = 2), type = "l")
```
