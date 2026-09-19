# Probability density function for the Lomax distribution, with Mean parametrization.

Probability density function for the Lomax distribution, with Mean
parametrization.

## Usage

``` r
dlomax(x, mu, alpha, log = FALSE)
```

## Arguments

- x:

  Model space, defined for x \>= 0

- mu:

  Mean parameter of pdf, mu \> 0

- alpha:

  Alpha parameter of pdf, alpha \> 1

- log:

  Optional log argument, if true, return log(pdf)

## Value

PDF of Lomax Distribution

## Details

\$\$f(y \| \mu, \alpha) = \frac{\alpha}{\mu(\alpha-1)} (1 +
\frac{y}{\mu(\alpha-1)})^{-\alpha-1}\$\$

## Examples

``` r
x <- seq(from = 0, to = 10, length.out = 100)
plot(x, dlomax(x, mu = 1, alpha = 2), type = "l")
```
