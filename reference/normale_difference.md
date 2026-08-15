# Uses euler metric for denominator

Uses euler metric for denominator

## Usage

``` r
normale_difference(va, vb)
```

## Arguments

- va:

  Numeric scalar or vector of entries

- vb:

  Numeric scalar or vector of entries If both va and vb are no scalars,
  their lengths have to be equal

## Value

Vector of normalized differences

## Examples

``` r
print(bayesfam:::normale_difference(c(1, 1, 1, 1, 1), c(-1, 0, 1, 2, 3)))
#> [1] 1.4142136 1.0000000 0.0000000 0.4472136 0.6324555
```
