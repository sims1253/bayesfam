# Uses euler metric for denominator

Finite operands are normalized by their largest absolute magnitude
before subtracting and squaring, which avoids squaring large raw values
to Inf (see issue \#35). Identical values, including equal infinities,
have distance 0. Pairs, where both magnitudes are below
`sqrt(.Machine$double.xmin)`, also count as distance 0 (the euler metric
is meaningless at that scale). Otherwise, an infinite value against a
finite one has distance 1 and opposite-sign infinities have distance
sqrt(2).

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
