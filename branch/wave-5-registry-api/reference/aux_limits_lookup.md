# Lookup for limits of family auxiliary parameters.

Returns a list with numeric lower and upper bound vectors of the
auxiliary parameters, matching the order of
aux_family_parameters_lookup(). Zero-auxiliary families return empty
numeric vectors. Unknown identifiers raise an error.

## Usage

``` r
aux_limits_lookup(family)
```

## Arguments

- family:

  The identifier string of a family.

## Value

List containing lower and upper bounds for the auxiliary parameter.

## Examples

``` r
aux_limits_lookup("beta")
#> $lb
#> [1] 0
#> 
#> $ub
#> [1] Inf
#> 
```
