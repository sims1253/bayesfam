# Lookup for limits of family auxiliary parameters.

Lookup for limits of family auxiliary parameters.

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
