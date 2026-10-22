# Lookup for limits of family auxiliary parameters.

Lookup for limits of family auxiliary parameters.

## Usage

``` r
aux_limits_lookup(family)
```

## Arguments

- family:

  The identifier string of a family. Unknown identifiers raise an error.

## Value

List containing lower and upper bounds for the auxiliary parameters, in
the order of
[`aux_family_parameters_lookup()`](https://sims1253.github.io/bayesfam/pr/51/reference/aux_family_parameters_lookup.md).

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
