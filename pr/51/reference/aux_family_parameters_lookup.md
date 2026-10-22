# Lookup function for the names of the auxiliary parameters of a likelihood

Lookup function for the names of the auxiliary parameters of a
likelihood

## Usage

``` r
aux_family_parameters_lookup(family)
```

## Arguments

- family:

  The identifier string of a family. Unknown identifiers raise an error.

## Value

Character vector of auxiliary parameter names (all distributional
parameters except `mu`), empty if there are none.

## Examples

``` r
aux_family_parameters_lookup("beta")
#> [1] "phi"
```
