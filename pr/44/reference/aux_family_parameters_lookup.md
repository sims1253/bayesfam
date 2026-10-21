# Lookup function for the names of the auxiliary parameters of a likelihood

Lookup function for the names of the auxiliary parameters of a
likelihood

## Usage

``` r
aux_family_parameters_lookup(family)
```

## Arguments

- family:

  The identifier string of a family.

## Value

Character vector of auxiliary parameter names.

## Examples

``` r
aux_family_parameters_lookup("beta")
#> [1] "phi"
```
