# Lookup function for the names of the auxiliary parameters of a likelihood

Returns the auxiliary distributional parameter names of a family, i.e.
all dpars except mu. Zero-auxiliary families return an empty character
vector. Unknown identifiers raise an error.

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
