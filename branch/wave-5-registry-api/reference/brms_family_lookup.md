# Lookup function for brms families via string identifier

Looks up the family in the bayesfam family registry, which covers all
families exported by bayesfam plus a set of common brms built-in
families. Unknown identifiers raise an error.

## Usage

``` r
brms_family_lookup(family, link = NULL)
```

## Arguments

- family:

  String identifier of the family.

- link:

  Link to be passed to the family function. If NULL (default), the
  constructor defaults are preserved.

## Value

A brmsfamily object matching the string identifier and using the link

## Examples

``` r
brms_family_lookup("weibull", "softplus")
#> 
#> Family: weibull 
#> Link function: softplus 
#> 
```
