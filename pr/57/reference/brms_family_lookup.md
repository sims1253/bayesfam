# Lookup function for brms families via string identifier

Lookup function for brms families via string identifier

## Usage

``` r
brms_family_lookup(family, link = NULL)
```

## Arguments

- family:

  String identifier of the family. Unknown identifiers raise an error.

- link:

  Link to be passed to the family function. If `NULL`, the family's
  default link is used.

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
