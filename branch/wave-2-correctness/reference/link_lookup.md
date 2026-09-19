# Lookup function for link and repsonse functions via string identifier.

If a transformed normal likelihood is passed, the respective built-in
link will be returned instead of the identity link that would commonly
be used with transformed normal likelihood families.

## Usage

``` r
link_lookup(link, family = NULL, inv = FALSE)
```

## Arguments

- link:

  String identifier for the link function of interest.

- family:

  If a transformed normal family is passed, returns the respective link
  instead of `link`

- inv:

  True to return the response function instead of the link function.

## Value

The respective link function.

## Examples

``` r

link_lookup("log", "gaussian", FALSE)
#> function (x, base = exp(1))  .Primitive("log")

link_lookup("identiy", "logitnormal", FALSE)
#> function (x) 
#> {
#>     if (any(x < 0 | x > 1)) {
#>         stop("The logit link is only defined between 0 and 1!")
#>     }
#>     return(qlogis(x))
#> }
#> <bytecode: 0x55f31728b998>
#> <environment: namespace:bayesfam>
```
