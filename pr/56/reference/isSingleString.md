# Check, if the input is a single string

Check, if the input is a single string

## Usage

``` r
isSingleString(input)
```

## Arguments

- input:

  String argument

## Value

Is a string and only one string

## Examples

``` r
bayesfam:::isSingleString("abc") # should be TRUE
#> [1] TRUE
bayesfam:::isSingleString(c("abc", "def")) # should be FALSE, not a single string
#> [1] FALSE
```
