# Length and data check function

Length and data check function

## Usage

``` r
lenEqual(
  list_of_vectors,
  scalars_allowed = FALSE,
  type_check = NULL,
  na_allowed = FALSE
)
```

## Arguments

- list_of_vectors:

  List of vectors of any type you want

- scalars_allowed:

  In many applications, mixing scalars with vectors may be fine, so this
  argument controls, if scalars (length == 1) would be fine to. Default
  = FALSE

- type_check:

  Function pointer argument of type to be checked. All vectors have to
  be of this type, if defined. Default = NULL

- na_allowed:

  True if test shall pass with NAs in the data.

## Value

boolean, given the arguments above

## Examples

``` r
va <- c(1, 2, 3)
vb <- c(4, 5, 6)
vc <- c(7, 8)
bayesfam:::lenEqual(list(va, vb)) # both got 3 entries
#> [1] TRUE
bayesfam:::lenEqual(list(va, vb, vc)) # not all vectors have the same number of entries
#> [1] FALSE
```
