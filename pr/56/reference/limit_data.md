# Data limit function

Data limit function

## Usage

``` r
limit_data(data, limits)
```

## Arguments

- data:

  Data to be limited

- limits:

  Limits to be used. Vector with 2 real entries,
  `limits[1] <= limits[2]` If the lower bound does not have to be
  restricted, set it to NA and vice versa. Sets data outside those
  bounds to those bounds.

## Value

data limited by the limits

## Examples

``` r
input <- c(1, 2, 3, 4)
print(bayesfam:::limit_data(input, c(2, 3))) # lower and upper bounds
#> [1] 2 2 3 3
print(bayesfam:::limit_data(input, c(2, NA))) # only lower bound
#> [1] 2 2 3 4
```
