# Whether Two Information Objects Select the Same Matrix

Returns `TRUE` when both objects are
[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md),
or both are
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
with the same `method` (two `NULL` values included).

## Usage

``` r
same_information(a, b)
```

## Arguments

- a, b:

  [Information](https://statmodels7.github.io/distributions7/reference/Information.md)
  objects.

## Value

A single logical.

## See also

[`default_information()`](https://statmodels7.github.io/distributions7/reference/default_information.md)
