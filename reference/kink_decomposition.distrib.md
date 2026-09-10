# The Base Method: No Kink

Returns `NULL`: the family is smooth in every parameter. Every family
inherits this except the three that declare an absolute value of the
residual in their log-density.

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

## Value

`NULL`.

## See also

[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md)
for the generic.

## Examples

``` r
kink_decomposition(gaussian1_distrib())
#> NULL
```
