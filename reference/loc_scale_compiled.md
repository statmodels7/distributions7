# The Families Whose Location-Scale Diagonal Is Compiled

Returns the class names of the families for which
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md)
takes the diagonal of the expected information, and the derivative of
each diagonal entry in its own parameter, from compiled code. The same
code serves the scalar C registry, so the two routes agree to the last
bit.

## Usage

``` r
loc_scale_compiled()
```

## Value

A character vector of S7 class names.
