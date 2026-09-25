# The Centered Skew Normal Answers for Its Parent

Returns what
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
returns, the expected information here being a chain onto that parent's.

## Arguments

- x:

  A `SkewNormal2Distrib` object.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A logical of length 1.

## See also

[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
for the generic and the rule it encodes, and
[`distrib_expected_hessian.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal2Distrib.md)
for the method it describes.

## Examples

``` r
eh <- distributions7:::expected_hessian_exact
c(centered = eh(skewnormal2_distrib()), direct = eh(skewnormal1_distrib()))
#> centered   direct 
#>     TRUE     TRUE 
```
