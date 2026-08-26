# The Centered Skew Normal Does Not Write Its Expected Information Out

Returns `FALSE`, by asking
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
the same question. The registration of
[`distrib_expected_hessian.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal2Distrib.md)
says where the arithmetic is assembled, not that it is closed form: it
is a chain onto the parent, whose own expected information is the base
class's quadrature.

## Arguments

- x:

  A `SkewNormal2Distrib` object.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

`FALSE`, a logical of length 1.

## Details

Reading the owning class would answer `TRUE` here and be wrong, which is
this predicate's whole reason for existing as a generic. The cost
separates the two cases by four orders of magnitude: measured at 100
observations, this family costs 5220 milliseconds and the parent it
chains onto 2230, where the families that do write their information out
answer in a median of 0.183 milliseconds.

Reported as exact, it made
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reject a legitimate `fisher_scoring(approx = )` here with a message
saying the family computes its expected information in closed form,
which is untrue.

## See also

[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
for the generic and the rule it encodes, and
[`distrib_expected_hessian.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal2Distrib.md)
for the method it describes.

## Examples

``` r
# Neither parametrization of the skew normal writes its information out.
eh <- distributions7:::expected_hessian_exact
c(centered = eh(skewnormal2_distrib()),
  direct = eh(skewnormal1_distrib()),
  gaussian = eh(gaussian1_distrib()))
#> centered   direct gaussian 
#>    FALSE    FALSE     TRUE 
```
