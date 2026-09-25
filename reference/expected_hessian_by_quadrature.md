# Is a Family's Expected Information Computed by Quadrature?

`TRUE` when the family's expected information is exact but obtained by a
quadrature over the standardized response for each distinct value of its
shape parameters, as
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md)
does; `FALSE` otherwise.

## Usage

``` r
expected_hessian_by_quadrature(x, ...)
```

## Arguments

- x:

  An object inheriting from class `"distrib"`.

- ...:

  Passed to methods.

## Value

A single logical.

## Details

The answer is a statement about cost, and it is read beside
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md),
which states accuracy. A family answering `TRUE` to both computes its
expected information exactly, at one integral per distinct shape:
nothing where the shape is the same for every observation, and one
integral per observation where it is modelled. Measured on smooth
regressions at 1000 observations with the shape developed over a
covariate, a fit that takes the scoring step on this information costs 7
to 21 times one that takes it on the observed information, at the same
estimate. statmodels7's `iwls(hessian = "auto")` therefore settles on
the observed information for such a family.

Four shipped families answer `TRUE`:
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md),
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
and
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).
The default answers for a wrapper's parent, read off its
`parent_distrib` property, and `FALSE` for a family without one.

## See also

[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md),
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md)

## Examples

``` r
expected_hessian_by_quadrature(skewt_distrib())
#> [1] TRUE
expected_hessian_by_quadrature(gaussian1_distrib())
#> [1] FALSE
expected_hessian_by_quadrature(fixed(skewt_distrib(), nu = 6))
#> [1] TRUE
```
