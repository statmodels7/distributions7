# Is a Family's Exact Expected Information Costly to Compute?

`TRUE` when the family's expected information is exact but costs far
more than its observed information, because it is obtained by an
integral or a sum for each observation or each distinct shape; `FALSE`
otherwise.

## Usage

``` r
expected_hessian_costly(x, ...)
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
which states accuracy. Two kinds of family answer `TRUE`.

The location-scale families
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md),
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
and
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
take one integral over the standardized response per distinct shape,
from
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md):
nothing where the shape is the same for every observation, one integral
per observation where it is modelled. Measured on smooth regressions at
1000 observations with the shape developed over a covariate, a fit that
takes the scoring step on this information costs 7 to 21 times one that
takes it on the observed information, at the same estimate.

The Poisson-inverse Gaussians
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
and
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
sum over the support for each observation, over a number of terms that
grows as \\\sigma\mu\\: about 25 microseconds an observation at \\\mu =
3\\, \\\sigma = 0.3\\ and 2 milliseconds at \\\mu = 30\\, \\\sigma =
3\\, against a few microseconds for the observed information. Measured
at 1000 observations, a fit on it costs 1 to 6 times one on the observed
information at the same estimate, and far more where the fit passes
through a large \\\sigma\mu\\.

statmodels7's `iwls(hessian = "auto")` therefore settles on the observed
information for these families. The default answers for a wrapper's
parent, read off its `parent_distrib` property, and `FALSE` for a family
without one.

## See also

[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md),
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md)

## Examples

``` r
expected_hessian_costly(skewt_distrib())
#> [1] TRUE
expected_hessian_costly(pig1_distrib())
#> [1] TRUE
expected_hessian_costly(gaussian1_distrib())
#> [1] FALSE
expected_hessian_costly(fixed(skewt_distrib(), nu = 6))
#> [1] TRUE
```
