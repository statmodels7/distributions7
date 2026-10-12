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
which concerns accuracy. Three kinds of family return `TRUE`.

The location-scale families
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md),
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md),
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
and
[`pseudohuber2_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber2_distrib.md)
take one integral over the standardized response per distinct shape,
from
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md):
one integral where the shape is the same for every observation, and one
per observation where it is modeled. A fit that takes the scoring step
on this information then costs several times one that takes it on the
observed information, at the same estimate.

The Poisson-inverse Gaussians
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
and
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
sum over the support for each observation, over a number of terms that
grows as \\\sigma\mu\\, whereas their observed information is a closed
form. A fit on it costs more than one on the observed information at the
same estimate, and far more where the fit passes through a large
\\\sigma\mu\\.

A
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
distribution also returns `TRUE`.

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
