# Method-of-Moments Estimates, Family by Family

The parameters a family's own first two moments imply for a sample, in
closed form where the inversion has one.

## Usage

``` r
moment_estimates(distrib, y)
```

## Arguments

- distrib:

  A univariate distribution.

- y:

  The response.

## Value

A named list of parameters, or `NULL` where this family has no entry.

## Details

A starting value should be an estimate, not a guess, and for most
families the moment estimate is one line: the mean and the variance of
the sample are set equal to the family's own and the pair is inverted.
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
then refines it by maximum likelihood, and statmodels7 takes the result
as the intercept of each equation.

The inversions are written against each family's
[`mean()`](https://rdrr.io/r/base/mean.html) and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
rather than from memory, and the tests check them that way: a family's
moment estimate applied to a large sample from a known parameter must
return that parameter.

Where the family has no moments the robust analogue is used instead –
the Cauchy takes the median and half the interquartile range, which is
what its location and scale are. Where the inversion needs a root the
standard approximation is used, the Weibull's shape from the coefficient
of variation and the von Mises's concentration from the mean resultant
length; a starting value is allowed to be approximate. Families whose
inversion is neither closed nor standard – the generalized gamma, the
skew normal in its direct parametrization, the skew t, the elastic net –
are not listed and fall back to
[`start_from_moments`](https://statmodels7.github.io/distributions7/reference/start_from_moments.md)'s
reading of `params_interpretation`.

## See also

[`distrib_start`](https://statmodels7.github.io/distributions7/reference/distrib_start.md),
[`start_from_moments`](https://statmodels7.github.io/distributions7/reference/start_from_moments.md)
