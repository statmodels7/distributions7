# Method-of-Moments Estimates, Family by Family

Returns the parameters a family's own first two moments imply for a
sample, in closed form where the inversion has one, and `NULL` where
this family has no entry. **37 of the 42 univariate families have one.**

A starting value should be an estimate. For most families the moment
estimate is one line: the sample mean and variance are set equal to the
family's own and the pair is inverted.
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
then refines it by maximum likelihood, and statmodels7 takes the result
as the intercept of each equation.

## Usage

``` r
moment_estimates(distrib, y)
```

## Arguments

- distrib:

  A univariate distribution, read for its `distrib_name` and its
  parameter names.

- y:

  The response, a numeric vector, already filtered to finite values by
  the caller.

## Value

A named list of parameters on the parameter scale, one component per
entry of `distrib@params`, or `NULL` where this family has no entry.

## How the inversions are written and checked

Each is written against the family's own
[`mean()`](https://rdrr.io/r/base/mean.html) and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
not from memory, and the tests check it the same way: a family's moment
estimate applied to a large sample drawn from a known parameter must
return that parameter.

## Where a moment does not exist, and where the inversion is not closed

A family with no moments takes the robust analogue. The Cauchy uses the
median and half the interquartile range, which are its location and its
scale. Where the inversion needs a root, the standard approximation
stands in: the Weibull's shape from the coefficient of variation, the
von Mises's concentration from the mean resultant length. A starting
value is allowed to be approximate.

Five families have no entry and fall through to
[`start_from_moments()`](https://statmodels7.github.io/distributions7/reference/start_from_moments.md)'s
reading of `params_interpretation`: `gengamma1`, `gengamma2`, `skewt`,
`pseudohuber` and `enet`. Their inversions are neither closed nor
standard: two gamma-function ratios in two shapes for the generalized
gamma, a moment that does not identify the pair for the skew \\t\\, no
elementary moments at all for the elastic net.

## A fixed constant in the name

A family carrying one announces it there, as in
`"beta-binomial [size=10]"`, so the bracketed part is stripped before
the name is used as a lookup key. Without that no beta-binomial of any
size would match its own entry, and silently: a missing key is a legal
fallback.

## See also

[`start_from_moments()`](https://statmodels7.github.io/distributions7/reference/start_from_moments.md),
the only caller and the fallback;
[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
for the generic; [`mean()`](https://rdrr.io/r/base/mean.html) and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
which the inversions are written against.
