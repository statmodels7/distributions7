# Variance of a Multinomial

Returns the covariance matrix of the counts by delegating to
[`mv_sigma.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MultinomialDistrib.md).
It is singular, the counts summing to the trial count, so a computation
that inverts a covariance works with the marginals or with the free
vector of the simplex.

## Arguments

- x:

  A `MultinomialDistrib` object, from
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md).
  The argument is named `x` for consistency with
  [`base::mean()`](https://rdrr.io/r/base/mean.html), whose signature
  the moment generics follow.

- theta:

  A named list of the simplex's free values on the parameter scale, each
  of length 1.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A symmetric \\p \times p\\ numeric matrix of rank \\p-1\\.

## See also

[`mean.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.MultinomialDistrib.md)
for the mean vector,
[`mv_sigma.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MultinomialDistrib.md),
which this calls, and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
for the generic.
