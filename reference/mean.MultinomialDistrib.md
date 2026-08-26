# Mean of a Multinomial

Returns the mean count vector \\\mathbb{E}\[Y\] = np\\ by delegating to
[`mv_location.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.MultinomialDistrib.md).
The result sums to the trial count, not to one.

## Arguments

- x:

  A `MultinomialDistrib` object, from
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md).
  The argument is named `x` because the generic is
  [`base::mean()`](https://rdrr.io/r/base/mean.html).

- theta:

  A named list of the simplex's free values on the parameter scale, each
  of length 1.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length \\p\\ summing to the object's `size`.

## See also

[`variance.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.MultinomialDistrib.md)
for the covariance,
[`mv_location.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.MultinomialDistrib.md),
which this calls, and
[`mean.distrib()`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md)
for the generic's default.
