# Mean of a Dirichlet

Returns the mean vector \\\mathbb{E}\[Y\] = \mu\\, which for this family
is a parameter and needs no computation: the method delegates to
[`mv_location.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.DirichletDistrib.md).
The result sums to one, being a point of the simplex.

## Arguments

- x:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).
  The argument is named `x` because the generic is
  [`base::mean()`](https://rdrr.io/r/base/mean.html).

- theta:

  A named list of parameters on the parameter scale: the mean's free
  values followed by `phi`, each of length 1.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length \\p\\, strictly positive and summing to one.

## See also

[`variance.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.DirichletDistrib.md)
for the covariance,
[`mv_location.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.DirichletDistrib.md),
which this calls, and
[`mean.distrib()`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md)
for the generic's default.
