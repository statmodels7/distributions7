# Variance of a Dirichlet

Returns the covariance matrix of the family by delegating to
[`mv_sigma.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.DirichletDistrib.md).
It is singular, the coordinates summing to one, so a caller wanting
something to invert should take the marginals or work on the simplex's
free vector.

## Arguments

- x:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).
  The argument is named `x` for consistency with
  [`base::mean()`](https://rdrr.io/r/base/mean.html), whose signature
  the moment generics follow.

- theta:

  A named list of parameters on the parameter scale: the mean's free
  values followed by `phi`, each of length 1.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A symmetric \\p \times p\\ numeric matrix of rank \\p-1\\.

## See also

[`mean.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.DirichletDistrib.md)
for the mean vector,
[`mv_sigma.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.DirichletDistrib.md),
which this calls, and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
for the generic.
