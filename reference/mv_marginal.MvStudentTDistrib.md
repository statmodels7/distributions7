# Marginal of a Multivariate Student t

A marginal of a \\t\\ is a \\t\\ with the same degrees of freedom, the
subvector of the location and the corresponding block of the scale
matrix. The degrees of freedom do not change with the dimension, which
is what makes the family closed under marginalisation at all.

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- theta:

  A named list of parameters.

- which:

  An integer vector of coordinates.

- ...:

  Unused.

## Value

A list with `distrib` and `theta`.
