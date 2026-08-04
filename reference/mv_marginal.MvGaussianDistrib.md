# Marginal of a Multivariate Gaussian

A marginal of a gaussian is a gaussian: the mean is the subvector and
the covariance the corresponding block, with no integration to perform.

## Arguments

- distrib:

  A
  [`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- theta:

  A named list of parameters.

- which:

  An integer vector of coordinates.

- ...:

  Unused.

## Value

A list with `distrib` and `theta`.

## Details

The marginal is returned on an unstructured covariance whatever the
parent carried, because a block of a structured matrix need not have the
parent's structure – the leading block of an AR(1) is AR(1), but a block
of a compound-symmetry matrix taken at scattered indices need not be,
and a precision block is not the inverse of the covariance block at all.
