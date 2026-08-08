# Laplace Fourth-Order Derivatives in Location and Rate

Closed form, almost everywhere. The only non-zero component is
\\\ell^{(\lambda\lambda\lambda\lambda)} = -6/\lambda^4\\; the observed
and the expected derivatives coincide.

## Arguments

- distrib:

  A `Laplace2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `lambda`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivatives.

## Value

A named list of fourth-derivative component vectors.

## See also

[`laplace2_distrib`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
