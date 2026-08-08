# Laplace Response Derivative in Location and Rate

\\\partial \ell / \partial y = -\lambda\\\mathrm{sign}(y-\mu)\\, almost
everywhere; the analytic form is provided because finite differences
would be inaccurate across the kink at \\y = \mu\\.

## Arguments

- distrib:

  A `Laplace2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `lambda`.

## Value

A numeric vector.

## See also

[`laplace2_distrib`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
