# Gamma Analytical Observed Hessian

Computes the analytical observed Hessian (second derivatives) of the
Gamma log-density with respect to the parameters \\\mu\\ and
\\\sigma^2\\.

## Arguments

- distrib:

  A `GammaDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

## Value

A list containing the vectors of second derivatives.

## See also

[`gamma_distrib`](https://statmodels7.github.io/distributions7/reference/gamma_distrib.md)
