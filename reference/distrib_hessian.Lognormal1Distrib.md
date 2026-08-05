# Lognormal Analytical Observed Hessian

Computes the analytical observed Hessian (second derivatives) of the
Lognormal log-density with respect to the parameters \\\mu\\ and
\\\sigma^2\\.

## Arguments

- distrib:

  A `Lognormal1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

## Value

A list containing the vectors of second derivatives.

## See also

[`lognormal1_distrib`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
