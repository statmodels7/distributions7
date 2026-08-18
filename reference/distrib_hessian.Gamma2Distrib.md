# Gamma Analytical Observed Hessian

Computes the analytical observed Hessian (second derivatives) of the
Gamma log-density with respect to the parameters \\\mu\\ and
\\\sigma^2\\.

## Arguments

- distrib:

  A `Gamma2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

- threads:

  How many threads the kernel may use; below the measured internal
  threshold it stays sequential whatever the count says.

## Value

A list containing the vectors of second derivatives.

## See also

[`gamma2_distrib`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md)
