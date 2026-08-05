# Beta Analytical Observed Hessian

Computes the analytical observed Hessian (second derivatives) of the
Beta log-density with respect to the parameters \\\mu\\ and \\\phi\\.

## Arguments

- distrib:

  A `Beta1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `phi`.

## Value

A list containing the vectors of second derivatives.

## See also

[`beta1_distrib`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
