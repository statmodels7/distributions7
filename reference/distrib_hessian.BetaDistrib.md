# Beta Analytical Observed Hessian

Computes the analytical observed Hessian (second derivatives) of the
Beta log-density with respect to the parameters \\\mu\\ and \\\phi\\.

## Arguments

- distrib:

  A `BetaDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `phi`.

## Value

A list containing the vectors of second derivatives.

## See also

[`beta_distrib`](https://statmodels7.github.io/distributions7/reference/beta_distrib.md)
