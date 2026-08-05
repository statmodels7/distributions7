# Lognormal Response Second Derivative

Closed-form \\\partial^2 \ell / \partial y^2 = (1 + (\log y - \mu -
1)/\sigma^2)/y^2\\.

## Arguments

- distrib:

  A `Lognormal1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

## Value

A numeric vector.

## See also

[`lognormal1_distrib`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
