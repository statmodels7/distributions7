# Logistic Response Second Derivative

Closed-form \\\partial^2 \ell / \partial y^2 =
-\mathrm{sech}^2(z/2)/(2\sigma^2)\\, \\z = (y-\mu)/\sigma\\.

## Arguments

- distrib:

  A `LogisticDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A numeric vector.

## See also

[`logistic_distrib`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
