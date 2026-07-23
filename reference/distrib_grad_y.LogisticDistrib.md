# Logistic Response Derivatives

Closed-form derivatives of the Logistic log-density with respect to the
response. With \\z = (y-\mu)/\sigma\\: \\\partial \ell / \partial y =
-\tanh(z/2)/\sigma\\ and \\\partial^2 \ell / \partial y^2 =
-\mathrm{sech}^2(z/2)/(2\sigma^2)\\.

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
