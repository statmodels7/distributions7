# Beta Response Derivatives

Closed-form derivatives of the Beta log-density with respect to the
response, with \\\alpha = \mu\phi\\, \\\beta = (1-\mu)\phi\\: \\\partial
\ell / \partial y = (\alpha-1)/y - (\beta-1)/(1-y)\\ and \\\partial^2
\ell / \partial y^2 = -(\alpha-1)/y^2 - (\beta-1)/(1-y)^2\\.

## Arguments

- distrib:

  A `Beta1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `phi`.

## Value

A numeric vector.

## See also

[`beta1_distrib`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
