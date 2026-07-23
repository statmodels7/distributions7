# Gamma Response Derivatives

Closed-form derivatives of the Gamma log-density with respect to the
response, using the shape/rate reparameterization \\\alpha =
\mu^2/\sigma^2\\, \\\lambda = \mu/\sigma^2\\: \\\partial \ell / \partial
y = (\alpha-1)/y - \lambda\\ and \\\partial^2 \ell / \partial y^2 =
-(\alpha-1)/y^2\\.

## Arguments

- distrib:

  A `GammaDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

## Value

A numeric vector.

## See also

[`gamma_distrib`](https://statmodels7.github.io/distributions7/reference/gamma_distrib.md)
