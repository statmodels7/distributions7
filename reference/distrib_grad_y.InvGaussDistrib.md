# Inverse-Gaussian Response Derivatives

Closed-form derivatives of the Inverse-Gaussian log-density with respect
to the response: \\\partial \ell / \partial y = -\dfrac{3}{2y} -
\dfrac{y^2 - \mu^2}{2\phi\mu^2 y^2}\\ and \\\partial^2 \ell / \partial
y^2 = \dfrac{3}{2y^2} - \dfrac{1}{\phi y^3}\\.

## Arguments

- distrib:

  An `InvGaussDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `phi`.

## Value

A numeric vector.

## See also

[`invgauss_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss_distrib.md)
