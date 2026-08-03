# Gaussian Mixed Derivatives

Closed form: with \\r = y - \mu\\, \\\partial^2 \ell / \partial y\\
\partial \mu = 1/\sigma^2\\ and \\\partial^2 \ell / \partial y\\
\partial \sigma = 2r/\sigma^3\\.

## Arguments

- distrib:

  A `GaussianDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `mu` and `sigma`.
