# Student's t Mixed Derivatives

Closed form: with \\r = y - \mu\\ and \\D = \nu\sigma^2 + r^2\\,
\\\partial^2 \ell / \partial y\\ \partial \mu = (\nu+1)(\nu\sigma^2 -
r^2)/D^2\\, \\\partial^2 \ell / \partial y\\ \partial \sigma =
2\nu\sigma(\nu+1)\\ r/D^2\\, \\\partial^2 \ell / \partial y\\ \partial
\nu = -r\\(r^2 - \sigma^2)/D^2\\.

## Arguments

- distrib:

  A `StudentT1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma` and `nu`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `mu`, `sigma` and `nu`.
