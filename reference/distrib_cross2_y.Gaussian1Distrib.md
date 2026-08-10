# Gaussian Mixed Second-Response Derivatives

Closed form: \\\ell^{(yy)} = -1/\sigma^2\\ does not depend on the
location, so \\\partial^3\ell/\partial y^2\partial\mu = 0\\ and
\\\partial^3\ell/\partial y^2\partial\sigma = 2/\sigma^3\\.

## Arguments

- distrib:

  A `Gaussian1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
