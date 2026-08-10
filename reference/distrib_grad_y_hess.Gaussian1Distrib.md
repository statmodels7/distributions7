# Gaussian Second-Order Mixed Derivatives

Closed forms. With \\r = y - \mu\\, the response gradient is
\\-r/\sigma^2\\ and the response curvature \\-1/\sigma^2\\, so
\$\$\partial^3\ell/\partial y\\\partial\mu^2 = 0, \quad
\partial^3\ell/\partial y\\\partial\mu\partial\sigma = -2/\sigma^3,
\quad \partial^3\ell/\partial y\\\partial\sigma^2 = -6r/\sigma^4,\$\$
and the only component of the fourth derivative that does not vanish is
\\\partial^4\ell/\partial y^2\partial\sigma^2 = -6/\sigma^4\\, the
curvature carrying no location at all.

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

A named list keyed by parameter pair.
