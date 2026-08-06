# Inverse Gaussian Mixed Derivatives

Closed form. With \\\ell^{(y)} = -3/(2y) - 1/(2\phi\mu^2) + 1/(2\phi
y^2)\\, the components are \\1/(\phi\mu^3)\\ and \\1/(2\phi^2\mu^2) -
1/(2\phi^2 y^2)\\.

## Arguments

- distrib:

  An `InvGauss1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `phi`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
