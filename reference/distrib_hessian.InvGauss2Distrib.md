# Inverse Gaussian Analytical Observed Hessian in Mean and Shape

\$\$\ell^{(\mu\mu)} = \dfrac{\lambda(2\mu-3y)}{\mu^4}, \qquad
\ell^{(\mu\lambda)} = \dfrac{y-\mu}{\mu^3}, \qquad
\ell^{(\lambda\lambda)} = -\dfrac{1}{2\lambda^2}\$\$

## Arguments

- distrib:

  An `InvGauss2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `lambda`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of second derivatives.

## See also

[`invgauss2_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md)
