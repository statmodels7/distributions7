# Inverse Gaussian Analytical Expected Hessian in Mean and Shape

The response enters every derivative linearly, so the expectations need
only \\\mathbb{E}\[Y\] = \mu\\: \$\$\mathbb{E}\[\ell^{(\mu\mu)}\] =
-\dfrac{\lambda}{\mu^3}, \qquad \mathbb{E}\[\ell^{(\mu\lambda)}\] = 0,
\qquad \mathbb{E}\[\ell^{(\lambda\lambda)}\] =
-\dfrac{1}{2\lambda^2}\$\$ The mean and the shape are orthogonal.

## Arguments

- distrib:

  An `InvGauss2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `lambda`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of expected second derivatives.

## See also

[`invgauss2_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md)
