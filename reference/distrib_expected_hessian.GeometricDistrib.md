# Geometric Analytical Expected Hessian

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\mu(1+\mu)}\$\$ the reciprocal of the variance, as it must be
for a family written in its mean.

## Arguments

- distrib:

  A `GeometricDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list with the `mu_mu` component.

## See also

[`geometric_distrib`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
