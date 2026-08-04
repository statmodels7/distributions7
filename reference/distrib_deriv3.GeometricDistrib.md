# Geometric Analytical Third-Order Derivative

\$\$\ell^{(\mu\mu\mu)} = 2\left(\dfrac{y}{\mu^3} -
\dfrac{y+1}{(1+\mu)^3}\right), \qquad \mathbb{E}\[\ell^{(\mu\mu\mu)}\] =
2\left(\dfrac{1}{\mu^2} - \dfrac{1}{(1+\mu)^2}\right)\$\$

## Arguments

- distrib:

  A `GeometricDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- expected:

  Logical; if `TRUE`, returns the expected derivative.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list with the `mu_mu_mu` component.

## See also

[`geometric_distrib`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
