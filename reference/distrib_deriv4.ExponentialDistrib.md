# Exponential Analytical Fourth-Order Derivative

\$\$\ell^{(\mu\mu\mu\mu)} = \dfrac{6}{\mu^4} - \dfrac{24y}{\mu^5},
\qquad \mathbb{E}\[\ell^{(\mu\mu\mu\mu)}\] = -\dfrac{18}{\mu^4}\$\$

## Arguments

- distrib:

  An `ExponentialDistrib` object.

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

A named list with the `mu_mu_mu_mu` component.

## See also

[`exponential_distrib`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
