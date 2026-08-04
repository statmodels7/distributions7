# Exponential Analytical Expected Hessian

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\mu^2}\$\$ obtained from the observed form by
\\\mathbb{E}\[y\] = \mu\\.

## Arguments

- distrib:

  An `ExponentialDistrib` object.

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

[`exponential_distrib`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
