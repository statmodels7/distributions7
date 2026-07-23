# Beta Analytical Expected Hessian

Computes the analytical expected Hessian of the Beta log-density with
respect to the parameters \\\mu\\ and \\\phi\\.

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\phi^2 \left\[ \psi_1(\mu\phi) + \psi_1((1-\mu)\phi) \right\]\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \phi^2}\right\] =
\psi_1(\phi) - \mu^2\psi_1(\mu\phi) - (1-\mu)^2\psi_1((1-\mu)\phi)\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu \partial
\phi}\right\] = -\phi \left\[ \mu\psi_1(\mu\phi) -
(1-\mu)\psi_1((1-\mu)\phi) \right\]\$\$

## Arguments

- distrib:

  A `BetaDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `phi`.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`beta_distrib`](https://statmodels7.github.io/distributions7/reference/beta_distrib.md)
