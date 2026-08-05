# Gaussian Analytical Expected Hessian

Computes the analytical expected Hessian of the Gaussian log-density
with respect to the parameters \\\mu\\ and \\\sigma\\.

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\sigma^2}\$\$ \$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \sigma^2}\right\] = -\dfrac{2}{\sigma^2}\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu \partial
\sigma}\right\] = 0\$\$

## Arguments

- distrib:

  A `Gaussian1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
