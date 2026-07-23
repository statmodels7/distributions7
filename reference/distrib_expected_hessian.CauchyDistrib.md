# Cauchy Analytical Expected Hessian

Computes the analytical expected Hessian of the Cauchy log-density with
respect to the parameters \\\mu\\ and \\\sigma\\.

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{2\sigma^2}\$\$ \$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \sigma^2}\right\] = -\dfrac{1}{2\sigma^2}\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu \partial
\sigma}\right\] = 0\$\$

## Arguments

- distrib:

  A `CauchyDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`cauchy_distrib`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
