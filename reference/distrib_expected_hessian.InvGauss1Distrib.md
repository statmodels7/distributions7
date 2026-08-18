# Inverse-Gaussian Analytical Expected Hessian

Computes the analytical expected Hessian of the Inverse-Gaussian
log-density with respect to the parameters \\\mu\\ and \\\phi\\.

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\phi\mu^3}\$\$ \$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \phi^2}\right\] = -\dfrac{1}{2\phi^2}\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu \partial
\phi}\right\] = 0\$\$

## Arguments

- distrib:

  An `InvGauss1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `phi`.

- threads:

  How many threads the kernel may use; below the measured internal
  threshold it stays sequential whatever the count says.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`invgauss1_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
