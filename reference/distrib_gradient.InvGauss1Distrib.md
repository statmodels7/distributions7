# Inverse-Gaussian Analytical Gradient

Computes the analytical gradient (first derivatives) of the
Inverse-Gaussian log-density with respect to the parameters \\\mu\\ and
\\\phi\\.

\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\phi\mu^3}\$\$
\$\$\dfrac{\partial \ell}{\partial \phi} = \dfrac{(y - \mu)^2 -
y\mu^2\phi}{2y\phi^2\mu^2}\$\$

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

A list containing the vectors of first derivatives.

## See also

[`invgauss1_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
