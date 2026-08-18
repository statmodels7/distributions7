# Inverse-Gaussian Analytical Observed Hessian

Computes the analytical observed Hessian (second derivatives) of the
Inverse-Gaussian log-density with respect to the parameters \\\mu\\ and
\\\phi\\.

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{3y -
2\mu}{\phi\mu^4}\$\$ \$\$\dfrac{\partial^2 \ell}{\partial \phi^2} =
\dfrac{\phi - 2(y-\mu)^2/(\mu^2 y)}{2\phi^3}\$\$ \$\$\dfrac{\partial^2
\ell}{\partial \mu \partial \phi} = -\dfrac{y - \mu}{\phi^2\mu^3}\$\$

**Note:** The observed Hessian with respect to \\\phi\\ is not
guaranteed to be negative for all observed values of \\y\\.
Specifically, \\\partial^2 \ell/\partial \phi^2 \< 0\\ only when \\\phi
\< 2(y-\mu)^2/(\mu^2 y)\\. This condition may be violated when
observations are far from the mean or when the dispersion parameter is
large, potentially causing numerical instability in optimization
algorithms that rely on the observed Hessian (e.g., Newton-Raphson). In
such cases, using the expected Hessian (`distrib_expected_hessian`) is
recommended for more stable convergence.

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

A list containing the vectors of second derivatives.

## See also

[`invgauss1_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
