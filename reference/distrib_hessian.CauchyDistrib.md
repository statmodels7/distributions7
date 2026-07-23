# Cauchy Analytical Observed Hessian

Computes the analytical observed Hessian (second derivatives) of the
Cauchy log-density with respect to the parameters \\\mu\\ and
\\\sigma\\.

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{2(y-\mu)^2 -
2\sigma^2}{(\sigma^2 + (y-\mu)^2)^2}\$\$ \$\$\dfrac{\partial^2
\ell}{\partial \sigma^2} = \dfrac{\sigma^4 - 4\sigma^2 (y-\mu)^2 -
(y-\mu)^4}{\sigma^2(\sigma^2 + (y-\mu)^2)^2}\$\$ \$\$\dfrac{\partial^2
\ell}{\partial \mu \partial \sigma} = -\dfrac{4\sigma
(y-\mu)}{(\sigma^2 + (y-\mu)^2)^2}\$\$

## Arguments

- distrib:

  A `CauchyDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A list containing the vectors of second derivatives.

## See also

[`cauchy_distrib`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
