# Cauchy Analytical Gradient

Computes the analytical gradient (first derivatives) of the Cauchy
log-density with respect to the parameters \\\mu\\ and \\\sigma\\.

\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{2(y-\mu)}{\sigma^2 +
(y-\mu)^2}\$\$ \$\$\dfrac{\partial \ell}{\partial \sigma} =
\dfrac{(y-\mu)^2 - \sigma^2}{\sigma(\sigma^2 + (y-\mu)^2)}\$\$

## Arguments

- distrib:

  A `CauchyDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A list containing the vectors of first derivatives.

## See also

[`cauchy_distrib`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
