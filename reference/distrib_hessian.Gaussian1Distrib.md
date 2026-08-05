# Gaussian Analytical Observed Hessian

Computes the analytical observed Hessian (second derivatives) of the
Gaussian log-density with respect to the parameters \\\mu\\ and
\\\sigma\\.

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{\sigma^2}\$\$
\$\$\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^2 - 3(y -
\mu)^2}{\sigma^4}\$\$ \$\$\dfrac{\partial^2 \ell}{\partial \mu \partial
\sigma} = -\dfrac{2(y - \mu)}{\sigma^3}\$\$

## Arguments

- distrib:

  A `Gaussian1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A list containing the vectors of second derivatives.

## See also

[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
