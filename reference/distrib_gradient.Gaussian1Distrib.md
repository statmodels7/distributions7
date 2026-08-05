# Gaussian Analytical Gradient

Computes the analytical gradient (first derivatives) of the Gaussian
log-density with respect to the parameters \\\mu\\ and \\\sigma\\.

\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\sigma^2}\$\$
\$\$\dfrac{\partial \ell}{\partial \sigma} = \dfrac{(y - \mu)^2 -
\sigma^2}{\sigma^3}\$\$

## Arguments

- distrib:

  A `Gaussian1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A list containing the vectors of first derivatives.

## See also

[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
