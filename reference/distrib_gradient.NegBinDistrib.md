# Negative Binomial Analytical Gradient

Computes the analytical gradient (first derivatives) of the Negative
Binomial log-probability with respect to the parameters \\\mu\\ and
\\\theta\\.

\$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{\theta}{\theta+\mu}\left(\dfrac{y}{\mu} - 1\right)\$\$
\$\$\dfrac{\partial \ell}{\partial \theta} = \psi(y+\theta) -
\psi(\theta) + \log\left(\dfrac{\theta}{\theta+\mu}\right) +
\dfrac{\mu - y}{\theta+\mu}\$\$

## Arguments

- distrib:

  A `NegBinDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `theta`.

## Value

A list containing the vectors of first derivatives.

## See also

[`negbin_distrib`](https://statmodels7.github.io/distributions7/reference/negbin_distrib.md)
