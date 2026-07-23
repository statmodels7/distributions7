# Negative Binomial Analytical Observed Hessian

Computes the analytical observed Hessian (second derivatives) of the
Negative Binomial log-probability with respect to the parameters \\\mu\\
and \\\theta\\.

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
\dfrac{y+\theta}{(\theta+\mu)^2} - \dfrac{y}{\mu^2}\$\$
\$\$\dfrac{\partial^2 \ell}{\partial \theta^2} = \psi_1(y+\theta) -
\psi_1(\theta) + \dfrac{\mu}{\theta(\theta+\mu)} +
\dfrac{y-\mu}{(\theta+\mu)^2}\$\$ \$\$\dfrac{\partial^2 \ell}{\partial
\mu \partial \theta} = \dfrac{y-\mu}{(\theta+\mu)^2}\$\$

## Arguments

- distrib:

  A `NegBinDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `theta`.

## Value

A list containing the vectors of second derivatives.

## See also

[`negbin_distrib`](https://statmodels7.github.io/distributions7/reference/negbin_distrib.md)
