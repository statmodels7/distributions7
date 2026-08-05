# Gamma Analytical Gradient

Computes the analytical gradient (first derivatives) of the Gamma
log-density with respect to the parameters \\\mu\\ and \\\sigma^2\\.

\$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{-2\mu\psi\left(\dfrac{\mu^2}{\sigma^2}\right) +
2\mu\log\left(\dfrac{\mu}{\sigma^2}\right) + \mu + 2\mu\log(y) -
y}{\sigma^2}\$\$ \$\$\dfrac{\partial \ell}{\partial \sigma^2} =
-\dfrac{\mu\left\[-\mu\psi\left(\dfrac{\mu^2}{\sigma^2}\right) + \mu +
\mu\left(\log\left(\dfrac{\mu}{\sigma^2}\right) + \log(y)\right) -
y\right\]}{(\sigma^2)^2}\$\$

## Arguments

- distrib:

  A `Gamma2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

## Value

A list containing the vectors of first derivatives.

## See also

[`gamma2_distrib`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md)
