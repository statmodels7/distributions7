# Beta Analytical Gradient

Computes the analytical gradient (first derivatives) of the Beta
log-density with respect to the parameters \\\mu\\ and \\\phi\\.

\$\$\dfrac{\partial \ell}{\partial \mu} = \phi \left\[
\log\left(\dfrac{y}{1-y}\right) - \psi(\mu\phi) + \psi((1-\mu)\phi)
\right\]\$\$ \$\$\dfrac{\partial \ell}{\partial \phi} = \psi(\phi) -
\mu\psi(\mu\phi) - (1-\mu)\psi((1-\mu)\phi) + \mu \log(y) + (1-\mu)
\log(1-y)\$\$

## Arguments

- distrib:

  A `Beta1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `phi`.

## Value

A list containing the vectors of first derivatives.

## See also

[`beta1_distrib`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
