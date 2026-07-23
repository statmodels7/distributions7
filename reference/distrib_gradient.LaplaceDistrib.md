# Laplace Analytical Gradient

Computes the analytical gradient of the Laplace log-density. Note that
the derivative with respect to \\\mu\\ exists only almost everywhere
(there is a kink at \\y = \mu\\, a set of probability zero); the
subgradient value 0 is returned there.

\$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{\mathrm{sign}(y-\mu)}{b}\$\$ \$\$\dfrac{\partial \ell}{\partial
b} = \dfrac{1}{b}\left(\dfrac{\|y-\mu\|}{b} - 1\right)\$\$

## Arguments

- distrib:

  A `LaplaceDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `b`.

## Value

A list containing the vectors of first derivatives.

## See also

[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
