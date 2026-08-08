# Laplace Analytical Gradient in Location and Rate

The derivative with respect to \\\mu\\ exists only almost everywhere
(there is a kink at \\y = \mu\\, a set of probability zero); the
subgradient value 0 is returned there.

\$\$\dfrac{\partial \ell}{\partial \mu} = \lambda\\\mathrm{sign}(y-\mu),
\qquad \dfrac{\partial \ell}{\partial \lambda} = \dfrac{1}{\lambda} -
\|y-\mu\|\$\$

## Arguments

- distrib:

  A `Laplace2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `lambda`.

## Value

A list containing the vectors of first derivatives.

## See also

[`laplace2_distrib`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
