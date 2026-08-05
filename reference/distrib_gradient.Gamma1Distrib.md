# Gamma Analytical Gradient in Mean and Dispersion

With \\s = 1/\phi\\ and \\z = y/\mu\\,
\$\$\dfrac{\partial\ell}{\partial\mu} = \dfrac{y-\mu}{\phi\mu^2}, \qquad
\dfrac{\partial\ell}{\partial\phi} = -s^2\left\\\log s + 1 - \psi(s) +
\log z - z\right\\\$\$ the first being the score of a gamma generalized
linear model.

## Arguments

- distrib:

  A `Gamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `phi`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of first derivatives.

## See also

[`gamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
