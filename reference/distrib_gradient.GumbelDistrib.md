# Gumbel Analytical Gradient

Closed-form first derivatives of the Gumbel log-density, written in \\z
= (y-\mu)/\sigma\\ and \\w = e^{-z}\\: \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{1 - w}{\sigma}, \qquad \dfrac{\partial
\ell}{\partial \sigma} = \dfrac{z(1 - w) - 1}{\sigma}\$\$

## Arguments

- distrib:

  A `GumbelDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- scale:

  Either `"parameter"` or `"link"`.

- ...:

  Unused.

## Value

A named list of first derivatives.

## See also

[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
