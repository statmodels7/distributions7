# Skew Normal Analytical Gradient

Closed-form first derivatives, with \\z = (y-\mu)/\sigma\\, \\t = \alpha
z\\ and \\R = \phi(t)/\Phi(t)\\: \$\$\dfrac{\partial \ell}{\partial \mu}
= \dfrac{z - \alpha R}{\sigma}, \qquad \dfrac{\partial \ell}{\partial
\sigma} = \dfrac{z^2 - 1 - \alpha z R}{\sigma}, \qquad \dfrac{\partial
\ell}{\partial \alpha} = z R.\$\$

## Arguments

- distrib:

  A `SkewNormal1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma` and `alpha`.

- scale:

  Either `"parameter"` or `"link"`.

- ...:

  Unused.

## Value

A named list of first derivatives.

## See also

[`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
