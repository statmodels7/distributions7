# Skew Normal Analytical Observed Hessian

Closed-form second derivatives, with \\R' = -R(t + R)\\:
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{\alpha^2 R' -
1}{\sigma^2}, \qquad \dfrac{\partial^2 \ell}{\partial \mu \\ \partial
\sigma} = \dfrac{\alpha^2 z R' - 2z + \alpha R}{\sigma^2}, \qquad
\dfrac{\partial^2 \ell}{\partial \mu \\ \partial \alpha} = -\dfrac{R +
\alpha z R'}{\sigma},\$\$ \$\$\dfrac{\partial^2 \ell}{\partial \sigma^2}
= \dfrac{1 - 3z^2 + 2\alpha z R + \alpha^2 z^2 R'}{\sigma^2}, \qquad
\dfrac{\partial^2 \ell}{\partial \sigma \\ \partial \alpha} = -\dfrac{z
R + \alpha z^2 R'}{\sigma}, \qquad \dfrac{\partial^2 \ell}{\partial
\alpha^2} = z^2 R'.\$\$

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

A named list of second derivatives.

## See also

[`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
