# Elastic-Net Density

\$\$f(y; \mu, \lambda, \alpha) = \frac{1}{Z} \exp\left\\-a\|y-\mu\| -
\tfrac{c}{2}(y-\mu)^2\right\\,\$\$ with \\a = \lambda\alpha\\, \\c =
\lambda(1-\alpha)\\ and \\Z = 2M(a/\sqrt{c})/\sqrt{c}\\, where \\M\\ is
the Mills ratio. The constant is evaluated through the log Mills ratio,
both factors of which underflow separately.

## Arguments

- distrib:

  An `EnetDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `lambda` and `alpha`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
