# Pseudo-Huber Probability Density Function

Computes the probability density function for the Pseudo-Huber
distribution: \$\$f(y; \mu, \sigma, \nu) = \dfrac{1}{2 \sigma \sqrt{\nu}
K_1(\sqrt{\nu})} \exp\left( - \sqrt{\nu +
\left(\dfrac{y-\mu}{\sigma}\right)^2} \right)\$\$ where \\K_1\\ is the
modified Bessel function of the second kind.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma` and `nu`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`pseudohuber_distrib`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
