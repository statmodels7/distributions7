# von Mises Density

\$\$f(y) = \dfrac{e^{\kappa \cos(y - \mu)}}{2\pi I_0(\kappa)}, \qquad y
\in \[-\pi, \pi)\$\$

## Arguments

- distrib:

  A `VonMises1Distrib` object.

- y:

  A numeric vector of angles.

- theta:

  A list containing `mu` and `kappa`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`vonmises1_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
