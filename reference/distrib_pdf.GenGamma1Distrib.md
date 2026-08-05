# Generalised Gamma Density

\$\$f(y) = \dfrac{p}{a^{d}\\\Gamma(d/p)}\\ y^{d-1} e^{-(y/a)^{p}},
\qquad y \> 0\$\$

## Arguments

- distrib:

  A `GenGamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `a`, `d` and `p`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`gengamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
