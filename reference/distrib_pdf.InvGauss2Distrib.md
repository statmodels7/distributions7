# Inverse Gaussian Density in Mean and Shape

\$\$f(y) = \sqrt{\dfrac{\lambda}{2\pi y^3}}
\exp\left\\-\dfrac{\lambda(y-\mu)^2}{2\mu^2 y}\right\\\$\$

## Arguments

- distrib:

  An `InvGauss2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `lambda`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector.

## See also

[`invgauss2_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md)
