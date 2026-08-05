# Gaussian Density in Mean and Variance

\$\$f(y) = \dfrac{1}{\sqrt{2\pi\sigma^2}}
\exp\left\\-\dfrac{(y-\mu)^2}{2\sigma^2}\right\\\$\$

## Arguments

- distrib:

  A `Gaussian2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `sigma2`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector.

## See also

[`gaussian2_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
