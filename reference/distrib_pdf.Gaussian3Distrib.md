# Gaussian Density in Mean and Precision

\$\$f(y) = \sqrt{\dfrac{\tau}{2\pi}}
\exp\left\\-\dfrac{\tau(y-\mu)^2}{2}\right\\\$\$

## Arguments

- distrib:

  A `Gaussian3Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `tau`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector.

## See also

[`gaussian3_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
