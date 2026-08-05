# Generalised Pareto Density

\$\$f(y) = \dfrac{1}{\sigma} \left(1 + \dfrac{\xi
y}{\sigma}\right)^{-1/\xi - 1}\$\$ with the exponential density as the
limit at \\\xi = 0\\.

## Arguments

- distrib:

  A `GPDDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `sigma` and `xi`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
