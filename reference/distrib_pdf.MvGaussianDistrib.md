# Multivariate Gaussian Density

\$\$\ell = -\frac{p}{2}\log 2\pi - \frac{1}{2}\log\|\Sigma\| -
\frac{1}{2}(y-\mu)^\top \Sigma^{-1} (y-\mu),\$\$ evaluated row by row.
The quadratic form goes through the matrix parameter's own factor rather
than through an explicit inverse.

## Arguments

- distrib:

  A
  [`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations, or a vector of length \\p\\
  for one observation.

- theta:

  A named list of parameters.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector with one value per observation.
