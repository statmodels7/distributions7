# Multivariate Student t Density

\$\$\ell = \log\Gamma\\\left(\tfrac{\nu+p}{2}\right) -
\log\Gamma\\\left(\tfrac{\nu}{2}\right) - \tfrac{p}{2}\log(\nu\pi) -
\tfrac{1}{2}\log\|\Sigma\| - \tfrac{\nu+p}{2}\log\\\left(1 +
\tfrac{q}{\nu}\right),\$\$ with \\q = (y-\mu)^\top \Sigma^{-1}(y-\mu)\\.
The logarithm is taken with `log1p`, which is the difference between a
number and a loss of every significant digit when \\q/\nu\\ is small.

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector with one value per observation.
