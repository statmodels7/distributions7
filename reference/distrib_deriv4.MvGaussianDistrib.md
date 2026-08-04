# Multivariate Gaussian Fourth Derivatives

Closed form; the same algebra as
[`distrib_deriv3.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MvGaussianDistrib.md)
one order up.

## Arguments

- distrib:

  A
  [`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- expected:

  Logical; the expectation is approximated by sampling.

- approx:

  Strategy label; sampling is the only multivariate route.

- nsim:

  Monte Carlo sample size.

## Value

A named list of fourth-derivative component vectors.
