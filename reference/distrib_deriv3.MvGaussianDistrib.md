# Multivariate Gaussian Third Derivatives

Closed form, built on the structure's own third derivatives from
parameters7. A component with three mean indices vanishes, the quadratic
form being quadratic; one mean index gives \\(P\_{klm} r)\_i\\; two give
\\-P\_{kl}\[i, j\]\\; none gives \\\mp\tfrac{1}{2}\\\partial^3
\log\|M\| - \tfrac{1}{2} r^\top P\_{klm} r\\. The precision's derivative
tensors \\P_t\\ come directly from the structure under a precision
parametrisation, and by the Leibniz recursion on \\P_k = -P A_k P\\
under a covariance one.

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

  Logical; the expectation is approximated by sampling, as for the
  Hessian.

- approx:

  Strategy label; sampling is the only multivariate route.

- nsim:

  Monte Carlo sample size.

## Value

A named list of third-derivative component vectors.
