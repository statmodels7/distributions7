# Multivariate Gaussian Score

Closed form. With \\w = \Sigma^{-1}(y - \mu)\\ and \\A_k\\ the
derivative of \\\Sigma\\ in the \\k\\-th free value of the matrix
parameter, \$\$\frac{\partial \ell}{\partial \mu} = w, \qquad
\frac{\partial \ell}{\partial \eta_k} = -\frac{1}{2}\frac{\partial
\log\|\Sigma\|}{\partial \eta_k} + \frac{1}{2} w^\top A_k w.\$\$ The
first term of the second expression is the matrix parameter's own
`param_dlogdet()`, so no trace is formed here.

## Arguments

- distrib:

  A
  [`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- scale:

  Handled by the generic; the two scales coincide here.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
