# Multivariate Gaussian Higher Mixed Response Derivatives

The three derivatives a marginal criterion reads when this family is a
prior. Writing \\B_k = \Sigma^{-1}A_k\Sigma^{-1}\\ with \\A_k =
\partial\Sigma/\partial\eta_k\\, \$\$\frac{\partial^3\ell}{\partial
y\\\partial y^\top\partial\eta_k} = B_k, \qquad
\frac{\partial^4\ell}{\partial y\\\partial y^\top
\partial\eta_k\partial\eta_l} = \Sigma^{-1}A\_{kl}\Sigma^{-1} -
\Sigma^{-1}\\\left(A_l\Sigma^{-1}A_k +
A_k\Sigma^{-1}A_l\right)\\\Sigma^{-1},\$\$
\$\$\frac{\partial^3\ell}{\partial y\\\partial\mu_j\partial\eta_k} =
-B_k e_j, \qquad \frac{\partial^3\ell}{\partial
y\\\partial\eta_k\partial\eta_l} = \Sigma^{-1}A\_{kl}w -
\Sigma^{-1}\\\left(A_l\Sigma^{-1}A_k + A_k\Sigma^{-1}A_l\right)\\w.\$\$

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

  Handled by the generic.

- ...:

  Unused.

## Value

A named list, keyed by parameter for `distrib_cross2_y` and by parameter
pair for the other two.

## Details

The response Hessian is \\-\Sigma^{-1}\\, which does not depend on the
observation and does not depend on the mean at all, so every component
of the first two involving a mean is exactly zero and the rest are one
matrix rather than one per row. Only the third carries the observation,
through \\w\\.
