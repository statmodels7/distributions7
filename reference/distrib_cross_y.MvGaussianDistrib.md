# Multivariate Gaussian Mixed Response-Parameter Derivatives

\\\partial^2 \ell / \partial y \partial \theta_k\\, one \\n \times p\\
matrix per parameter. With \\w = \Sigma^{-1}(y - \mu)\\ the response
gradient is \\-w\\, so differentiating it in the mean and in the free
values of the matrix parameter gives \$\$\frac{\partial^2 \ell}{\partial
y \partial \mu_j} = \Sigma^{-1}e_j, \qquad \frac{\partial^2
\ell}{\partial y \partial \eta_k} = \Sigma^{-1}A_k w,\$\$ with \\A_k =
\partial\Sigma/\partial\eta_k\\. The mean block is the same at every
observation, the matrix block is not.

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

A named list, one \\n \times p\\ matrix per parameter.

## Details

The shape is the one a consumer needs: a penalty whose prior is this
family reads the block of
\\\partial^2\rho/\partial\beta\\\partial\theta_k\\ for one
hyperparameter at a time, and the coefficients of one group are the row
of \\y\\ the density is read at.

The link scale is the parameter scale here: the mean components carry
the identity link and so do the matrix parameter's free values, which
are unconstrained by construction.
