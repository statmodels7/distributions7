# Multivariate Student t Mixed Response-Parameter Derivatives

\\\partial^2 \ell / \partial y \partial \theta_k\\, one \\n \times p\\
matrix per parameter. The response gradient is \\-c\\w\\, so every
component carries the derivative of the weight beside the gaussian term
it multiplies: with \\A_k = \partial\Sigma/\partial\eta_k\\,
\$\$\frac{\partial^2\ell}{\partial y\\\partial\mu_j} =
c\\\Sigma^{-1}e_j - \frac{\partial c}{\partial\mu_j}\\w, \qquad
\frac{\partial c}{\partial\mu_j} = \frac{2(\nu+p)\\w_j}{(\nu+q)^2},\$\$
\$\$\frac{\partial^2\ell}{\partial y\\\partial\eta_k} =
c\\\Sigma^{-1}A_k w - \frac{\partial c}{\partial\eta_k}\\w, \qquad
\frac{\partial c}{\partial\eta_k} = \frac{(\nu+p)\\w^\top A_k
w}{(\nu+q)^2},\$\$ \$\$\frac{\partial^2\ell}{\partial y\\\partial\nu} =
-\frac{(q-p)\\w}{(\nu+q)^2}.\$\$

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
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

Nothing here is obstructed: the log-density carries no distribution
function, only `lgamma`, a logarithm and a quadratic form, each
elementary in \\\nu\\. As \\\nu \to \infty\\ the weight and its
derivatives go to one and to zero, and every component becomes the
gaussian's, which is what the tests compare it against.
