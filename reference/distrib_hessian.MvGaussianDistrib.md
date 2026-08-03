# Multivariate Gaussian Observed Hessian

Closed form. With \\w = \Sigma^{-1}(y-\mu)\\, \\A_k\\ and \\A\_{kl}\\
the first and second derivatives of \\\Sigma\\, \$\$\ell^{(\mu_a \mu_b)}
= -(\Sigma^{-1})\_{ab}, \qquad \ell^{(\mu_a \eta_k)} = -(\Sigma^{-1} A_k
w)\_a,\$\$ \$\$\ell^{(\eta_k \eta_l)} = -\tfrac{1}{2}\frac{\partial^2
\log\|\Sigma\|}{\partial\eta_k \partial\eta_l} + \tfrac{1}{2} w^\top
A\_{kl} w - w^\top A_l \Sigma^{-1} A_k w.\$\$

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

A named list keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)`(distrib@params)`.
