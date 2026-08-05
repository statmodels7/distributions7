# Dirichlet Analytical Gradient

With \\g_j = \log y_j - \psi(\alpha_j)\\ and \\A =
\partial\mu/\partial\eta\\, \$\$\dfrac{\partial\ell}{\partial\eta_k} =
\phi \sum_j g_j A\_{jk}, \qquad \dfrac{\partial\ell}{\partial\phi} =
\psi(\phi) + \sum_j g_j \mu_j\$\$

## Arguments

- distrib:

  A `DirichletDistrib` object.

- y:

  A matrix with one row per observation.

- theta:

  A named list of parameters.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list, one component per parameter.

## See also

[`dirichlet_distrib`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
