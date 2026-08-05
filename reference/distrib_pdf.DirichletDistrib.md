# Dirichlet Density

\$\$f(y) = \dfrac{\Gamma(\phi)}{\prod_j \Gamma(\alpha_j)} \prod_j
y_j^{\alpha_j - 1}, \qquad \alpha = \phi\mu\$\$ evaluated on the open
simplex, one row of `y` per observation.

## Arguments

- distrib:

  A `DirichletDistrib` object.

- y:

  A matrix with one row per observation, each summing to one.

- theta:

  A named list of parameters.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector, one entry per row of `y`.

## See also

[`dirichlet_distrib`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
