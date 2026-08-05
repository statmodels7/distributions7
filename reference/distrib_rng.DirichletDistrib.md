# Dirichlet Random Generation

Independent Gamma draws with the shapes \\\alpha_j\\, normalized by
their sum, which is the representation the family is defined by.

## Arguments

- distrib:

  A `DirichletDistrib` object.

- n:

  The number of draws.

- theta:

  A named list of parameters.

## Value

A matrix with `n` rows, each summing to one.

## See also

[`dirichlet_distrib`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
