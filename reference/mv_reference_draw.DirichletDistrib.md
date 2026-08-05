# A Uniform Proposal on the Simplex

Draws from the uniform distribution on the simplex, which is the
Dirichlet with every shape equal to one and has the constant density
\\\Gamma(p)\\ with respect to the same dominating measure the family's
density is written against. The base class's gaussian proposal lives in
\\\mathbb{R}^p\\ and would place no mass at all on the simplex.

## Arguments

- distrib:

  A `DirichletDistrib` object.

- theta:

  A named list of parameters.

- n:

  The number of draws.

- ...:

  Unused.

## Value

A list with the draws `y` and their log-density `logd`.

## See also

[`dirichlet_distrib`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
