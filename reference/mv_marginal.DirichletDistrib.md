# Dirichlet Marginal

A single coordinate is \\\mathrm{Beta}(\alpha_j, \phi-\alpha_j)\\, which
in this package's mean-and-precision parametrization of the Beta is
simply mean \\\mu_j\\ and precision \\\phi\\: the concentration is
shared by every marginal. A group of coordinates is again Dirichlet, but
only after the remaining mass is collapsed into one of its own, so that
case is rejected.

## Arguments

- distrib:

  A `DirichletDistrib` object.

- theta:

  A named list of parameters.

- which:

  The coordinate wanted.

- ...:

  Unused.

## Value

A list with the marginal `distrib` and its `theta`.

## See also

[`dirichlet_distrib`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
