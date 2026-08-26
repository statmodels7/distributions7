# A Uniform Proposal on the Simplex

Supplies the importance-sampling proposal
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
integrates the density against: the uniform distribution on the simplex,
which is the Dirichlet with every shape equal to one. Its density is the
constant \\\Gamma(p)\\ with respect to the same dominating measure the
family's own density is written against, so the normalization check is a
plain average.

The base class's proposal is an inflated gaussian on \\\mathbb{R}^p\\,
which places no mass on the simplex at all. It does not fail loudly
there: [`base::chol()`](https://rdrr.io/r/base/chol.html) accepts the
singular covariance, and the estimate of an integral that is 1 comes
back at about `2e-08`. Overriding the proposal is what makes that check
mean anything for this family.

## Arguments

- distrib:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md),
  read for its `n_dim` alone.

- theta:

  A named list of parameters. Ignored: the proposal is uniform and does
  not depend on where the family sits.

- n:

  A single positive integer, the number of draws.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `y`, an `n` by \\p\\ matrix whose rows are uniform on
the simplex, and `logd`, a numeric vector of length `n` holding the
constant \\\log\Gamma(p)\\.

## See also

[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
which consumes this,
[`distrib_rng.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.DirichletDistrib.md)
for draws from the family itself, and
[`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md)
for the generic.
