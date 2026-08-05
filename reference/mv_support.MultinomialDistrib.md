# The Support Points of a Multinomial

Every vector of non-negative integers summing to the size, enumerated as
a matrix with one row per point.

## Arguments

- distrib:

  A `MultinomialDistrib` object.

- theta:

  Ignored; the support does not depend on the parameters.

- ...:

  Unused.

## Value

A matrix with `choose(n + p - 1, p - 1)` rows.

## See also

[`multinomial_distrib`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
