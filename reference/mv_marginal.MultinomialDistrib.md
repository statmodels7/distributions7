# Multinomial Marginal

A single coordinate is \\\mathrm{Binomial}(n, p_j)\\, the other outcomes
collapsing into a single failure.

## Arguments

- distrib:

  A `MultinomialDistrib` object.

- theta:

  A named list of parameters.

- which:

  The coordinate wanted.

- ...:

  Unused.

## Value

A list with the marginal `distrib` and its `theta`.

## See also

[`multinomial_distrib`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
