# Multinomial Probability Mass Function

\$\$P(Y = y) = \dfrac{n!}{\prod_j y_j!}\prod_j p_j^{y_j}\$\$ for a row
of non-negative integers summing to \\n\\.

## Arguments

- distrib:

  A `MultinomialDistrib` object.

- y:

  A matrix with one row per observation, each summing to the size.

- theta:

  A named list of parameters.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector, one entry per row of `y`.

## See also

[`multinomial_distrib`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
