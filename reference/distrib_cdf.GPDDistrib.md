# Generalized Pareto Distribution Function

\$\$F(q) = 1 - \left(1 + \dfrac{\xi q}{\sigma}\right)^{-1/\xi}\$\$

## Arguments

- distrib:

  A `GPDDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `sigma` and `xi`.

- lower.tail:

  Logical; if `TRUE` (default), \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logarithms.

## Value

A numeric vector of cumulative probabilities.

## See also

[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
