# Beta-Binomial Cumulative Distribution Function

The cumulative sum of the mass function over the finite support, which
is exact rather than approximated.

## Arguments

- distrib:

  A `BetaBinom1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `sigma`.

- lower.tail:

  Logical; if `TRUE` (default), \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logarithms.

## Value

A numeric vector of cumulative probabilities.

## See also

[`betabinom1_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
