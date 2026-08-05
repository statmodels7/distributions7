# Beta-Binomial Quantile Function

The generalized inverse \\Q(p) = \min\\y : F(y) \ge p\\\\, obtained from
the exact cumulative sum.

## Arguments

- distrib:

  A `BetaBinom1Distrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing `mu` and `sigma`.

- lower.tail:

  Logical; if `TRUE` (default), `p` is \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, `p` is given as its logarithm.

## Value

A numeric vector of quantiles.

## See also

[`betabinom1_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
