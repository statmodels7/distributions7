# Zero-Adjusted Continuous Quantile Function

Inverts the mixed CDF, handling the jump of size \\\pi\\ at 0: quantiles
falling in the jump are 0; on either side the parent's quantile function
is used on the rescaled probability.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list with the parent's parameters followed by `za`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\.

- log.p:

  Logical; if `TRUE`, probabilities are given as logs.

## Value

A numeric vector of quantiles.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
