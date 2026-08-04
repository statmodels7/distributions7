# Zero-Adjusted Continuous Cumulative Distribution Function

\$\$F_Y(q) = (1-\pi)F_W(q;\theta) + \pi\\\mathbb{I}(q \ge 0)\$\$

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list with the parent's parameters followed by `za`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logs.

## Value

A numeric vector of cumulative probabilities.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
