# Zero-Inflated Cumulative Distribution Function

\$\$F\_{ZI}(q) = (1-\zeta) F(q; \theta) + \zeta\\\mathbb{I}(q \ge 0)\$\$

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list with the parent's parameters followed by `zi`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logs.

## See also

[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
