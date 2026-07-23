# Zero-Inflated Quantile Function

Inverts the mixture CDF: the quantile is 0 for \\p \le \zeta +
(1-\zeta)F(0;\theta)\\, otherwise \\Q\left(\dfrac{p-\zeta}{1-\zeta};
\theta\right)\\.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list with the parent's parameters followed by `zi`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\.

- log.p:

  Logical; if `TRUE`, probabilities are given as logs.

## See also

[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
