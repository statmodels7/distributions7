# Zero-Adjusted Discrete Quantile Function

Inverts the hurdle CDF: 0 for \\p \le \pi\\, otherwise the parent
quantile at \\u(1-f(0)) + f(0)\\ with \\u = (p-\pi)/(1-\pi)\\.

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list with the parent's parameters followed by `za`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\.

- log.p:

  Logical; if `TRUE`, probabilities are given as logs.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
