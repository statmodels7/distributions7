# Zero-Adjusted Discrete Cumulative Distribution Function

\$\$F\_{ZA}(q) = \pi + (1-\pi)\dfrac{F(q;\theta) - f(0;\theta)}{1 -
f(0;\theta)} \quad (q \ge 0)\$\$

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list with the parent's parameters followed by `za`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logs.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
