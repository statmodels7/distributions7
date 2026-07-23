# Default Numerical Quantile Function for Discrete Distributions

Fallback method: discrete distributions that do not implement an
analytical quantile function get one by inverting the cumulative pmf
table: the quantile is the smallest support point \\k\\ with \\F(k) \ge
p\\.

## Arguments

- distrib:

  An object inheriting from class `"discrete_distrib"`.

- p:

  A numeric vector of probabilities.

- theta:

  A named list of parameters.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\.

- log.p:

  Logical; if `TRUE`, probabilities are given as logs.
