# Generalized Gamma Distribution Function

\\F(q) = P(d/p,\\ (q/a)^{p})\\, the regularized lower incomplete gamma
function, since \\(Y/a)^{p}\\ is Gamma with shape \\d/p\\ and unit rate.

## Arguments

- distrib:

  A `GenGamma1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `a`, `d` and `p`.

- lower.tail:

  Logical; if `TRUE` (default), \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logarithms.

## Value

A numeric vector of cumulative probabilities.

## See also

[`gengamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
