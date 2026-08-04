# Geometric Cumulative Distribution Function

\$\$F(q; \mu) = 1 - \left(\dfrac{\mu}{1+\mu}\right)^{\lfloor q \rfloor +
1}\$\$

## Arguments

- distrib:

  A `GeometricDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameter `mu`.

- lower.tail:

  Logical; if `TRUE` (default), \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logarithms.

## Value

A numeric vector of cumulative probabilities.

## See also

[`geometric_distrib`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
