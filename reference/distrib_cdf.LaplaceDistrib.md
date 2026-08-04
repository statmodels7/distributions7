# Laplace Cumulative Distribution Function

Computes the cumulative distribution function for the Laplace
distribution: \$\$F(q; \mu, b) = \begin{cases}
\dfrac{1}{2}\exp\left(\dfrac{q-\mu}{b}\right) & q \< \mu \\ 1 -
\dfrac{1}{2}\exp\left(-\dfrac{q-\mu}{b}\right) & q \ge \mu
\end{cases}\$\$

## Arguments

- distrib:

  A `LaplaceDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu` and `b`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of cumulative probabilities.

## See also

[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
