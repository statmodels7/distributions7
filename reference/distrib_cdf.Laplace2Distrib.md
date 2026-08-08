# Laplace Distribution Function in Location and Rate

\$\$F(q; \mu, \lambda) = \begin{cases}
\dfrac{1}{2}\exp\left(\lambda(q-\mu)\right) & q \< \mu \\ 1 -
\dfrac{1}{2}\exp\left(-\lambda(q-\mu)\right) & q \ge \mu \end{cases}\$\$

## Arguments

- distrib:

  A `Laplace2Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu` and `lambda`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of cumulative probabilities.

## See also

[`laplace2_distrib`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
