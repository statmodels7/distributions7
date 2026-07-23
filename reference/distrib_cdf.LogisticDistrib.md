# Logistic Cumulative Distribution Function

Computes the cumulative distribution function for the Logistic
distribution: \$\$F(q; \mu, \sigma) = \dfrac{1}{1 +
\exp\left(-\dfrac{q-\mu}{\sigma}\right)}\$\$

## Arguments

- distrib:

  A `LogisticDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu` and `sigma`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## See also

[`logistic_distrib`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
