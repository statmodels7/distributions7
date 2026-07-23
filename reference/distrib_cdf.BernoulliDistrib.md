# Bernoulli Cumulative Distribution Function

Computes the cumulative distribution function for the Bernoulli
distribution: \$\$F(q; \mu) = \begin{cases} 0 & q \< 0 \\ 1-\mu & 0 \le
q \< 1 \\ 1 & q \ge 1 \end{cases}\$\$

## Arguments

- distrib:

  A `BernoulliDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameter `mu`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## See also

[`bernoulli_distrib`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md)
