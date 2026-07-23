# Binomial Cumulative Distribution Function

Computes the cumulative distribution function for the Binomial
distribution: \$\$F(q; \mu, n) = \sum\_{k=0}^{\lfloor q \rfloor}
\dbinom{n}{k} \mu^k (1-\mu)^{n-k}\$\$

## Arguments

- distrib:

  A `BinomialDistrib` object.

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

[`binomial_distrib`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
