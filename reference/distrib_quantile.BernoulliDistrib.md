# Bernoulli Quantile Function

Computes the quantile function for the Bernoulli distribution, the
generalized inverse of the CDF: \$\$Q(p; \mu) = \begin{cases} 0 & p \le
1-\mu \\ 1 & p \> 1-\mu \end{cases}\$\$

## Arguments

- distrib:

  A `BernoulliDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameter `mu`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of quantiles.

## See also

[`bernoulli_distrib`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md)
