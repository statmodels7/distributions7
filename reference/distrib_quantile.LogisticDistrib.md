# Logistic Quantile Function

Computes the quantile function (inverse CDF) for the Logistic
distribution: \$\$Q(p; \mu, \sigma) = \mu + \sigma
\log\left(\dfrac{p}{1-p}\right)\$\$

## Arguments

- distrib:

  A `LogisticDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameters `mu` and `sigma`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## See also

[`logistic_distrib`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
