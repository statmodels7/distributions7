# Student's t Quantile Function

Computes the quantile function (inverse CDF) for the (location-scale)
Student's t distribution: \$\$Q(p; \mu, \sigma, \nu) = \mu + \sigma\\
T\_\nu^{-1}(p)\$\$ where \\T\_\nu^{-1}\\ is the standard Student's t
quantile function with \\\nu\\ degrees of freedom.

## Arguments

- distrib:

  A `StudentT1Distrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameters `mu`, `sigma`, and `nu`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of quantiles.

## See also

[`student_t1_distrib`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
