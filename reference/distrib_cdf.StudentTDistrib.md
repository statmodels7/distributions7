# Student's t Cumulative Distribution Function

Computes the cumulative distribution function for the (location-scale)
Student's t distribution: \$\$F(q; \mu, \sigma, \nu) =
T\_\nu\\\left(\dfrac{q-\mu}{\sigma}\right)\$\$ where \\T\_\nu\\ is the
CDF of the standard Student's t with \\\nu\\ degrees of freedom.

## Arguments

- distrib:

  A `StudentTDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu`, `sigma`, and `nu`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## See also

[`student_t_distrib`](https://statmodels7.github.io/distributions7/reference/student_t_distrib.md)
