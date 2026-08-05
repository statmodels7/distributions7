# Beta Cumulative Distribution Function

Computes the cumulative distribution function for the Beta distribution,
using the mean/precision reparameterization \\\alpha = \mu\phi\\,
\\\beta = (1-\mu)\phi\\: \$\$F(q; \mu, \phi) = I_q(\alpha, \beta)\$\$
where \\I_q(\cdot, \cdot)\\ is the regularized incomplete beta function.

## Arguments

- distrib:

  A `Beta1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu` and `phi`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of cumulative probabilities.

## See also

[`beta1_distrib`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
