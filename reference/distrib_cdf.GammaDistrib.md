# Gamma Cumulative Distribution Function

Computes the cumulative distribution function for the Gamma
distribution, using the shape/rate reparameterization \\\alpha =
\mu^2/\sigma^2\\, \\\lambda = \mu/\sigma^2\\: \$\$F(q; \mu, \sigma^2) =
\dfrac{\gamma(\alpha, \lambda q)}{\Gamma(\alpha)}\$\$ where
\\\gamma(\cdot, \cdot)\\ is the lower incomplete gamma function.

## Arguments

- distrib:

  A `GammaDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu` and `sigma2`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## See also

[`gamma_distrib`](https://statmodels7.github.io/distributions7/reference/gamma_distrib.md)
