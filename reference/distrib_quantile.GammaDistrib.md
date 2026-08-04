# Gamma Quantile Function

Computes the quantile function for the Gamma distribution as the inverse
of the CDF, \\Q(p; \mu, \sigma^2) = F^{-1}(p; \mu, \sigma^2)\\. There is
no elementary closed form; it is obtained numerically (via
[`qgamma`](https://rdrr.io/r/stats/GammaDist.html)) on the shape/rate
reparameterization \\\alpha = \mu^2/\sigma^2\\, \\\lambda =
\mu/\sigma^2\\.

## Arguments

- distrib:

  A `GammaDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameters `mu` and `sigma2`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of quantiles.

## See also

[`gamma_distrib`](https://statmodels7.github.io/distributions7/reference/gamma_distrib.md)
