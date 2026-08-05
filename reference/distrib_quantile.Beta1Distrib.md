# Beta Quantile Function

Computes the quantile function for the Beta distribution as the inverse
of the CDF, \\Q(p; \mu, \phi) = F^{-1}(p; \mu, \phi)\\. There is no
elementary closed form; it is obtained numerically (via
[`qbeta`](https://rdrr.io/r/stats/Beta.html)) on the mean/precision
reparameterization \\\alpha = \mu\phi\\, \\\beta = (1-\mu)\phi\\.

## Arguments

- distrib:

  A `Beta1Distrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameters `mu` and `phi`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of quantiles.

## See also

[`beta1_distrib`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
