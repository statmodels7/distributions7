# Beta Quantile Function

Computes the quantile function for the Beta distribution as the inverse
of the CDF, \\Q(p; \mu, \phi) = F^{-1}(p; \mu, \phi)\\. There is no
elementary closed form; it is obtained numerically (via
[`qbeta`](https://rdrr.io/r/stats/Beta.html)) on the mean/precision
reparameterization \\\alpha = \mu\phi\\, \\\beta = (1-\mu)\phi\\.

## Arguments

- distrib:

  A `BetaDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameters `mu` and `phi`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\,
  otherwise \\P(Y \> p)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## See also

[`beta_distrib`](https://statmodels7.github.io/distributions7/reference/beta_distrib.md)
