# Inverse-Gaussian Quantile Function

Computes the quantile function for the Inverse-Gaussian distribution as
the inverse of the CDF, \\Q(p; \mu, \phi) = F^{-1}(p; \mu, \phi)\\.
There is no closed form; it is obtained numerically via
[`qinvgauss`](https://rdrr.io/pkg/statmod/man/invgauss.html).

## Arguments

- distrib:

  An `InvGauss1Distrib` object.

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

[`invgauss1_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
