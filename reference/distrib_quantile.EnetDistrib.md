# Elastic-Net Quantile Function

The distribution function of
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.EnetDistrib.md)
inverted in closed form, each half being a truncated Gaussian.

## Arguments

- distrib:

  An `EnetDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing `mu`, `lambda` and `alpha`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, `p` is on the log scale.

## Value

A numeric vector of quantiles.

## See also

[`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
