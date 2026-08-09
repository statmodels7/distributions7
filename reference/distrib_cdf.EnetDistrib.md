# Elastic-Net Distribution Function

Each half of the density is a truncated Gaussian, so \\F(q) =
\Phi(\sqrt{c}\\z - x) / (2\Phi(-x))\\ for \\z \le 0\\ and its reflection
above, with \\z = q-\mu\\ and \\x = a/\sqrt{c}\\.

## Arguments

- distrib:

  An `EnetDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu`, `lambda` and `alpha`.

- lower.tail:

  Logical; if `TRUE` (default), \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector of probabilities.

## See also

[`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
