# Generalised Pareto Quantile Function

\$\$Q(p) = \dfrac{\sigma}{\xi}\left((1-p)^{-\xi} - 1\right)\$\$ with
\\-\sigma\log(1-p)\\ at \\\xi = 0\\.

## Arguments

- distrib:

  A `GPDDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing `sigma` and `xi`.

- lower.tail:

  Logical; if `TRUE` (default), `p` is \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, `p` is given as its logarithm.

## Value

A numeric vector of quantiles.

## See also

[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
