# Generalised Gamma Quantile Function

\\Q(u) = a\\\\Q\_{\Gamma}(u; d/p)\\^{1/p}\\, inverting the same
representation the distribution function uses.

## Arguments

- distrib:

  A `GenGammaDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing `a`, `d` and `p`.

- lower.tail:

  Logical; if `TRUE` (default), `p` is \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, `p` is given as its logarithm.

## Value

A numeric vector of quantiles.

## See also

[`gengamma_distrib`](https://statmodels7.github.io/distributions7/reference/gengamma_distrib.md)
