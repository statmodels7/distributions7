# Log-CDF Gradient of a Reparametrized Distribution

The chain rule on the parent's cdf derivatives, which is exact whenever
the parent's are; when the parent differences its own cdf, so does this,
the chain having nothing closed to carry.

## Arguments

- distrib:

  A reparametrized distribution.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of the new parameters.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
