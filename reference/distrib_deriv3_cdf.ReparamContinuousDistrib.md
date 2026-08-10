# Third and Fourth Log-CDF Derivatives of a Reparametrized Distribution

The chain rule on the parent's, exact whenever the parent's are and the
stencil otherwise, as at the two orders below.

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

- ...:

  Unused.

## Value

A named list, one vector per component.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
