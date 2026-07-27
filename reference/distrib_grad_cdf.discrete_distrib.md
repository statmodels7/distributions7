# Log-CDF Gradient for Discrete Distributions

Exact: the partial expectation of the score is a finite sum over the
support up to \\q\\, so nothing is differenced.

## Arguments

- distrib:

  A `discrete_distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.
