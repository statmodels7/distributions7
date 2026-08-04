# Default Numerical CDF for Discrete Distributions

Fallback method: discrete distributions that do not implement an
analytical CDF get one by summing the pmf from the (finite) lower bound
of the support up to \\\lfloor q \rfloor\\.

## Arguments

- distrib:

  An object inheriting from class `"discrete_distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logs.

## Value

A numeric vector of cumulative probabilities.
