# Logistic Log-CDF Second Derivatives

Closed form, from the location-scale structure.

## Arguments

- distrib:

  A `LogisticDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `sigma`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## See also

[`logistic_distrib`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
