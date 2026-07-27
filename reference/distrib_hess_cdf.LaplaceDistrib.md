# Laplace Log-CDF Second Derivatives

Closed form, exact away from \\y = \mu\\. At the kink the second
derivative of \\F\\ in \\\mu\\ genuinely does not exist — it jumps
between \\\pm 1/(2b^{2})\\ — so the value returned there is the
one-sided limit the sign convention picks out, and is reported rather
than smoothed.

## Arguments

- distrib:

  A `LaplaceDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `b`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## See also

[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
