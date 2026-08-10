# Generalized Pareto Log-CDF Derivatives

Closed form at every order, from the survival function \\S = (1 + \xi
q/\sigma)^{-1/\xi}\\. Its logarithm is written \\-(q/\sigma)\Lambda(\xi
q/\sigma)\\ with \\\Lambda(u) = \log(1+u)/u\\, which carries no division
by the shape, so the exponential limit \\\xi \to 0\\ is an ordinary
point of the formula rather than a branch.

## Arguments

- distrib:

  A `GPDDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `sigma` and `xi`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

- ...:

  Unused.

## Value

A named list, one vector per component.

## See also

[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
