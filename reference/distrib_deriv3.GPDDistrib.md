# Generalized Pareto Third and Fourth Derivatives

Closed form at both orders, from
[`gpd_components`](https://statmodels7.github.io/distributions7/reference/gpd_components.md):
the log-density splits into \\-\log\sigma\\, \\-\log t\\ and
\\-\log(t)/\xi\\, and the last is taken from its series where the
Leibniz form's terms cancel.

## Arguments

- distrib:

  A `GPDDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `sigma` and `xi`.

- expected:

  Logical; if `TRUE`, the expected derivatives.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  The approximation used when `expected` is `TRUE`.

- nsim:

  Monte Carlo draws when `approx = "mc"`.

- ...:

  Unused.

## Value

A named list of third-derivative components.

A named list of fourth-derivative components.

## See also

[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
