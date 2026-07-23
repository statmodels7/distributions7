# Beta Analytical Fourth-Order Derivatives

Closed-form fourth-order derivatives of the Beta log-density. As for the
third order, they do not depend on \\y\\, so observed and expected
coincide.

## Arguments

- distrib:

  A `BetaDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `phi`.

- expected:

  Logical; ignored (observed and expected coincide).

## Value

A named list of fourth-derivative component vectors.

## See also

[`beta_distrib`](https://statmodels7.github.io/distributions7/reference/beta_distrib.md)
