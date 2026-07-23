# Beta Analytical Third-Order Derivatives

Closed-form third-order derivatives of the Beta log-density. Because the
\\y\\-terms of the log-density are linear in the parameters, these
derivatives do not depend on \\y\\; the observed and expected values
therefore coincide.

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

A named list of third-derivative component vectors.

## See also

[`beta_distrib`](https://statmodels7.github.io/distributions7/reference/beta_distrib.md)
