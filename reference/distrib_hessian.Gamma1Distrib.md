# Gamma Analytical Observed Hessian in Mean and Dispersion

Closed form. The derivatives in \\\phi\\ are those in \\s = 1/\phi\\
carried across by the one-variable chain rule, which is what keeps the
polygamma functions to one evaluation each.

## Arguments

- distrib:

  A `Gamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu` and `phi`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of second derivatives.

## See also

[`gamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
