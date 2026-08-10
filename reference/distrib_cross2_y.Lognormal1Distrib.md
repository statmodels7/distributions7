# Lognormal Second-Response and Second-Order Mixed Derivatives

The gaussian's own components at \\t = \log y\\: the transformation
carries no parameter, so \\t\\ does not move with \\\theta\\ and only
the Jacobian of the response derivatives enters.

## Arguments

- distrib:

  A `Lognormal1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma2`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list, keyed by parameter or by parameter pair.
