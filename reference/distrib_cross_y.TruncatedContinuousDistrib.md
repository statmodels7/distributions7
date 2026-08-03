# Mixed Derivatives of a Truncated Distribution

The parent's mixed derivatives, unchanged: the truncated log-density is
the parent's minus \\\log Z(\theta)\\, and the normalising constant does
not depend on \\y\\, so it vanishes from any derivative that involves
the response.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
