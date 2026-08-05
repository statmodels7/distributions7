# Skew Normal Third-Order Derivatives in the Centred Parametrisation

The partition sum at order three.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object.

- y:

  The response.

- theta:

  A list with `mu`, `sigma` and `gamma1`.

- expected:

  Logical; if `TRUE`, carries the expected derivatives.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Passed to the parent.

- nsim:

  Passed to the parent.

- ...:

  Unused.

## Value

A named list of third-derivative components.

## See also

[`skewnormal2_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
