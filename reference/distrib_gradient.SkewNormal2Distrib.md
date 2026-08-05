# Skew Normal Gradient in the Centred Parametrisation

The parent's score carried by the Jacobian of the centred-to-direct map.
The components in \\\gamma_1\\ stay of order one however small
\\\gamma_1\\ is, although the Jacobian itself grows without bound: the
divergent parts cancel, which is what the centred parametrisation is
for.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object.

- y:

  The response.

- theta:

  A list with `mu`, `sigma` and `gamma1`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of first derivatives.

## See also

[`skewnormal2_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
