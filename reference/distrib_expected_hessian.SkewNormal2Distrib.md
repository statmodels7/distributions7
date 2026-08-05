# Skew Normal Expected Hessian in the Centered Parametrization

The parent's expected information carried by the same congruence. It is
**non-singular at zero skewness**, which the direct parametrization's is
not: there the score for \\\alpha\\ is proportional to the score for the
location and the information loses a rank.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object.

- y:

  The response.

- theta:

  A list with `mu`, `sigma` and `gamma1`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Passed to the parent.

- nsim:

  Passed to the parent.

- ...:

  Unused.

## Value

A named list of expected second derivatives.

## See also

[`skewnormal2_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
