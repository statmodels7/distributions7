# Beta-Binomial Fourth-Order Derivatives in Its Shapes

Closed form, with the expectation an exact sum over the support.

## Arguments

- distrib:

  A `BetaBinom2Distrib` object.

- y:

  A numeric vector of counts.

- theta:

  A list with `alpha` and `beta`.

- expected:

  Logical; if `TRUE`, returns the expected derivatives.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of fourth-derivative components.

## See also

[`betabinom2_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
