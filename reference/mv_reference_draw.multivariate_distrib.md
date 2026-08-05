# An Inflated Gaussian Proposal

The default: a gaussian with the family's mean and twice its covariance.
Requires the covariance to be non-singular, which is what a support of
full dimension gives.

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters.

- n:

  The number of draws.

- ...:

  Unused.

## Value

A list with the draws `y` and their log-density `logd`.

## See also

[`mv_reference_draw`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md)
