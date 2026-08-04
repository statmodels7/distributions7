# Folded Analytical Fourth-Order Derivatives

Fourth-order derivatives from the same partition sums as the Hessian;
see
[`distrib_hessian.FoldedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.FoldedDistrib.md).

## Arguments

- distrib:

  A `FoldedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

- expected:

  Logical; if `TRUE`, the expectation is approximated.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  How the expectation is approximated when requested.

- nsim:

  Monte Carlo sample size.

- ...:

  Unused.

## Value

A named list of fourth-order components.

## See also

[`folded`](https://statmodels7.github.io/distributions7/reference/folded.md)
