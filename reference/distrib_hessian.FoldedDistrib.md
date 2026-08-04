# Folded Analytical Observed Hessian

The moment-to-cumulant relation applied to the ratios \\d^B L / L =
w\\(d^B f(x)/f(x)) + (1-w)\\(d^B f(-x)/f(-x))\\, which at second order
is the familiar mixture Hessian.

## Arguments

- distrib:

  A `FoldedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list, one component per pair of parameters.

## See also

[`folded`](https://statmodels7.github.io/distributions7/reference/folded.md)
