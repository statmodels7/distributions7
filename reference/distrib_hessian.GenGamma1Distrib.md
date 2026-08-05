# Generalised Gamma Analytical Observed Hessian

The second derivatives of the same expressions. The mixed \\a\\–\\d\\
component is \\-1/a\\, free of the data, the scale and the first shape
entering the log-density through \\-d\log a\\ alone.

## Arguments

- distrib:

  A `GenGamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `a`, `d` and `p`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of second-derivative components.

## See also

[`gengamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
