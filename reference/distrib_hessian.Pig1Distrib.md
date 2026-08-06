# Poisson-Inverse Gaussian Analytical Observed Hessian

The exact second derivatives in \\(\mu, \sigma)\\, from the compiled
fourth-order kernel described in
[`pig_hd_block`](https://statmodels7.github.io/distributions7/reference/pig_hd_block.md).

## Arguments

- distrib:

  A `Pig1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of second-derivative components.

## See also

[`pig1_distrib`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
