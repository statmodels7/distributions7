# Poisson-Inverse Gaussian Analytical Third Derivatives

The exact third derivatives in \\(\mu, \sigma)\\, from the compiled
fourth-order kernel described in
[`pig_hd_block`](https://statmodels7.github.io/distributions7/reference/pig_hd_block.md).

## Arguments

- distrib:

  A `Pig1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma`.

- expected:

  Logical; the expected version goes through the generic's strategies.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx, nsim:

  Passed on when `expected` is `TRUE`.

- ...:

  Unused.

## Value

A named list of third-derivative components.

## See also

[`pig1_distrib`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
