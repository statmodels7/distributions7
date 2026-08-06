# Orthogonal Poisson-Inverse Gaussian Analytical Fourth Derivatives

The exact fourth derivatives in \\(\mu, \alpha)\\, from the compiled
fourth-order kernel described in
[`pig_hd_block`](https://statmodels7.github.io/distributions7/reference/pig_hd_block.md).

## Arguments

- distrib:

  A `Pig2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `alpha`.

- expected:

  Logical; the expected version goes through the generic's strategies.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx, nsim:

  Passed on when `expected` is `TRUE`.

- ...:

  Unused.

## Value

A named list of fourth-derivative components.

## See also

[`pig2_distrib`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
