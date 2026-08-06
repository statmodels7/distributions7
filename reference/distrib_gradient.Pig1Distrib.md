# Poisson-Inverse Gaussian Analytical Gradient

The exact score in \\(\mu, \sigma)\\, from the compiled fourth-order
kernel described in
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

A named list with the `mu` and `sigma` components.

## See also

[`pig1_distrib`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
