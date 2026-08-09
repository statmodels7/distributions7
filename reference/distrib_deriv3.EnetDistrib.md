# Elastic-Net Third and Fourth Derivatives

The third and fourth derivatives of the log-density, assembled in the
two rates and carried onto \\(\lambda, \alpha)\\ by the bilinear map.

## Arguments

- distrib:

  An `EnetDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `lambda` and `alpha`.

- expected:

  Logical; whether to return the expected derivatives.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  How the expectation is approximated.

- nsim:

  Monte Carlo sample size, used when `approx = "mc"`.

- ...:

  Unused.

## Value

A named list of third-derivative components.

A named list of fourth-derivative components.

## Details

Written in \\(\mu, a, c)\\ the log-density is quadratic in \\z\\ and
linear in each rate, so at these orders only the normalizing constant
contributes, apart from \\\partial^{3}\ell/\partial\mu^{2}\partial c =
-1\\; see `.enet_ac_derivs`. Running the same assembly at orders one and
two reproduces the hand-written score and Hessian, which is what
licenses it at the orders where there is nothing to compare against.

The expected derivatives have no closed form here and go through
[`expected_derivative`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md).

## See also

[`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
