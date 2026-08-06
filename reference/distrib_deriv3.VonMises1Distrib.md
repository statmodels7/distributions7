# von Mises Third and Fourth Derivatives

Closed form at both orders. The log-density is \\\kappa\cos(y-\mu) -
\log(2\pi I_0(\kappa))\\, linear in \\\kappa\\ apart from the
normalizing constant, so every component with one \\\mu\\ and two or
more \\\kappa\\ vanishes exactly. The pure \\\mu\\ components cycle
through \\\kappa\\\sin, -\cos, -\sin, \cos\\(y-\mu)\\, and the pure
\\\kappa\\ ones are minus the derivatives of \\A(\kappa)\\, which
numericals7 supplies from the Riccati recursion \\A' = 1 - A/\kappa -
A^2\\.

## Arguments

- distrib:

  A `VonMises1Distrib` object.

- y:

  A numeric vector of angles.

- theta:

  A list containing `mu` and `kappa`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- expected:

  Logical; if `TRUE`, the expected derivatives.

- approx:

  The approximation used when `expected` is `TRUE`.

- nsim:

  Monte Carlo draws when `approx = "mc"`.

- ...:

  Unused.

## Value

A named list of third-derivative components.

A named list of fourth-derivative components.

## See also

[`vonmises1_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
