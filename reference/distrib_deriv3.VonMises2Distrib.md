# von Mises Third and Fourth Derivatives in the Resultant Length

Closed form at both orders. Every derivative that carries at least one
\\\mu\\ is \\D_a \kappa^{(b)}\\, where \\D_a\\ is the \\a\\-th
\\\mu\\-derivative of \\\cos(y-\mu)\\ and \\\kappa^{(b)}\\ the \\b\\-th
derivative of \\A^{-1}\\: the concentration parametrization's
\\\mu\\-derivatives are linear in \\\kappa\\, so the composition
collapses to a single term. The pure \\\rho\\ components carry the full
one-variable chain rule on \\\log I_0\\.

## Arguments

- distrib:

  A `VonMises2Distrib` object.

- y:

  A numeric vector of angles.

- theta:

  A list containing `mu` and `rho`.

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

[`vonmises2_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md)
