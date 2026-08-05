# von Mises Analytical Gradient in the Resultant Length

\$\$\dfrac{\partial\ell}{\partial\mu} = \kappa\sin(y-\mu), \qquad
\dfrac{\partial\ell}{\partial\rho} = \left\\\cos(y-\mu) -
A(\kappa)\right\\\kappa'(\rho)\$\$ The map touches only the second
parameter, so the chain rule is the one-variable one and no cancellation
is involved.

## Arguments

- distrib:

  A `VonMises2Distrib` object.

- y:

  A numeric vector of angles.

- theta:

  A list with `mu` and `rho`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of first derivatives.

## See also

[`vonmises2_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md)
