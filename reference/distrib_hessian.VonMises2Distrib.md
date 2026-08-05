# von Mises Analytical Observed Hessian in the Resultant Length

The concentration parametrisation's second derivatives carried through
the one-variable chain rule, \$\$\ell^{(\rho\rho)} =
\ell^{(\kappa\kappa)}(\kappa')^2 + \ell^{(\kappa)}\kappa'',\$\$ with
\\\ell^{(\kappa\kappa)} = -A'(\kappa)\\ and \\\ell^{(\mu\kappa)} =
\sin(y-\mu)\\.

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

A named list of second derivatives.

## See also

[`vonmises2_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md)
