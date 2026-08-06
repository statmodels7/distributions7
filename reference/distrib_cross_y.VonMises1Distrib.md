# von Mises Mixed Derivatives

Closed form. With \\\ell^{(y)} = -\kappa\sin(y-\mu)\\, the location
component is \\\kappa\cos(y-\mu)\\ and the concentration one
\\-\sin(y-\mu)\\.

## Arguments

- distrib:

  A `VonMises1Distrib` object.

- y:

  A numeric vector of angles.

- theta:

  A list containing `mu` and `kappa`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
