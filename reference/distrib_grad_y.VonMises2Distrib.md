# von Mises Response and Mixed Derivatives in the Resultant Length

Obtains the derivatives from those of the concentration parametrization,
at \\\kappa = A^{-1}(\rho)\\. The response derivatives are the parent's
unchanged; the mixed block is the parent's carried by \\\kappa'(\rho)\\,
and the location passes through.

## Arguments

- distrib:

  A `VonMises2Distrib` object.

- y:

  A numeric vector of angles.

- theta:

  A list containing `mu` and `rho`.

- ...:

  Unused.

- scale:

  Handled by the generic after dispatch; this method returns the
  parameter scale.

## Value

A numeric vector.
