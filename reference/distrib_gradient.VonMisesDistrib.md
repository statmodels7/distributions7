# von Mises Analytical Gradient

\$\$\dfrac{\partial \ell}{\partial \mu} = \kappa \sin(y - \mu), \qquad
\dfrac{\partial \ell}{\partial \kappa} = \cos(y - \mu) - A(\kappa)\$\$
with \\A = I_1/I_0\\, the derivative of \\\log I_0\\.

## Arguments

- distrib:

  A `VonMisesDistrib` object.

- y:

  A numeric vector of angles.

- theta:

  A list containing `mu` and `kappa`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list with the `mu` and `kappa` components.

## See also

[`vonmises_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises_distrib.md)
