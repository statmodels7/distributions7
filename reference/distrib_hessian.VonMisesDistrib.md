# von Mises Analytical Observed Hessian

\$\$\ell^{(\mu\mu)} = -\kappa\cos(y-\mu), \qquad \ell^{(\mu\kappa)} =
\sin(y-\mu), \qquad \ell^{(\kappa\kappa)} = -A'(\kappa)\$\$ the last
free of the data, \\\log I_0\\ being the only place \\\kappa\\ appears
alone.

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

A named list of second-derivative components.

## See also

[`vonmises_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises_distrib.md)
