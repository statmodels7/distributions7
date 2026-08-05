# Beta Analytical Observed Hessian in Its Shapes

\$\$\ell^{(\alpha\alpha)} = \psi'(\alpha+\beta) - \psi'(\alpha), \qquad
\ell^{(\alpha\beta)} = \psi'(\alpha+\beta), \qquad \ell^{(\beta\beta)} =
\psi'(\alpha+\beta) - \psi'(\beta)\$\$ free of the data.

## Arguments

- distrib:

  A `Beta2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with `alpha` and `beta`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of second derivatives.

## See also

[`beta2_distrib`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md)
