# Beta Analytical Gradient in Its Shapes

\$\$\dfrac{\partial\ell}{\partial\alpha} = \log y - \psi(\alpha) +
\psi(\alpha+\beta), \qquad \dfrac{\partial\ell}{\partial\beta} =
\log(1-y) - \psi(\beta) + \psi(\alpha+\beta)\$\$ The data enter only
through \\\log y\\ and \\\log(1-y)\\, which is what makes every higher
derivative free of them.

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

A named list of first derivatives.

## See also

[`beta2_distrib`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md)
