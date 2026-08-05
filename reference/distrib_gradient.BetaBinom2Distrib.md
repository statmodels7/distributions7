# Beta-Binomial Analytical Gradient in Its Shapes

\$\$\dfrac{\partial\ell}{\partial\alpha} = \psi(y+\alpha) -
\psi(n+\alpha+\beta) - \psi(\alpha) + \psi(\alpha+\beta)\$\$ and the
same with \\n-y\\ and \\\beta\\.

## Arguments

- distrib:

  A `BetaBinom2Distrib` object.

- y:

  A numeric vector of counts.

- theta:

  A list with `alpha` and `beta`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of first derivatives.

## See also

[`betabinom2_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
