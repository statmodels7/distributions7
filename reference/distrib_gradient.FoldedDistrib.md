# Folded Analytical Gradient

\$\$\dfrac{\partial \ell}{\partial \theta_i} = w\\ s_i(x) + (1-w)\\
s_i(-x), \qquad w = \dfrac{f(x)}{L(x)}\$\$ the score of a two-component
mixture, the components being the two preimages.

## Arguments

- distrib:

  A `FoldedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list, one component per parameter.

## See also

[`folded`](https://statmodels7.github.io/distributions7/reference/folded.md)
