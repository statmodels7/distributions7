# Folded Response Gradient

\$\$\dfrac{\partial \ell}{\partial x} = w\\ g(x) - (1-w)\\ g(-x)\$\$
with \\g = \partial \log f/\partial y\\. The minus sign is the chain
rule through the reflected preimage.

## Arguments

- distrib:

  A `FoldedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

- ...:

  Unused.

## Value

A numeric vector.

## See also

[`folded`](https://statmodels7.github.io/distributions7/reference/folded.md)
