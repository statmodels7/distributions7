# Folded Response Hessian

\$\$\dfrac{\partial^2 \ell}{\partial x^2} = w\left(h(x) +
g(x)^2\right) + (1-w)\left(h(-x) + g(-x)^2\right) - \left(w g(x) - (1-w)
g(-x)\right)^2\$\$ with \\g\\ and \\h\\ the parent's first and second
response derivatives.

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
