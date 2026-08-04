# Folded Density

\$\$L(x; \theta) = f(x; \theta) + f(-x; \theta), \qquad x \ge 0\$\$ the
two preimages of \\x\\ under the absolute value added together.

## Arguments

- distrib:

  A `FoldedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`folded`](https://statmodels7.github.io/distributions7/reference/folded.md)
