# Zero-Inflated Analytical Observed Hessian

Observed Hessian of the zero-inflated model, combining the parent's
observed Hessian with rank-one corrections at \\y=0\\ (see the old
package documentation for the full derivation).

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `zi`.

## Value

A list containing the vectors of second derivatives.

## See also

[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
