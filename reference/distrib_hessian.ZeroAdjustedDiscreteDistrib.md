# Zero-Adjusted Discrete Analytical Observed Hessian

Observed Hessian of the hurdle model. The mixed blocks are identically
zero; the truncation adds the correction \\H\_{corr} =
\dfrac{(1-f(0))f''(0) + f'(0)^2}{(1-f(0))^2}\\ for \\y \> 0\\.

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `za`.

## Value

A list containing the vectors of second derivatives.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
