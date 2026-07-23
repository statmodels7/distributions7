# Zero-Adjusted Discrete Analytical Expected Hessian

Expected Hessian of the hurdle model: \\E\[H\_{\pi\pi}\] =
-\dfrac{1}{\pi(1-\pi)}\\, mixed blocks are 0, and
\\E\[H\_{\theta\theta}\] = (1-\pi)\left(\dfrac{E\[H\] -
f(0)H(0)}{1-f(0)} + H\_{corr}\right)\\.

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `za`.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
