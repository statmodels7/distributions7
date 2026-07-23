# Zero-Adjusted Continuous Analytical Expected Hessian

Expected Hessian: \\E\[H\_{\pi\pi}\] = -\dfrac{1}{\pi(1-\pi)}\\, mixed
blocks are 0, and \\E\[H\_{\theta\theta}\] = (1-\pi) E\[H_W\]\\.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `za`.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
