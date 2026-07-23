# Zero-Adjusted Continuous Analytical Gradient

The likelihood separates completely: for the parent's parameters the
score is the parent's score at \\y \neq 0\\ and 0 at \\y=0\\; for
\\\pi\\: \$\$\dfrac{\partial \ell}{\partial \pi} =
\mathbb{I}(y=0)\dfrac{1}{\pi} - \mathbb{I}(y \neq 0)\dfrac{1}{1-\pi}\$\$

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `za`.

## Value

A list containing the vectors of first derivatives.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
