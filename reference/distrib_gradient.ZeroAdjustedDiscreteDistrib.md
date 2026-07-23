# Zero-Adjusted Discrete Analytical Gradient

The hurdle likelihood separates: for the parent's parameters the score
at \\y\>0\\ is the parent's score plus the truncation correction
\\\dfrac{f(0)}{1-f(0)} S(0)\\ (and 0 at \\y=0\\); for \\\pi\\:
\$\$\dfrac{\partial \ell}{\partial \pi} =
\mathbb{I}(y=0)\dfrac{1}{\pi} - \mathbb{I}(y\>0)\dfrac{1}{1-\pi}\$\$

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `za`.

## Value

A list containing the vectors of first derivatives.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
