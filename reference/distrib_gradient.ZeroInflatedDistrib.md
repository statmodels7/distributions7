# Zero-Inflated Analytical Gradient

Score function of the zero-inflated model. For the parent's parameters
the score is the parent's score weighted by \\w = (1-\zeta)f(0)/L_0\\ at
\\y=0\\ (and 1 otherwise), where \\L_0 = \zeta + (1-\zeta)f(0)\\. For
\\\zeta\\: \$\$\dfrac{\partial \ell}{\partial \zeta} =
\mathbb{I}(y=0)\dfrac{1-f(0)}{L_0} -
\mathbb{I}(y\>0)\dfrac{1}{1-\zeta}\$\$

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `zi`.

## Value

A list containing the vectors of first derivatives.

## See also

[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
