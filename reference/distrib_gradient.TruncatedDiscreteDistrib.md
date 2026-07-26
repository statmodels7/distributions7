# Truncated Analytical Gradient (Discrete)

\$\$\dfrac{\partial \ell_T}{\partial\theta_i} = s_i(y) - m_i, \qquad m_i
= \mathbb{E}\_T\[s_i\]\$\$

## Arguments

- distrib:

  A `TruncatedDiscreteDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

## Value

A list containing the vectors of first derivatives.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
