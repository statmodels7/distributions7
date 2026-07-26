# Truncated Analytical Observed Hessian (Continuous)

\$\$\dfrac{\partial^2 \ell_T}{\partial\theta_i\partial\theta_j} =
H\_{ij}(y) - M\_{ij} + m_i m_j\$\$ with \\M\_{ij} =
\mathbb{E}\_T\[H\_{ij} + s_i s_j\]\\.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

## Value

A list containing the vectors of second derivatives.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
