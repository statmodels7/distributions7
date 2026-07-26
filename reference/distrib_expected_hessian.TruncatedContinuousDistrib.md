# Truncated Analytical Expected Hessian (Continuous)

\$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell_T}{\partial\theta_i\partial\theta_j}\right\] =
-\operatorname{Cov}\_T(s_i, s_j)\$\$ the covariance of the parent's
score under the truncated distribution.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
