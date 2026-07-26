# Truncated Probability Mass Function

\$\$P_T(Y = y) = \dfrac{f(y;\theta)}{Z(\theta)}\\ \\ \\ (\ell \le y \le
u), \qquad 0 \text{ otherwise}\$\$ with \\Z(\theta) = F(u;\theta) -
F(\ell;\theta) + f(\ell;\theta)\\, the mass at the lower endpoint being
added back because the endpoint is included.

## Arguments

- distrib:

  A `TruncatedDiscreteDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

- log:

  Logical; if `TRUE`, returns the log-probability.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
