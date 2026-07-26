# Truncated Probability Density Function

\$\$f_T(y) = \dfrac{f(y;\theta)}{Z(\theta)}\\ \\ \\ (\ell \le y \le u),
\qquad 0 \text{ otherwise}\$\$ with \\Z(\theta) = F(u;\theta) -
F(\ell;\theta)\\.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

- log:

  Logical; if `TRUE`, returns the log-density.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
