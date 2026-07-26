# Truncated Quantile Function (Continuous)

\$\$Q_T(p) = Q\\\left(F(\ell^-;\theta) + p\\Z(\theta)\right)\$\$

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A named list of the parent's parameters.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\.

- log.p:

  Logical; if `TRUE`, probabilities are given as logs.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
