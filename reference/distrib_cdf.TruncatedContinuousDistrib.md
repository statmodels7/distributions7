# Truncated Cumulative Distribution Function (Continuous)

\$\$F_T(q) = \dfrac{F(q;\theta) - F(\ell^-;\theta)}{Z(\theta)}\$\$
clamped to \\\[0,1\]\\.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of the parent's parameters.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logs.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
