# Truncated Third Derivatives (Continuous)

\\\ell_T = \ell - \log Z\\, and the derivatives of \\\log Z\\ follow
from the truncated expectations \\\mathbb{E}\_T\[\partial^B f / f\]\\
through the moment-to-cumulant expansion. Each distinct block costs one
quadrature.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list of the parent's parameters.

- expected:

  Logical; if `TRUE`, the expected derivatives.

## Value

A named list of derivative components.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
