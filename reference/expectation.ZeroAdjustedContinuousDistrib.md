# Expectation for Zero-Adjusted Continuous Distributions

The zero-adjusted continuous distribution has a point mass at 0 that
plain numerical integration would miss. The expectation decomposes
exactly as \$\$E\[f(Y)\] = \pi f(0) + (1-\pi) E_W\[f(W)\]\$\$

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object.

- f:

  A function `f(y, theta, ...)` (receives the full theta, including
  `za`).

- theta:

  A list with the parent's parameters followed by `za`.

- ...:

  Additional arguments passed to `f`.

## Value

A numeric scalar, the expectation of `f` under the distribution.
