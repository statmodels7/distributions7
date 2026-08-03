# Expectation of a Continuous Distribution

Evaluates \\E\[f(Y)\] = \int f(y)\\p(y;\theta)\\dy\\ by adaptive
quadrature ([`integrate`](https://rdrr.io/r/stats/integrate.html)). The
domain is split at the 0.1, 0.5 and 0.9 quantiles of the distribution
and each panel is integrated separately, which anchors the quadrature
nodes on the probability mass wherever it sits; the panels are then
summed.

## Arguments

- distrib:

  A `continuous_distrib`.

- f:

  The function whose expectation is taken.

- theta:

  A named list of parameters.

- ...:

  Further arguments passed to `f`.

## Value

A numeric vector of expected values.
