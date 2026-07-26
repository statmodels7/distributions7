# Truncated Continuous Response Gradient

\\Z\\ does not depend on \\y\\, so inside \\\[\ell, u\]\\ the response
derivative is the parent's. Outside, the log-density is \\-\infty\\ and
no derivative exists, so `NaN` is returned.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

## Value

A numeric vector.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
