# Expectation of a Discrete Distribution

Evaluates \\E\[f(Y)\] = \sum_y f(y)\\P(Y = y;\theta)\\ by direct
summation over the support, truncating the series once the accumulated
tail contribution falls below tolerance.

## Arguments

- distrib:

  A `discrete_distrib`.

- f:

  The function whose expectation is taken.

- theta:

  A named list of parameters.

- ...:

  Further arguments passed to `f`.

## Value

A numeric vector of expected values.
