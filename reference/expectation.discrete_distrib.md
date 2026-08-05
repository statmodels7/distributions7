# Expectation of a Discrete Distribution

Evaluates \\E\[f(Y)\] = \sum_y f(y)\\P(Y = y;\theta)\\ by summation over
the support, every parameter combination in one batched pass: a finite
support is summed exactly in a single matrix evaluation, an infinite one
through
[`series_vec`](https://statmodels7.github.io/numericals7/reference/series_vec.html),
whose rows retire as they converge.

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
