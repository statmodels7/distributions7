# Beta-Binomial Analytical Expected Hessian

The observed Hessian averaged against the mass over \\\\0, \dots, n\\\\.
The support being finite, the expectation is an exact sum rather than a
quadrature, and `approx` is therefore ignored.

## Arguments

- distrib:

  A `BetaBinomDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is an exact sum.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of expected second-derivative components.

## See also

[`betabinom_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom_distrib.md)
