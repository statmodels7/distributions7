# Skew t Fourth-Order Derivatives

Fourth-order derivatives assembled with the discipline of
[`distrib_deriv3.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewTDistrib.md):
the generic construction serves every component whose Hessian entry is
closed form, and the ones it would nest are replaced by one stencil each
– \\(i, \nu, \nu, \nu)\\ by a third-difference of the closed-form score
component \\i\\, and \\(\nu, \nu, \nu, \nu)\\ by a fourth-difference of
the log-density. The pure-\\\nu\\ component is the least accurate
quantity the family reports, at roughly four significant digits.

## Arguments

- distrib:

  A `SkewTDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma`, `alpha` and `nu`.

- expected:

  Logical; if `TRUE`, the expectation is approximated numerically.

- approx:

  Strategy for the expectation; see
  [`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md).

- nsim:

  Monte Carlo sample size when `approx = "mc"`.

## Value

A named list of fourth-derivative component vectors.

## See also

[`skewt_distrib`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
