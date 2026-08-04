# Skew t Third-Order Derivatives

Third-order derivatives assembled so that no finite difference is ever
applied to another finite difference. Components whose Hessian entry is
closed form – both indices in \\(\mu, \sigma, \alpha)\\, or one index
equal to \\\nu\\ with the stencil taken along a different variable –
come from the generic construction, one stencil on an analytic quantity.
The components the generic construction would nest are replaced: \\(i,
\nu, \nu)\\ is one five-point second-difference of the closed-form score
component \\i\\, and \\(\nu, \nu, \nu)\\ is one five-point
third-difference of the log-density itself. The derivative of a Student
t distribution function in its degrees of freedom has no elementary
form, so this is the same obstruction, and the same remedy, as the
Hessian's.

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
  [`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md).

- nsim:

  Monte Carlo sample size when `approx = "mc"`.

## Value

A named list of third-derivative component vectors.

## See also

[`skewt_distrib`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
