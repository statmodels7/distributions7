# Response Derivative of a Zero-Adjusted Distribution

Evaluates a response derivative of the PARENT away from the atom and
returns `NaN` at it. Away from zero the \\1-\pi\\ factor is a constant
in \\y\\, so the parent's own derivative is exact and needs no
correction; at zero the log-density jumps, \\\log\pi\\ on one side and
\\\log\\(1-\pi)f(y)\\\\ on the other, and no derivative exists.

## Usage

``` r
za_y_deriv(distrib, y, theta, fun)
```

## Arguments

- distrib:

  A zero-adjusted distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, the parent's followed by the atom
  probability.

- fun:

  The parent's response-derivative function, such as
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
  or
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md).
  It is called on the parent with the non-zero observations and the
  parent's parameters, subset to those observations where a parameter
  varies.

## Value

A numeric vector as long as `y`, `NaN` wherever `y == 0` and the
parent's derivative elsewhere.

## Details

The `NaN` is the point of the function. The finite-difference default
inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
would straddle the jump and return a NUMBER for it, which is worse than
an error: nothing downstream would notice.

## Notation

\\f\\ is the parent's density, \\\pi\\ the probability of the atom and
\\\ell\\ the log-density of one observation.

## See also

[`distrib_grad_y.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.ZeroAdjustedContinuousDistrib.md)
and
[`distrib_hess_y.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.ZeroAdjustedContinuousDistrib.md),
the two methods it serves.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
theta <- list(mu = 1, sigma = 2, za = 0.3)

# NaN at the atom, the parent's derivative elsewhere.
distributions7:::za_y_deriv(d, c(-1, 0, 2), theta, distrib_grad_y)
#> [1]  0.50   NaN -0.25
distrib_grad_y(gaussian1_distrib(), c(-1, 2), theta[c("mu", "sigma")])
#> [1]  0.50 -0.25

# The second order takes the same route.
distributions7:::za_y_deriv(d, c(0, 2), theta, distrib_hess_y)
#> [1]   NaN -0.25
```
