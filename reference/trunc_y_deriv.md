# Response Derivative of a Truncated Distribution

Evaluates a response derivative of the PARENT at the observations inside
\\\[L, U\]\\ and returns `NaN` outside. The normalizing constant does
not depend on \\y\\, so inside the interval the parent's derivative is
exact and needs no correction; outside, the log-density is \\-\infty\\
and no derivative exists.

## Usage

``` r
trunc_y_deriv(distrib, y, theta, fun)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

- fun:

  The parent's response-derivative function, such as
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
  or
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md).

## Value

A numeric vector as long as `y`, `NaN` outside \\\[L, U\]\\ and the
parent's derivative inside.

## Details

The `NaN` is the point of the function. The finite-difference default
inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
would straddle the boundary and return a NUMBER for a point outside the
support, which nothing downstream would question.

A parameter varying by observation is subset to the retained
observations before the parent is called, so the parent sees a `theta`
matching the `y` it is given.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`distrib_grad_y.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.TruncatedContinuousDistrib.md)
and
[`distrib_hess_y.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.TruncatedContinuousDistrib.md),
the two methods it serves, and
[`trunc_inside()`](https://statmodels7.github.io/distributions7/reference/trunc_inside.md)
for the test it applies.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

distributions7:::trunc_y_deriv(tn, c(-2, 0, 3), theta, distrib_grad_y)
#> [1]       NaN 0.2083333       NaN
distrib_grad_y(gaussian1_distrib(), 0, theta)
#> [1] 0.2083333

# The second order takes the same route.
distributions7:::trunc_y_deriv(tn, c(0, 3), theta, distrib_hess_y)
#> [1] -0.6944444        NaN
```
