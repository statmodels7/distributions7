# Truncated Continuous Response Gradient

Returns \\\partial\ell_T/\partial y\\, which inside \\\[L, U\]\\ is the
PARENT's, \\Z\\ not depending on \\y\\. Outside the interval the
log-density is \\-\infty\\ and no derivative exists, so `NaN` is
returned where a finite difference would invent a number.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- y:

  A numeric vector of observations. A point outside \\\[L, U\]\\ gives
  `NaN`.

- theta:

  A named list of the parent's parameters.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector as long as `y`.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`trunc_y_deriv()`](https://statmodels7.github.io/distributions7/reference/trunc_y_deriv.md)
for the shared body,
[`distrib_hess_y.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.TruncatedContinuousDistrib.md)
for the second order, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

distrib_grad_y(tn, c(-2, 0, 1, 3), theta)
#> [1]        NaN  0.2083333 -0.4861111        NaN

# Inside the interval it is the parent's own.
distrib_grad_y(gaussian1_distrib(), c(0, 1), theta)
#> [1]  0.2083333 -0.4861111
```
