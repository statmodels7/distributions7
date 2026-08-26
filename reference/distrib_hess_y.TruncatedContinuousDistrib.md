# Truncated Continuous Response Hessian

Returns \\\partial^2\ell_T/\partial y^2\\, which inside \\\[L, U\]\\ is
the PARENT's and is `NaN` outside, for the reason
[`distrib_grad_y.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.TruncatedContinuousDistrib.md)
gives: the log-density is \\-\infty\\ there and has no derivative.

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
[`distrib_grad_y.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.TruncatedContinuousDistrib.md)
for the first order, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

distrib_hess_y(tn, c(-2, 0, 1, 3), theta)
#> [1]        NaN -0.6944444 -0.6944444        NaN

# For a gaussian it is -1 / sigma^2 wherever it exists.
c(truncated = distrib_hess_y(tn, 1, theta), parent = -1 / 1.2^2)
#>  truncated     parent 
#> -0.6944444 -0.6944444 
```
