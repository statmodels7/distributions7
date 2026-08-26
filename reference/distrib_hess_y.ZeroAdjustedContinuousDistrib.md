# Zero-Adjusted Continuous Response Hessian

Returns \\\partial^2\ell/\partial y^2\\, which equals the PARENT's for
\\y \ne 0\\ and is `NaN` at the atom, for the reason
[`distrib_grad_y.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.ZeroAdjustedContinuousDistrib.md)
gives: the log-density jumps there and no derivative exists.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- y:

  A numeric vector of observations. Exactly zero gives `NaN`.

- theta:

  A named list with the parent's parameters followed by `za`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector as long as `y`.

## Notation

\\\pi\\ is the probability of the atom at zero and \\\ell\\ the
log-density of one observation.

## See also

[`distrib_grad_y.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.ZeroAdjustedContinuousDistrib.md)
for the first order,
[`za_y_deriv()`](https://statmodels7.github.io/distributions7/reference/za_y_deriv.md)
for the shared body, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
theta <- list(mu = 1, sigma = 2, za = 0.3)

distrib_hess_y(d, c(-1, 0, 2), theta)
#> [1] -0.25   NaN -0.25

# Away from the atom it is the parent's, which for a gaussian is
# -1 / sigma^2 everywhere.
c(mixed = distrib_hess_y(d, 2, theta), parent = -1 / 2^2)
#>  mixed parent 
#>  -0.25  -0.25 
```
