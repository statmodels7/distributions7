# Zero-Adjusted Continuous Response Gradient

Returns \\\partial\ell/\partial y\\, which equals the PARENT's for \\y
\ne 0\\, the factor \\1-\pi\\ not depending on \\y\\. At the atom the
log-density jumps and no derivative exists, so `NaN` is returned. That
is deliberate: the numerical fallback would straddle the jump and hand
back a number nothing downstream would question.

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

[`distrib_hess_y.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.ZeroAdjustedContinuousDistrib.md)
for the second order,
[`za_y_deriv()`](https://statmodels7.github.io/distributions7/reference/za_y_deriv.md)
for the shared body, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
theta <- list(mu = 1, sigma = 2, za = 0.3)

distrib_grad_y(d, c(-1, 0, 2), theta)
#> [1]  0.50   NaN -0.25

# Away from the atom it is the parent's own.
distrib_grad_y(gaussian1_distrib(), c(-1, 2), theta[c("mu", "sigma")])
#> [1]  0.50 -0.25
```
