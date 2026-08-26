# Which Observations Lie in the Truncated Support

Tests `y >= L & y <= U`, with both endpoints included.
[`trunc_pdf()`](https://statmodels7.github.io/distributions7/reference/trunc_pdf.md)
uses it to set the log-density to \\-\infty\\ outside the interval,
[`trunc_y_deriv()`](https://statmodels7.github.io/distributions7/reference/trunc_y_deriv.md)
to return `NaN` there, and
[`distrib_atoms.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.TruncatedContinuousDistrib.md)
to keep the parent's atoms that survive.

## Usage

``` r
trunc_inside(distrib, y)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- y:

  A numeric vector of observations.

## Value

A logical vector as long as `y`.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`trunc_pdf()`](https://statmodels7.github.io/distributions7/reference/trunc_pdf.md)
and
[`trunc_y_deriv()`](https://statmodels7.github.io/distributions7/reference/trunc_y_deriv.md),
its two callers.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

# Both endpoints count as inside.
distributions7:::trunc_inside(tn, c(-2, -1, 0, 2, 3))
#> [1] FALSE  TRUE  TRUE  TRUE FALSE

# Which is what makes the density positive at them and zero beyond.
distrib_pdf(tn, c(-2, -1, 2, 3), theta)
#> [1] 0.0000000 0.2363006 0.1557790 0.0000000
```
