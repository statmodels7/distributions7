# Variance of the Generalized Pareto Distribution

Closed form: \\\operatorname{Var}(Y) = \sigma^2/\\(1-\xi)^2(1-2\xi)\\\\
for \\\xi \< 1/2\\, and `Inf` at or above one half. The threshold is
half the mean's, the second moment being the first to fail as the tail
grows heavy, and a fitted object at \\\xi = 0.6\\ reports a finite mean
beside an infinite variance.

## Arguments

- x:

  A `GPDDistrib`, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- theta:

  A named list with components `sigma` (positive) and `xi` (any real
  value), each a numeric vector of length 1 or `n`. Settings with \\\xi
  \ge 1/2\\ give `Inf`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$sigma), length(theta$xi))`, `Inf` wherever \\\xi \ge
1/2\\.

## Notation

\\\sigma \> 0\\ is the scale and \\\xi\\ the shape, of either sign.

## See also

[`mean.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.GPDDistrib.md),
whose threshold is 1;
[`skewness.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.GPDDistrib.md),
whose threshold is \\1/3\\;
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

## Examples

``` r
d <- gpd_distrib()

# Finite below shape one half and infinite at or above it.
variance(d, list(sigma = 1, xi = c(0, 0.4, 0.5, 0.6)))
#> [1]  1.00000 13.88889      Inf      Inf

# A mean without a variance, which is what the two thresholds allow.
c(mean = mean(d, list(sigma = 1, xi = 0.6)),
  variance = variance(d, list(sigma = 1, xi = 0.6)))
#>     mean variance 
#>      2.5      Inf 
```
