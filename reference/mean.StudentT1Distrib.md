# Mean of the Student t Distribution

The mean has the closed form \\E\[Y\] = \mu\\ for \\\nu \> 1\\, and is
`NaN` at or below one, where the first moment does not exist. At \\\nu =
1\\ the family is the Cauchy, whose divergence is the reason for the
threshold; above it the density is symmetric about \\\mu\\ and the
location is the mean.

## Arguments

- x:

  A `StudentT1Distrib`, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- theta:

  A named list with components `mu` (the location), `sigma` (the scale,
  positive) and `nu` (the degrees of freedom, positive), in which each
  component is a numeric vector of length 1 or `n`. Settings with \\\nu
  \le 1\\ give `NaN` in their own positions and do not affect the
  others.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, whose length is the one that the recycled
parameters imply, with `NaN` wherever \\\nu \le 1\\.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale and \\\nu \> 0\\ the
degrees of freedom.

## See also

[`variance.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.StudentT1Distrib.md),
whose threshold is 2;
[`mean.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.CauchyDistrib.md),
the \\\nu = 1\\ case;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

## Examples

``` r
d <- student_t1_distrib()

# The location is the mean above one degree of freedom, and NaN at or below.
mean(d, list(mu = 3, sigma = 1, nu = c(0.5, 1, 2)))
#> [1] NaN NaN   3
```
