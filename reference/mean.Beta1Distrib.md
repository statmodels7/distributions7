# Mean of the Beta Distribution

Closed form: \\E\[Y\] = \mu\\. This parametrization carries the mean and
a precision, so the mean is a read; the two shapes of the standard
parametrization are recovered as \\a = \mu\phi\\ and \\b =
(1-\mu)\phi\\, and they sum to \\\phi\\.

## Arguments

- x:

  A `Beta1Distrib`, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- theta:

  A named list with components `mu` (the mean, strictly between 0 and 1)
  and `phi` (the precision, positive), each a numeric vector of length 1
  or `n`. Aligned and validated by name, so a mean at or outside the
  unit interval throws.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$phi))`, strictly inside \\(0,1)\\.

## Notation

\\\mu \in (0,1)\\ is the mean, \\\phi \> 0\\ the precision, and \\a =
\mu\phi\\, \\b = (1-\mu)\phi\\ the two shapes.

## See also

[`variance.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Beta1Distrib.md),
where the precision does enter;
[`skewness.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Beta1Distrib.md);
[`mean.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Beta2Distrib.md)
for the two-shape parametrization;
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

## Examples

``` r
d <- beta1_distrib()

# The first parameter is the mean, and the precision does not move it.
mean(d, list(mu = c(0.2, 0.5, 0.8), phi = 5))
#> [1] 0.2 0.5 0.8
```
