# Mean of the Pseudo-Huber Distribution

Closed form, replacing the numerical default: \\E\[Y\] = \mu\\. The
density is symmetric about \\\mu\\ at every scale and shape, so the
location parameter is the mean, the median and the mode at once.

## Arguments

- x:

  A `PseudoHuberDistrib`, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- theta:

  A named list with components `mu` (the location), `sigma` (the scale,
  positive) and `nu` (the shape, positive), each a numeric vector of
  length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length equal to the longest of the three
components.

## See also

[`variance.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.PseudoHuberDistrib.md),
which does depend on the shape;
[`skewness.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.PseudoHuberDistrib.md),
zero by the same symmetry;
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
for the family.

## Examples

``` r
d <- pseudohuber_distrib()

# The location is the mean, whatever the scale and the shape.
mean(d, list(mu = c(-1, 0, 3), sigma = 2, nu = 4))
#> [1] -1  0  3
```
