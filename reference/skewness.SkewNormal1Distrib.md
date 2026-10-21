# Skewness of the Skew Normal Distribution

The closed form is \$\$\gamma_1 = \frac{4-\pi}{2}\\
\frac{(b\delta)^3}{(1 - b^2\delta^2)^{3/2}},\$\$ with \\\delta =
\alpha/\sqrt{1+\alpha^2}\\ and \\b = \sqrt{2/\pi}\\. It takes the sign
of the shape and vanishes at zero. Because \\\delta\\ is bounded by 1,
the skewness that the family can reach is bounded too, by about 0.9953
in absolute value. A sample skewness beyond that value lies outside the
range that a skew normal can have at any shape; the skewness of the skew
t has no such bound.

## Arguments

- x:

  A `SkewNormal1Distrib`, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- theta:

  A named list with components `mu`, `sigma` (positive) and `alpha` (any
  sign); each is a numeric vector of length 1 or `n`. Only `alpha`
  enters the value.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of the length that the recycled parameters imply,
always inside \\(-0.9953, 0.9953)\\. Only `alpha` enters the value, so a
setting that varies only the location or the scale repeats one number.

## Notation

\\\alpha\\ is the shape, \\\delta = \alpha/\sqrt{1+\alpha^2}\\ and \\b =
\sqrt{2/\pi}\\. The location and the scale do not enter a standardized
moment.

## See also

[`kurtosis.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewNormal1Distrib.md),
bounded for the same reason;
[`skewness.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewTDistrib.md),
which is not bounded;
[`skewnormal_delta()`](https://statmodels7.github.io/distributions7/reference/skewnormal_delta.md),
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

## Examples

``` r
d <- skewnormal1_distrib()

# Zero at shape zero, and it takes the sign of the shape.
skewness(d, list(mu = 0, sigma = 1, alpha = c(-3, 0, 3)))
#> [1] -0.6670236  0.0000000  0.6670236

# The reachable range stops short of one.
skewness(d, list(mu = 0, sigma = 1, alpha = 1e6))
#> [1] 0.9952717
```
