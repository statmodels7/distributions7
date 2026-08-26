# Skewness of the Skew Normal Distribution

Closed form: \$\$\gamma_1 = \frac{4-\pi}{2}\\ \frac{(b\delta)^3}{(1 -
b^2\delta^2)^{3/2}},\$\$ with \\\delta = \alpha/\sqrt{1+\alpha^2}\\ and
\\b = \sqrt{2/\pi}\\. It takes the sign of the shape and vanishes at
zero. Because \\\delta\\ is bounded by 1, the skewness the family can
reach is bounded too, by about 0.9953 in absolute value: a sample skewed
more than that cannot be fitted by a skew normal at any shape, and the
skew t is the family to reach for.

## Arguments

- x:

  A `SkewNormal1Distrib`, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- theta:

  A named list with components `mu`, `sigma` (positive) and `alpha` (any
  sign), each a numeric vector of length 1 or `n`. Only `alpha` enters
  the value.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of the length the recycled parameters imply, always
inside \\(-0.9953, 0.9953)\\. Only `alpha` enters the value, so a
setting that varies the location or the scale alone repeats one number.

## Notation

\\\alpha\\ is the shape, \\\delta = \alpha/\sqrt{1+\alpha^2}\\ and \\b =
\sqrt{2/\pi}\\. Neither the location nor the scale enters a standardized
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
