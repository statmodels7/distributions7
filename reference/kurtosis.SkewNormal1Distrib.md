# Excess Kurtosis of the Skew Normal Distribution

Closed form: \$\$\gamma_2 = 2(\pi - 3)\\ \frac{(b\delta)^4}{(1 -
b^2\delta^2)^{2}},\$\$ with \\\delta = \alpha/\sqrt{1+\alpha^2}\\ and
\\b = \sqrt{2/\pi}\\. The fourth power makes it even in the shape and
non-negative, so a skew normal is never lighter-tailed than a Gaussian,
and \\\delta\\'s bound caps it at about 0.8692. A sample needing more
excess kurtosis than that is outside the family whatever the shape.

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

A numeric vector of excess kurtoses, of the length the recycled
parameters imply, always in \\\[0, 0.8692)\\. Only `alpha` enters the
value, so a setting that varies the location or the scale alone repeats
one number.

## Notation

\\\alpha\\ is the shape, \\\delta = \alpha/\sqrt{1+\alpha^2}\\ and \\b =
\sqrt{2/\pi}\\.

## See also

[`skewness.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewNormal1Distrib.md),
bounded for the same reason;
[`kurtosis.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewTDistrib.md),
which is not bounded;
[`skewnormal_delta()`](https://statmodels7.github.io/distributions7/reference/skewnormal_delta.md),
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

## Examples

``` r
d <- skewnormal1_distrib()

# Even in the shape, and zero only at zero.
kurtosis(d, list(mu = 0, sigma = 1, alpha = c(-3, 0, 3)))
#> [1] 0.5097701 0.0000000 0.5097701

# The reachable range stops short of 0.8692.
kurtosis(d, list(mu = 0, sigma = 1, alpha = 1e6))
#> [1] 0.8691773
```
