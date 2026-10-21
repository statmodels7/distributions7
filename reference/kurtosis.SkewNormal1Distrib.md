# Excess Kurtosis of the Skew Normal Distribution

The closed form is \$\$\gamma_2 = 2(\pi - 3)\\ \frac{(b\delta)^4}{(1 -
b^2\delta^2)^{2}},\$\$ with \\\delta = \alpha/\sqrt{1+\alpha^2}\\ and
\\b = \sqrt{2/\pi}\\. The fourth power makes it even in the shape and
non-negative, so the excess kurtosis of a skew normal is never below the
Gaussian's, and the bound on \\\delta\\ caps it at about 0.8692. A
sample excess kurtosis beyond that value lies outside the range that a
skew normal can have at any shape.

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

A numeric vector of excess kurtoses, of the length that the recycled
parameters imply, always in \\\[0, 0.8692)\\. Only `alpha` enters the
value, so a setting that varies only the location or the scale repeats
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

# The kurtosis is even in the shape, and zero only at zero.
kurtosis(d, list(mu = 0, sigma = 1, alpha = c(-3, 0, 3)))
#> [1] 0.5097701 0.0000000 0.5097701

# The reachable range stops short of 0.8692.
kurtosis(d, list(mu = 0, sigma = 1, alpha = 1e6))
#> [1] 0.8691773
```
