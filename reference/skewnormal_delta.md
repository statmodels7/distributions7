# The Shape a Skew Normal's Moments Depend On

Returns \\\delta = \alpha/\sqrt{1+\alpha^2}\\ and the product
\\b\delta\\ with \\b = \sqrt{2/\pi}\\. Every moment of a skew normal is
a function of \\b\delta\\, so the four moment methods of the family
share this helper and assemble different combinations of one number.

## Usage

``` r
skewnormal_delta(alpha)
```

## Arguments

- alpha:

  The shape parameter, a numeric vector of any sign. Zero gives \\\delta
  = 0\\ and the symmetric Gaussian case; the magnitude saturates at one
  as the shape grows.

## Value

A named list with `delta` and `bd`, each a numeric vector the length of
`alpha`.

## Details

\\\delta\\ is bounded: it runs over \\(-1, 1)\\ as \\\alpha\\ runs over
the whole line, and it saturates quickly, reaching 0.995 at \\\alpha =
10\\. That bound is what caps the skewness and the kurtosis the family
can reach, and is why the skew t exists.

## See also

[`mean.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.SkewNormal1Distrib.md),
[`variance.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.SkewNormal1Distrib.md),
[`skewness.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewNormal1Distrib.md)
and
[`kurtosis.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewNormal1Distrib.md)
for the four consumers.

## Examples

``` r
# delta saturates at one as the shape grows.
distributions7:::skewnormal_delta(c(0, 1, 10, 1e6))$delta
#> [1] 0.0000000 0.7071068 0.9950372 1.0000000
```
