# Raw Moments of a Generalized Gamma

Returns the first four raw moments \\E\[Y^k\] =
a^k\\\Gamma\\(d+k)/p\\/\Gamma(d/p)\\ for \\k = 1, \ldots, 4\\. The four
moment methods of the family share this helper and assemble different
combinations of the four values.

## Usage

``` r
gengamma_raw_moments(a, d, p)
```

## Arguments

- a:

  The scale parameter, a positive numeric vector.

- d:

  The first shape parameter, a positive numeric vector.

- p:

  The second shape parameter, a positive numeric vector. Small `p`
  pushes \\(d+k)/p\\ to large arguments, where the log scale earns its
  keep.

## Value

A list of four numeric vectors, the raw moments of order 1 to 4, each
recycled to the longest of `a`, `d` and `p`.

## Details

The gamma ratio is formed on the log scale, as
`exp(k * log(a) + lgamma((d + k) / p) - lgamma(d / p))`, so it stays
finite at shapes where either gamma function on its own would overflow.

## See also

[`mean.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.GenGamma1Distrib.md),
[`variance.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.GenGamma1Distrib.md),
[`skewness.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.GenGamma1Distrib.md)
and
[`kurtosis.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.GenGamma1Distrib.md)
for the four consumers.

## Examples

``` r
# At d = p the family is Weibull, and the first raw moment is its mean.
distributions7:::gengamma_raw_moments(2, 3, 3)[[1]]
#> [1] 1.785959
mean(weibull1_distrib(), list(mu = 2, sigma = 3))
#> [1] 1.785959
```
