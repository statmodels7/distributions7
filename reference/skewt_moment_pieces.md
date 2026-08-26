# The Quantities a Skew t's Moments Are Built From

Returns the four pieces the family's four moment methods share: \\\delta
= \alpha/\sqrt{1+\alpha^2}\\, the constant \\b\_\nu =
\sqrt{\nu/\pi}\\\Gamma\\(\nu-1)/2\\/\Gamma(\nu/2)\\, the mean \\\mu_z =
\delta b\_\nu\\ of the standardized variable and its variance
\\\sigma_z^2 = \nu/(\nu-2) - \mu_z^2\\.

## Usage

``` r
skewt_moment_pieces(alpha, nu)
```

## Arguments

- alpha:

  The shape parameter, a numeric vector of any sign.

- nu:

  The degrees of freedom, a positive numeric vector. Values at or below
  1 give `NaN` in `bnu` and `mz`; values at or below 2 give `NaN` in
  `vz` as well.

## Value

A named list with `delta`, `bnu`, `mz` and `vz`, each a numeric vector
recycled to the longer of `alpha` and `nu`.

## Details

The gamma ratio in \\b\_\nu\\ is formed as
`exp(lgamma((nu - 1) / 2) - lgamma(nu / 2))`, so it stays finite at
degrees of freedom where the two gamma functions themselves would
overflow.

Existence is enforced here, once for all four callers: \\b\_\nu\\ is
`NaN` at \\\nu \le 1\\ and \\\sigma_z^2\\ is `NaN` at \\\nu \le 2\\, so
a moment built from either inherits the `NaN` and no method has to test
the threshold twice.

## See also

[`mean.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.SkewTDistrib.md),
[`variance.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.SkewTDistrib.md),
[`skewness.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewTDistrib.md)
and
[`kurtosis.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewTDistrib.md)
for the four consumers;
[`skewnormal_delta()`](https://statmodels7.github.io/distributions7/reference/skewnormal_delta.md)
for the same shape factor in the Gaussian-tailed family.

## Examples

``` r
# Below two degrees of freedom the variance piece is NaN.
distributions7:::skewt_moment_pieces(alpha = 2, nu = c(0.5, 1.5, 5))$vz
#> [1]       NaN       NaN 0.9461605
```
