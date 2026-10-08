# The Skewness Below Which the Expected Information Is a Series

Returns the bound on \\\|\gamma_1\|\\ below which
[`distrib_expected_hessian.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal2Distrib.md)
and its derivatives come from the series in \\r = (\gamma_1/c)^{1/3}\\
rather than from quadrature. It matches `SN2_GE` in
`src/pt_skewnormal2.h`.

## Usage

``` r
sn2_ge()
```

## Value

A single number.
