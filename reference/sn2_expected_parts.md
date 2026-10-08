# The Expected Information of the Centered Skew Normal, by Region

The expected information (`k = 0`) or its derivatives (`k = 1, 2`): from
the series kernels where \\\|\gamma_1\| \<\\
[`sn2_ge()`](https://statmodels7.github.io/distributions7/reference/sn2_ge.md),
from
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md)
elsewhere, each observation by its own skewness.

## Usage

``` r
sn2_expected_parts(distrib, theta, k, n, threads = 1L)
```

## Arguments

- distrib:

  A `SkewNormal2Distrib` object.

- theta:

  The parameters.

- k:

  The order: 0, 1 or 2.

- n:

  The number of observations.

- threads:

  Passed to
  [`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md).

## Value

A named list of numeric vectors of length `n`.
