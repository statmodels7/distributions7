# The Center of the Location-Scale Rule

Returns, for each distinct shape, the standardized point at which
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md)
centers its rule, the point where the nodes cluster. It is zero except
for the centered skew normal, whose density has its step at \\w = -r\\
with \\r = \\2\|\gamma_1\|/(4-\pi)\\^{1/3}\\ \mathrm{sign}(\gamma_1)\\;
at a skewness near the ceiling the step is too sharp for nodes clustered
at zero. The compiled diagonal of that family uses the same center.

## Usage

``` r
loc_scale_center(distrib, U)
```

## Arguments

- distrib:

  A distribution object.

- U:

  The matrix of distinct shapes, one row per shape.

## Value

A numeric vector with one entry per row of `U`.
