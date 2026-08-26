# The Moment Estimates a Multivariate Family Starts From

Returns the sample mean and the sample covariance of a multivariate
response, with the covariance made safely positive definite. The
covariance divides by \\n\\ rather than \\n - 1\\, which is the maximum
likelihood estimator and therefore the point a Gaussian fit is looking
for.

## Usage

``` r
mv_moment_start(y, p)
```

## Arguments

- y:

  The response, an \\n \times p\\ matrix, or anything
  [`base::as.matrix()`](https://rdrr.io/r/base/matrix.html) turns into
  one.

- p:

  The dimension, a single positive integer. Used only for the one-row
  fallback.

## Value

A list with two components: `mu`, a numeric vector of length \\p\\, and
`sigma`, a \\p \times p\\ positive definite matrix.

## Details

A sample covariance is singular when there are fewer observations than
coordinates, and nearly singular when two coordinates almost coincide.
Either way a `parameters7` structure cannot be inverted onto it, so the
eigenvalues are floored at a small multiple of the largest before the
matrix is handed on. The floor moves a starting value, which is allowed
to be approximate; it would not be allowed anywhere the answer is
reported.

A response with one row has no covariance at all, and the identity is
returned in its place.

## See also

[`distrib_start.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_start.MvGaussianDistrib.md)
and
[`distrib_start.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_start.MvStudentTDistrib.md),
the two callers;
[`param_free_or_fit()`](https://statmodels7.github.io/distributions7/reference/param_free_or_fit.md),
which carries the matrix onto a structure.
