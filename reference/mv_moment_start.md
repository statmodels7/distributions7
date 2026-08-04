# The Moment Estimates a Multivariate Family Starts From

Returns the sample mean and the sample covariance of a multivariate
response, with the covariance made safely positive definite.

## Usage

``` r
mv_moment_start(y, p)
```

## Arguments

- y:

  The response, an \\n \times p\\ matrix.

- p:

  The dimension.

## Value

A list with `mu` and `sigma`.

## Details

A sample covariance is singular when there are fewer observations than
coordinates, and nearly singular when two coordinates almost coincide.
Either way a structure cannot be inverted onto it, so the eigenvalues
are floored at a small multiple of the largest before the matrix is
handed on. The floor moves a starting value, which is allowed to be
approximate; it would not be allowed anywhere the answer is reported.
