# Build One of the Two Multivariate Gaussian Families

Implements the body that
[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
and
[`mvgaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
share. The two differ in the class that they instantiate, in the
`inverted` property, in the prefix of their free names, in the tag of
`distrib_name` (`sigma=` or `omega=`) and in what
`params_interpretation` calls the matrix. Every method is registered on
[MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md),
the parent of both.

## Usage

``` r
mvgaussian_build(n_dim, s, inverted, cls)
```

## Arguments

- n_dim:

  The dimension of one observation.

- s:

  The matrix parametrization, or `NULL` for
  `parameters7::log_cholesky(n_dim)`.

- inverted:

  `TRUE` where `s` carries the precision.

- cls:

  The S7 class to instantiate.

## Value

An object of class `cls`.

## See also

[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md),
the constructors that it serves.
