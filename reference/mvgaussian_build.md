# Build One of the Two Multivariate Gaussian Families

The body
[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
and
[`mvgaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
share. The two differ in the class they instantiate, in the prefix their
free names carry and in what `params_interpretation` calls the matrix,
and in nothing else: every method is registered on
[MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md),
the parent both inherit.

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
the constructors this serves.
