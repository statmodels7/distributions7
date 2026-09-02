# Build One of the Two Multivariate Student t Families

The body
[`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
and
[`mvstudent_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
share. The two differ in the class they instantiate, in the prefix their
free names carry and in what `params_interpretation` calls the matrix,
and in nothing else: every method is registered on
[MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md),
the parent both inherit.

## Usage

``` r
mvstudent_t_build(n_dim, s, link_nu, inverted, cls)
```

## Arguments

- n_dim:

  The dimension of one observation.

- s:

  The matrix parametrization, or `NULL` for
  `parameters7::log_cholesky(n_dim)`.

- link_nu:

  The link on the degrees of freedom.

- inverted:

  `TRUE` where `s` carries the inverse scale matrix.

- cls:

  The S7 class to instantiate.

## Value

An object of class `cls`.

## See also

[`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md),
the constructors this serves.
