# Build One of the Two Multivariate Student t Families

Builds the object for both
[`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
and
[`mvstudent_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md),
which share this body. The two differ in the class that they
instantiate, in the `inverted` property, in the prefix of their free
names, in the tag of `distrib_name` (`sigma=` or `omega=`) and in what
`params_interpretation` calls the matrix. Every method is registered on
[MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md),
the parent of both.

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
the constructors that it serves.
