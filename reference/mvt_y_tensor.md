# Response Tensors of the Multivariate Student t

`mvt_y_tensor()` assembles the third or fourth response derivative of
[`distrib_deriv3_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.MvStudentTDistrib.md)
from the pieces of
[`mvt_dpieces()`](https://statmodels7.github.io/distributions7/reference/mvt_dpieces.md),
one observation at a time. `mvt_y_tensor3_d()` assembles the derivative
of the third one in a parameter, given that parameter's derivatives of
\\N\\, \\s\\, \\w\\ and \\\Omega\\.

## Usage

``` r
mvt_y_tensor(z, order)

mvt_y_tensor3_d(z, N, Na, sa, wa, Oa)
```

## Arguments

- z:

  The result of
  [`mvt_dpieces()`](https://statmodels7.github.io/distributions7/reference/mvt_dpieces.md).

- order:

  `3L` or `4L`.

- N:

  The single number \\\nu + p\\.

- Na:

  The derivative of \\N\\ in the parameter.

- sa:

  A numeric vector of length \\n\\, the derivative of \\s\\.

- wa:

  An \\n \times p\\ matrix, the derivative of \\w\\.

- Oa:

  A \\p \times p\\ matrix, the derivative of \\\Omega\\.

## Value

A numeric array with \\n\\ as its last dimension.

## See also

[`distrib_deriv3_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.MvStudentTDistrib.md),
[`distrib_cross3_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.MvStudentTDistrib.md).

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0, nu = 4)
z <- distributions7:::mvt_dpieces(d, matrix(c(1, -1), 1, 2), theta)
dim(distributions7:::mvt_y_tensor(z, 3L))
#> [1] 2 2 2 1
```
