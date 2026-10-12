# The Weights of a Multivariate t, Scaled for a Distant Observation

Computes \\w/s\\, \\c = (\nu+p)/s\\ and \\1/s\\, with \\w =
\Sigma^{-1}(y-\mu)\\ and \\s = \nu + q\\, without forming \\q\\, which
overflows for \\\lvert y\rvert\\ above about \\10^{154}\\. Each residual
is divided by \\m = \max(1, \max_j\lvert r_j\rvert)\\ first, and the
powers of \\m\\ are carried separately.

## Usage

``` r
mvt_yscaled(y, pc)
```

## Arguments

- y:

  An \\n \times p\\ numeric matrix.

- pc:

  The result of
  [`mvt_pieces()`](https://statmodels7.github.io/distributions7/reference/mvt_pieces.md).

## Value

A named list with `ws` (an \\n \times p\\ matrix), `cw` and `is`
(numeric vectors of length \\n\\).

## See also

[`distrib_hess_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvStudentTDistrib.md)
and
[`mvt_y_tensor()`](https://statmodels7.github.io/distributions7/reference/mvt_y_tensor.md),
the consumers, and
[`mvt_weights()`](https://statmodels7.github.io/distributions7/reference/mvt_weights.md)
for the unscaled quantities.

## Examples

``` r
d <- mvstudent_t1_distrib(1)
pc <- distributions7:::mvt_pieces(d, list(mu1 = 0, sigma_log_L1 = 0, nu = 5))
distributions7:::mvt_yscaled(matrix(c(1, 1e200)), pc)$cw
#> [1] 1 0
```
