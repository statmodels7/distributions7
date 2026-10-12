# The Quadratic Form of a Residual, Through the Chart's Factor

Computes \\q_i = r_i^\top \Sigma^{-1} r_i\\ for each row of `r` as the
squared norm of a triangular solve or product with the matrix
parametrization's own factor: \\\lVert L^{-1} r_i\rVert^2\\ where the
parametrization carries \\\Sigma = L L^\top\\, and \\\lVert L^\top
r_i\rVert^2\\ where it carries the precision \\\Omega = L L^\top\\.

## Usage

``` r
mv_quad(r, s, eta, inverted)
```

## Arguments

- r:

  An \\n \times p\\ numeric matrix of residuals.

- s:

  A parameters7 matrix parametrization of full rank.

- eta:

  Its free vector.

- inverted:

  Logical of length 1, `TRUE` when `s` carries the precision.

## Value

A numeric vector of length \\n\\.

## Details

Near a singular \\\Sigma\\ the inverse has large entries of nearly rank
one, and \\r^\top \Sigma^{-1} r\\ formed with the explicit inverse loses
the digits that an observation close to the degenerate direction has: at
a log-Cholesky free value \\\log L\_{22} = -11.5\\ the relative error
was of order \\10^{-7}\\. The squared norm adds positive terms only.

## See also

[`mv_sampling_factor()`](https://statmodels7.github.io/distributions7/reference/mv_sampling_factor.md)
for the same factor used to draw, and
[`mvt_weights()`](https://statmodels7.github.io/distributions7/reference/mvt_weights.md)
and
[`distrib_pdf.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.MvGaussianDistrib.md)
for consumers.

## Examples

``` r
s <- parameters7::log_cholesky(2)
eta <- c(0.1, -0.2, 0.4)
r <- rbind(c(1, -1), c(0.5, 2))
distributions7:::mv_quad(r, s, eta, FALSE)
#> [1] 3.585867 5.140951
rowSums((r %*% solve(parameters7::param_value(s, eta))) * r)
#> [1] 3.585867 5.140951
```
