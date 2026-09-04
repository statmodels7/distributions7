# Multivariate Student t Third Derivative in Two Responses and One Parameter

Computes \\\partial^3\ell/\partial y\\\partial y^\top
\partial\theta_a\\, one \\p \times p \times n\\ array per parameter.
Writing \\M = \ell^{(yy)} = -c\\\Sigma^{-1} + 2d\\ww^\top\\ with \\c =
(\nu+p)/s\\, \\d = (\nu+p)/s^2\\ and \\s = \nu+q\\, \$\$\partial_a M =
-c_a\Sigma^{-1} - c\\\partial_a\Sigma^{-1} + 2d_a ww^\top +
2d\left(w_aw^\top + ww_a^\top\right),\$\$ every piece coming from
[`mvt_dpieces()`](https://statmodels7.github.io/distributions7/reference/mvt_dpieces.md).

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. The two differ here.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of \\p \times p \times n\\ numeric arrays, one per
parameter, in `distrib@params` order.

## Details

The shape differs from the gaussian's, and the difference is the
family's defining property. There the response Hessian is
\\-\Sigma^{-1}\\, one matrix for the whole sample, so its derivative is
one matrix per parameter. Here it depends on the observation through the
weight, so this method and
[`distrib_hess_y_hess.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.MvStudentTDistrib.md)
return ONE MATRIX PER ROW. Nor does any component vanish: the gaussian's
location components are exactly zero and none of these is.

This is one of the three derivatives a marginal criterion reads when the
family stands as a prior over a coefficient block.

## Notation

\\\Sigma\\ is the scale matrix, \\\nu\\ the degrees of freedom, \\q\\
the squared Mahalanobis distance, \\s = \nu+q\\, \\c\\ and \\d\\ the two
weights, \\w = \Sigma^{-1}(y-\mu)\\, and a subscript \\a\\ denotes a
derivative in parameter \\a\\.

## See also

[`distrib_hess_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvStudentTDistrib.md),
whose derivative in the parameters this is,
[`distrib_grad_y_hess.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.MvStudentTDistrib.md)
and
[`distrib_hess_y_hess.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.MvStudentTDistrib.md)
for the other two a marginal criterion reads,
[`mvt_dpieces()`](https://statmodels7.github.io/distributions7/reference/mvt_dpieces.md)
for the assembly, and
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(1)
y <- distrib_rng(d, 4, theta)

c2 <- distrib_cross2_y(d, y, theta)
dim(c2$sigma_L2.1)
#> [1] 2 2 4

# One matrix per row, so the location component is not zero as it is for a
# gaussian.
c2$mu1[, , 1]
#>            [,1]       [,2]
#> [1,]  0.7080609 -0.3831362
#> [2,] -0.3831362  0.5149049
distrib_cross2_y(mvgaussian1_distrib(2), y, theta[1:5])$mu1
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0

# Against a difference of the response Hessian.
h <- 1e-5
vapply(seq_along(d@params), function(k) {
  tp <- theta; tp[[k]] <- tp[[k]] + h
  tm <- theta; tm[[k]] <- tm[[k]] - h
  max(abs(c2[[k]] -
          (distrib_hess_y(d, y, tp) - distrib_hess_y(d, y, tm)) / (2 * h)))
}, numeric(1))
#> [1] 2.428202e-11 6.885559e-11 1.211347e-10 1.811258e-10 2.061795e-10
#> [6] 1.042630e-11
```
