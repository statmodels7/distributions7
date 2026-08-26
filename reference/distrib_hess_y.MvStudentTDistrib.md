# Multivariate Student t Response Hessian

Computes the second derivative of the log-density in the response,
\$\$\frac{\partial^2 \ell}{\partial y\\ \partial y^\top} =
-c\\\Sigma^{-1} + 2d\\ww^\top, \qquad c = \frac{\nu+p}{\nu+q}, \quad d =
\frac{\nu+p}{(\nu+q)^2},\$\$ one \\p \times p\\ matrix PER OBSERVATION.
The gaussian's response Hessian is one matrix for the whole sample; this
one moves with the observation, and far enough out its rank-one term
makes the curvature positive along the direction of the residual, which
is the redescending score seen at second order.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters, each component a single number.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A \\p \times p \times n\\ numeric array, slice \\i\\ holding
\\\partial^2\ell_i/\partial y_i\partial y_i^\top\\. Compare
[`distrib_hess_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvGaussianDistrib.md),
which returns a \\p \times p\\ matrix and no third dimension at all.

## Notation

\\\Sigma\\ is the scale matrix, \\\nu\\ the degrees of freedom, \\p\\
the dimension, \\q\\ the squared Mahalanobis distance, \\c\\ and \\d\\
the two weights and \\w = \Sigma^{-1}(y-\mu)\\.

## See also

[`distrib_grad_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.MvStudentTDistrib.md)
for the first derivative,
[`distrib_cross2_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.MvStudentTDistrib.md)
for its derivative in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- mvstudent_t_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(1)
y <- distrib_rng(d, 4, theta)

hy <- distrib_hess_y(d, y, theta)
dim(hy)
#> [1] 2 2 4

# Against a numerical Hessian at one observation.
max(abs(hy[, , 1] -
        numDeriv::hessian(function(z) distrib_pdf(d, z, theta, log = TRUE),
                          y[1, ])))
#> [1] 3.413825e-11

# Far out along one coordinate the curvature turns positive, where a
# gaussian's never does.
far <- cbind(0.5 + c(1, 3, 12), -0.3)
vapply(1:3, function(i) distrib_hess_y(d, far, theta)[1, 1, i], numeric(1))
#> [1] -0.82220937  0.11087892  0.04915131
distrib_hess_y(mvgaussian_distrib(2), far, theta[1:5])[1, 1]
#> [1] -1.014155
```
