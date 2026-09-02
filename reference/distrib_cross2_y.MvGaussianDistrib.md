# Multivariate Gaussian Third Derivative in Two Responses and One Parameter

Computes \\\partial^3\ell/\partial y\\\partial y^\top
\partial\theta_k\\, one \\p \times p\\ matrix per parameter. The
response Hessian is \\-\Sigma^{-1}\\, which does not involve the mean,
so every mean component is exactly the zero matrix; the matrix
components are \$\$\frac{\partial^3\ell}{\partial y\\\partial
y^\top\partial\eta_k} = B_k = \Sigma^{-1}A_k\Sigma^{-1},\$\$ with \\A_k
= \partial\Sigma/\partial\eta_k\\. No component depends on the
observation, so one matrix is returned per parameter and `y` is not read
at all.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. Not read: no
  component of this derivative depends on the response.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. The two coincide here.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of \\p \times p\\ numeric matrices, one per parameter, in
`distrib@params` order. The \\p\\ mean components are the zero matrix.

## Details

This is one of the three derivatives a marginal criterion reads when the
family stands as a prior over a coefficient block, the other two being
[`distrib_hess_y_hess.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.MvGaussianDistrib.md)
and
[`distrib_grad_y_hess.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.MvGaussianDistrib.md).
Every link of this family is the identity, so `scale = "link"` and
`scale = "parameter"` give the same numbers.

## Notation

\\\Sigma\\ is the covariance, \\\eta\\ the free vector of the matrix
parametrization, \\A_k = \partial\Sigma/\partial\eta_k\\ and \\\ell\\
the log-density of one observation.

## See also

[`distrib_hess_y_hess.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.MvGaussianDistrib.md)
and
[`distrib_grad_y_hess.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.MvGaussianDistrib.md)
for the two fourth- and third-order siblings,
[`distrib_hess_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvGaussianDistrib.md),
whose derivative in the parameters this is, and
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
y <- matrix(0, 4, 2)

c2 <- distrib_cross2_y(d, y, theta)
names(c2)
#> [1] "mu1"          "mu2"          "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  

# Every mean component is the zero matrix, the response Hessian carrying
# no mean at all.
c2$mu1
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0

# A matrix component against a difference of the response Hessian.
h <- 1e-5
tp <- theta; tp$sigma_L2.1 <- tp$sigma_L2.1 + h
tm <- theta; tm$sigma_L2.1 <- tm$sigma_L2.1 - h
max(abs(c2$sigma_L2.1 -
        (distrib_hess_y(d, y, tp) - distrib_hess_y(d, y, tm)) / (2 * h)))
#> [1] 3.452794e-12
```
