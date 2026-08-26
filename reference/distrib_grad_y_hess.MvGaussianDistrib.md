# Multivariate Gaussian Third Derivative in One Response and Two Parameters

Computes \\\partial^3\ell/\partial
y\\\partial\theta_a\partial\theta_b\\, one \\n \times p\\ matrix per
unordered pair of parameters. The response gradient
\\-\Sigma^{-1}(y-\mu)\\ is linear in the mean, so a pair naming two mean
parameters gives the zero matrix. The mixed pairs and the pure matrix
pairs are \$\$\frac{\partial^3\ell}{\partial
y\\\partial\mu_j\partial\eta_k} = -B_k e_j, \qquad
\frac{\partial^3\ell}{\partial y\\\partial\eta_k\partial\eta_l} =
\Sigma^{-1}A\_{kl}w - \Sigma^{-1}\\\left(A_l\Sigma^{-1}A_k +
A_k\Sigma^{-1}A_l\right)\\w,\$\$ with \\B_k =
\Sigma^{-1}A_k\Sigma^{-1}\\ and \\w = \Sigma^{-1}(y-\mu)\\. Only the
pure matrix pairs carry the observation, through \\w\\; the mixed pairs
repeat one row.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. The two coincide here.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of \\n \times p\\ numeric matrices, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).
Every key naming two mean parameters holds a matrix of zeros.

## Notation

\\\mu\\ is the mean, \\\Sigma\\ the covariance, \\\eta\\ the free vector
of the matrix parametrization, \\A_k\\ and \\A\_{kl}\\ its first and
second derivative arrays, \\e_j\\ the \\j\\-th standard basis vector and
\\\ell\\ the log-density of one observation.

## See also

[`distrib_cross2_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.MvGaussianDistrib.md)
and
[`distrib_hess_y_hess.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.MvGaussianDistrib.md)
for the other two derivatives a marginal criterion reads, and
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
for the generic.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
set.seed(1)
y <- distrib_rng(d, 4, theta)

gh <- distrib_grad_y_hess(d, y, theta)
dim(gh$sigma_log_L1_sigma_L2.1)
#> [1] 4 2

# Two mean parameters give exactly zero, the response gradient being linear
# in the mean.
gh$mu1_mu2
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0
#> [3,]    0    0
#> [4,]    0    0

# A mixed pair repeats one row; a matrix pair does not.
gh$mu1_sigma_L2.1
#>           [,1]      [,2]
#> [1,] 0.9771222 -1.349859
#> [2,] 0.9771222 -1.349859
#> [3,] 0.9771222 -1.349859
#> [4,] 0.9771222 -1.349859
round(gh$sigma_log_L1_sigma_L2.1, 4)
#>         [,1]    [,2]
#> [1,] -1.3789  0.9346
#> [2,]  1.2042 -0.2740
#> [3,] -1.8923  1.2466
#> [4,]  1.7681 -2.3799

# Against a second difference of the response gradient.
h <- 1e-4
f <- function(a, b) {
  t2 <- theta
  t2$sigma_log_L1 <- t2$sigma_log_L1 + a
  t2$sigma_L2.1 <- t2$sigma_L2.1 + b
  distrib_grad_y(d, y, t2)
}
max(abs(gh$sigma_log_L1_sigma_L2.1 -
        (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h)))
#> [1] 2.464924e-08
```
