# Multivariate Gaussian Mixed Response-Parameter Derivatives

Computes \\\partial^2\ell/\partial y\\ \partial\theta_k\\, one \\n
\times p\\ matrix per parameter. With \\w = \Sigma^{-1}(y - \mu)\\ the
response gradient is \\-w\\, and differentiating it in the mean and in
the free values of the matrix parametrization gives \$\$\frac{\partial^2
\ell}{\partial y \\\partial \mu_j} = \Sigma^{-1}e_j, \qquad
\frac{\partial^2 \ell}{\partial y \\\partial \eta_k} = -P_k (y - \mu) =
\Sigma^{-1}A_k w,\$\$ with \\A_k\\ and \\P_k\\ the derivatives of
\\\Sigma\\ and of \\\Sigma^{-1}\\. The mean block is column \\j\\ of
\\\Sigma^{-1}\\ at every observation; the matrix block carries the
observation, and is computed from \\P_k\\ as
[`mv_matrix_pieces()`](https://statmodels7.github.io/distributions7/reference/mv_matrix_pieces.md)
supplies it.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`. The method returns the
  parameter scale and the generic applies the link scale to its result.
  Every link of this family is the identity, so the two scales coincide.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of \\n \times p\\ numeric matrices, one per parameter, in
`distrib@params` order.

## Details

The result has one \\n \times p\\ matrix per parameter, the shape of the
response, so the derivative in \\\theta_k\\ of each row of the response
gradient is read row by row.

Every link of this family is the identity, so `scale = "link"` and
`scale = "parameter"` give the same numbers.

## Notation

\\\mu\\ is the mean, \\\Sigma\\ the covariance, \\\eta\\ the free vector
of the matrix parametrization, \\A_k = \partial\Sigma/\partial\eta_k\\,
\\P_k = \partial\Sigma^{-1}/\partial\eta_k\\, \\e_j\\ the \\j\\-th
standard basis vector and \\\ell\\ the log-density of one observation.

## See also

[`distrib_grad_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.MvGaussianDistrib.md),
whose derivative in the parameters this is,
[`distrib_cross2_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.MvGaussianDistrib.md)
for the next order, and
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
set.seed(1)
y <- distrib_rng(d, 4, theta)

cy <- distrib_cross_y(d, y, theta)
names(cy)
#> [1] "mu1"          "mu2"          "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
dim(cy$sigma_L2.1)
#> [1] 4 2

# The mean block is a row of the inverse covariance, the same at every
# observation.
cy$mu1
#>          [,1]       [,2]
#> [1,] 1.014155 -0.5399435
#> [2,] 1.014155 -0.5399435
#> [3,] 1.014155 -0.5399435
#> [4,] 1.014155 -0.5399435
solve(mv_sigma(d, theta))[1, ]
#>         v1         v2 
#>  1.0141552 -0.5399435 

# The matrix block against a difference of the response gradient.
h <- 1e-5
tp <- theta; tp$sigma_L2.1 <- tp$sigma_L2.1 + h
tm <- theta; tm$sigma_L2.1 <- tm$sigma_L2.1 - h
max(abs(cy$sigma_L2.1 -
        (distrib_grad_y(d, y, tp) - distrib_grad_y(d, y, tm)) / (2 * h)))
#> [1] 1.055289e-11
```
