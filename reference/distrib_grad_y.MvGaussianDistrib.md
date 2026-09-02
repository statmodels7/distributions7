# Multivariate Gaussian Response Gradient

Computes the derivative of the log-density in the response,
\$\$\frac{\partial \ell}{\partial y} = -\Sigma^{-1}(y - \mu),\$\$ one
row per observation. This is the whitened residual with a minus sign, so
it is exactly the negative of the score in the mean: for a location
family the two derivatives differ only in sign, and the method returns
the same numbers
[`distrib_gradient.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvGaussianDistrib.md)
gives for the mean components.

The shape is what separates this from the univariate case. There the
derivative in a scalar response is a numeric vector of length \\n\\;
here it is an \\n \times p\\ matrix, and a consumer written for a vector
will recycle it silently.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation and gives a one-row result.

- theta:

  A named list of parameters, each component a single number.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

An \\n \times p\\ numeric matrix, row \\i\\ holding
\\\partial\ell_i/\partial y_i\\.

## See also

[`distrib_hess_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvGaussianDistrib.md)
for the second derivative,
[`distrib_cross_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.MvGaussianDistrib.md)
for the mixed block, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
set.seed(1)
y <- distrib_rng(d, 4, theta)

distrib_grad_y(d, y, theta)
#>            [,1]       [,2]
#> [1,]  0.7125038 -0.4024617
#> [2,] -0.5288705  1.0021223
#> [3,]  0.9715850 -0.5953472
#> [4,] -1.1170798 -0.9017918

# Against a numerical derivative taken row by row.
num <- t(apply(y, 1, function(r)
  numDeriv::grad(function(z) distrib_pdf(d, z, theta, log = TRUE), r)))
max(abs(distrib_grad_y(d, y, theta) - num))
#> [1] 5.809908e-11

# And it is minus the score in the mean, this being a location family.
g <- distrib_gradient(d, y, theta)
all.equal(distrib_grad_y(d, y, theta), -cbind(g$mu1, g$mu2),
          check.attributes = FALSE)
#> [1] TRUE
```
