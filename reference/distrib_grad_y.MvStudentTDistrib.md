# Multivariate Student t Response Gradient

Computes the derivative of the log-density in the response,
\$\$\frac{\partial \ell}{\partial y} = -c\\\Sigma^{-1}(y-\mu), \qquad c
= \frac{\nu+p}{\nu+q},\$\$ one row per observation: the gaussian's
expression with the family's weight in front of it. The weight is what
bounds the influence of a distant row. A gaussian's response gradient
grows without limit as the observation moves away; this one rises, turns
and decays like \\1/\lVert y\rVert\\, which is the redescending score of
a heavy-tailed family.

It is also minus the score in the location, this being a location
family, so the same numbers appear in
[`distrib_gradient.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvStudentTDistrib.md).

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

An \\n \times p\\ numeric matrix, row \\i\\ holding
\\\partial\ell_i/\partial y_i\\.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom, \\p\\ the dimension, \\q\\ the squared Mahalanobis
distance and \\c\\ the weight.

## See also

[`distrib_hess_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvStudentTDistrib.md)
for the second derivative,
[`distrib_cross_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.MvStudentTDistrib.md)
for the mixed block,
[`mvt_weights()`](https://statmodels7.github.io/distributions7/reference/mvt_weights.md)
for the weight, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- mvstudent_t_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(1)
y <- distrib_rng(d, 4, theta)

distrib_grad_y(d, y, theta)
#>           [,1]       [,2]
#> [1,]  0.821471 -0.4640124
#> [2,] -0.728796  1.3809482
#> [3,]  1.422240 -0.8714900
#> [4,] -1.037485 -0.8375372

# Against a numerical derivative taken row by row.
num <- t(apply(y, 1, function(r)
  numDeriv::grad(function(z) distrib_pdf(d, z, theta, log = TRUE), r)))
max(abs(distrib_grad_y(d, y, theta) - num))
#> [1] 8.781598e-11

# It is minus the score in the location.
g <- distrib_gradient(d, y, theta)
all.equal(distrib_grad_y(d, y, theta), -cbind(g$mu1, g$mu2),
          check.attributes = FALSE)
#> [1] TRUE

# It redescends: the gaussian's grows with the distance and the t's does
# not, along the first coordinate from the location.
far <- cbind(0.5 + c(1, 3, 10, 40), -0.3)
g0 <- mvgaussian_distrib(2)
rbind(t = distrib_grad_y(d, far, theta)[, 1],
      gaussian = distrib_grad_y(g0, far, theta[1:5])[, 1])
#>               [,1]      [,2]        [,3]        [,4]
#> t        -1.156695 -1.608983  -0.7553137  -0.1992632
#> gaussian -1.014155 -3.042466 -10.1415519 -40.5662078
```
