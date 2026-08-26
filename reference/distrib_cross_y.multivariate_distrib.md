# No Numerical Mixed Response Derivative in Several Dimensions

Signals an error. Unlike the two refusals beside it, the obstruction
here is not that a line differences wrongly: it is that no consumer has
fixed the shape. The univariate generic returns one number per
observation per parameter; the multivariate quantity is one \\n \times
p\\ matrix per parameter, and a family that supplies it registers a
method returning that, as the gaussian and the Student t do.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object with no method of its own.

- y:

  An \\n \times p\\ numeric matrix of observations. Not examined.

- theta:

  A named list of parameters. Not examined.

- scale:

  One of `"parameter"` or `"link"`, handled by the generic before
  dispatch. Not examined.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Never returns: it always signals an error naming the family and the
shape a method should produce.

## See also

[`distrib_cross_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.MvGaussianDistrib.md)
and
[`distrib_cross_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.MvStudentTDistrib.md)
for the two families that answer, and
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.

## Examples

``` r
d <- dirichlet_distrib(3)
theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
set.seed(1)
y <- distrib_rng(d, 3, theta)
try(distrib_cross_y(d, y, theta))
#> Error : distrib_cross_y() is not defined for 'dirichlet [3d, mean=simplex]': register a closed form on the family, returning one n-by-p matrix per parameter as the gaussian and the Student t do.

# The shape a method returns: one n by p matrix per parameter.
g <- mvgaussian_distrib(2)
th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
           sigma_L2.1 = 0.5)
cy <- distrib_cross_y(g, rbind(c(1, -1), c(0, 0)), th)
c(entries = length(cy), rows = nrow(cy[[1]]), cols = ncol(cy[[1]]))
#> entries    rows    cols 
#>       5       2       2 
```
