# Exponential First Derivative in the Response

Returns \\-1/\mu\\ for every observation. The log-density is
\\-\log\mu - y/\mu\\, which is **linear** in \\y\\, so its derivative in
the response is the constant slope \\-1/\mu\\ and carries no information
about where the observation fell.

That linearity is the memorylessness of the family stated on the log
scale: the log survival function is \\-q/\mu\\, a straight line, so the
hazard is constant.

## Arguments

- distrib:

  An `ExponentialDistrib` object, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. `mu` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `length(y)`, every entry \\-1/\mu\\.

## See also

[`distrib_hess_y.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.ExponentialDistrib.md)
for the second derivative, which is zero;
[`distrib_gradient.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ExponentialDistrib.md)
for the score in the mean;
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- exponential_distrib()
y <- c(0.3, 1.1, 4.0)
th <- list(mu = 2)

distrib_grad_y(d, y, th)
#> [1] -0.5 -0.5 -0.5
-1 / 2
#> [1] -0.5

# It is the derivative of the log-density, so a central difference in y
# reproduces it.
eps <- 1e-6
(distrib_pdf(d, y + eps, th, log = TRUE) -
  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
#> [1] -0.5 -0.5 -0.5

# The same slope appears in the log survival function: a constant hazard.
diff(distrib_cdf(d, c(1, 2, 3), th, lower.tail = FALSE, log.p = TRUE))
#> [1] -0.5 -0.5
```
