# Elastic-Net Response Derivative

Computes \\\partial\ell/\partial y = -a\\\mathrm{sgn}(y-\mu) -
c(y-\mu)\\, with \\a = \lambda\alpha\\ and \\c = \lambda(1-\alpha)\\. It
is minus the derivative in \\\mu\\, the response and the location
entering only through their difference.

It is **undefined at the location**, where the sign jumps. `sign(0)` is
0 in R, so the value returned there is \\-c(y-\mu) = 0\\, one point of
the subdifferential.

The score is bounded below by \\-a - c(y-\mu)\\ on either side, so the
Laplace part contributes a jump of \\2a\\ across the location and the
Gaussian part the slope.

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `lambda` and `alpha`, each a
  numeric vector of length 1 or of the length of `y`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the length of the recycled inputs.

## Notation

\\\ell\\ is the log-density of one observation, \\a = \lambda\alpha\\
and \\c = \lambda(1-\alpha)\\.

## See also

[`distrib_hess_y.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.EnetDistrib.md)
for the second derivative,
[`distrib_cross_y.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.EnetDistrib.md)
for the mixed one,
[`distrib_gradient.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.EnetDistrib.md)
for the parameter derivatives, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- enet_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, lambda = 2, alpha = 0.5)

# It is minus the score in the location, exactly.
all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#> [1] TRUE

# Against a central difference of the log-density, away from the kink.
eps <- 1e-6
rbind(analytic = distrib_grad_y(d, y, th),
      numeric = (distrib_pdf(d, y + eps, th, log = TRUE) -
                 distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps))
#>          [,1] [,2] [,3] [,4]
#> analytic  2.5  1.3 -1.4 -3.1
#> numeric   2.5  1.3 -1.4 -3.1

# The jump across the location is 2a, which is the Laplace part.
c(just_below = distrib_grad_y(d, -1e-9, th),
  just_above = distrib_grad_y(d, 1e-9, th),
  two_a = 2 * 2 * 0.5)
#> just_below just_above      two_a 
#>          1         -1          2 
```
