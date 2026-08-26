# Pseudo-Huber First Derivative in the Response

Computes \\\partial \ell / \partial y\\ in closed form. With \\r = y -
\mu\\ and \\D = \sqrt{\nu + (r/\sigma)^2}\\, \$\$\dfrac{\partial
\ell}{\partial y} = -\dfrac{r}{\sigma^2 D}.\$\$ The family is a location
family in \\\mu\\, so this is exactly the negative of the location score
[`distrib_gradient.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuberDistrib.md)`$mu`,
and like it, it is **bounded**: it tends to \\\mp 1/\sigma\\ as the
residual runs away. This quantity is what a quantile residual's
delta-method standard error and a change of variable in the response
both need.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length
`max(length(y), length(mu), length(sigma), length(nu))`, one value per
observation.

## See also

[`distrib_hess_y.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.PseudoHuberDistrib.md)
for the second derivative in the response,
[`distrib_cross_y.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.PseudoHuberDistrib.md)
for the mixed derivative,
[`distrib_gradient.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuberDistrib.md)
for the derivatives in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 2)

# The closed form, written out.
r <- y - 0.4; D <- sqrt(2 + (r / 1.2)^2)
all.equal(distrib_grad_y(d, y, th), -r / (1.2^2 * D))
#> [1] TRUE

# A location family, so this is minus the location score.
all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#> [1] TRUE

# A central difference of the log-density in y reproduces it.
eps <- 1e-6
all.equal((distrib_pdf(d, y + eps, th, log = TRUE) -
           distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps),
          distrib_grad_y(d, y, th), tolerance = 1e-6)
#> [1] TRUE

# Bounded by 1 / sigma however far out the response is.
c(distrib_grad_y(d, 0.4 + c(10, 1000, 1e6), th), bound = -1 / 1.2)
#>                                       bound 
#> -0.8215865 -0.8333321 -0.8333333 -0.8333333 
```
