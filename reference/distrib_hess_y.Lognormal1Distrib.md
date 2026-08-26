# Lognormal Second Derivative in the Response

Computes the second derivative of the lognormal log-density with respect
to the response, in closed form. With \\r = \log y - \mu\\,
\$\$\dfrac{\partial^2 \ell}{\partial y^2} = \dfrac{1}{y^2}\left(1 +
\dfrac{r - 1}{\sigma^2}\right).\$\$ The sign changes at \\r = 1 -
\sigma^2\\, that is at \\y = e^{\mu + 1 - \sigma^2}\\: the log-density
is **concave** in the response below that point and **convex** above it,
so a lognormal log-likelihood is not a concave function of an
observation over the whole support.

## Arguments

- distrib:

  A `Lognormal1Distrib` object, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- y:

  A numeric vector of strictly positive observations. The value diverges
  as `y` approaches zero.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma2` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu), length(sigma2))`,
one value per observation.

## See also

[`distrib_grad_y.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Lognormal1Distrib.md)
for the first derivative in the response,
[`distrib_hessian.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Lognormal1Distrib.md)
for the curvature in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- lognormal1_distrib()
y <- c(0.5, 1.6, 4)
th <- list(mu = 0.5, sigma2 = 0.36)

r <- log(y) - 0.5
all.equal(distrib_hess_y(d, y, th), (1 + (r - 1) / 0.36) / y^2)
#> [1] TRUE

# Concave below exp(mu + 1 - sigma2) and convex above it.
cut <- exp(0.5 + 1 - 0.36)
c(cut = cut,
  below = distrib_hess_y(d, cut / 2, th),
  at = distrib_hess_y(d, cut, th),
  above = distrib_hess_y(d, 2 * cut, th))
#>           cut         below            at         above 
#>  3.126768e+00 -7.877557e-01  3.406748e-17  4.923473e-02 
```
