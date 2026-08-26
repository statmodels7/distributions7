# Lognormal First Derivative in the Response

Computes the first derivative of the lognormal log-density with respect
to the response, in closed form. With \\r = \log y - \mu\\,
\$\$\dfrac{\partial \ell}{\partial y} = -\dfrac{1}{y}\left(1 +
\dfrac{r}{\sigma^2}\right).\$\$ This is where the family parts company
with the Gaussian: the Jacobian \\1/y\\ of the log transformation
carries no parameter and so leaves every derivative in \\\mu\\ and
\\\sigma^2\\ alone, but it is a function of the response and does enter
here.

The derivative vanishes at \\r = -\sigma^2\\, that is at \\y = e^{\mu -
\sigma^2}\\, which is the mode of the density and lies below both the
median \\e^{\mu}\\ and the mean \\e^{\mu + \sigma^2/2}\\.

Quantile residuals and the mixed derivatives of
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
read it.

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

[`distrib_hess_y.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Lognormal1Distrib.md)
for the second derivative in the response,
[`distrib_gradient.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal1Distrib.md)
for the score in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- lognormal1_distrib()
y <- c(0.5, 1.6, 4)
th <- list(mu = 0.5, sigma2 = 0.36)

# Written out.
r <- log(y) - 0.5
all.equal(distrib_grad_y(d, y, th), -(1 + r / 0.36) / y)
#> [1] TRUE

# Zero at the mode exp(mu - sigma2), which lies below the median exp(mu).
mode <- exp(0.5 - 0.36)
c(mode = mode, median = exp(0.5), at_mode = distrib_grad_y(d, mode, th))
#>          mode        median       at_mode 
#>  1.150274e+00  1.648721e+00 -9.651815e-17 

# It is the derivative of the log-density, so a central difference in y
# reproduces it.
eps <- 1e-7
(distrib_pdf(d, y + eps, th, log = TRUE) -
  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
#> [1]  4.6285954 -0.5729230 -0.8654822
```
