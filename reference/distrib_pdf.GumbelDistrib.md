# Gumbel Probability Density Function

Computes the Gumbel density, with \\z = (y - \mu)/\sigma\\, \$\$f(y;
\mu, \sigma) = \dfrac{1}{\sigma} \exp\left\\-z - e^{-z}\right\\,\$\$
written out rather than delegated, base R carrying no Gumbel. The two
tails are of very different weight: the density falls like \\e^{-z}\\ to
the right and like \\e^{-e^{-z}}\\ to the left, so the left one is
doubly exponential and dies far faster.

With `log = TRUE` the exponent is returned directly and stays finite
where the density itself underflows.

## Arguments

- distrib:

  A `GumbelDistrib` object, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

- y:

  A numeric vector of observations. Every real value is in the support,
  so no value is rejected.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma))`, one value per observation.

## See also

[`distrib_cdf.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GumbelDistrib.md)
for the distribution function,
[`distrib_gradient.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GumbelDistrib.md)
for the derivatives of the log-density,
[`distrib_pdf.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Weibull1Distrib.md)
for the family this becomes under \\e^{-Y}\\, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- gumbel_distrib()
y <- c(-1, 0, 1)
th <- list(mu = 0, sigma = 1)

# The closed form, written out.
z <- (y - 0) / 1
all.equal(distrib_pdf(d, y, th), exp(-z - exp(-z)))
#> [1] TRUE

# The mode is at mu, where the density peaks at 1/(sigma e).
c(at_mode = distrib_pdf(d, 0, th), one_over_e = 1 / exp(1))
#>    at_mode one_over_e 
#>  0.3678794  0.3678794 

# The two tails are of very different weight.
distrib_pdf(d, c(-5, 5), th)
#> [1] 5.205427e-63 6.692700e-03

# Far to the left the density underflows and its logarithm does not.
distrib_pdf(d, -6, th)
#> [1] 2.505348e-173
distrib_pdf(d, -6, th, log = TRUE)
#> [1] -397.4288
```
