# Skew Normal Density

Computes the skew normal density, with \\z = (y-\mu)/\sigma\\: \$\$f(y;
\mu, \sigma, \alpha) = \dfrac{2}{\sigma}\\\phi(z)\\\Phi(\alpha z).\$\$
The Gaussian density is multiplied by \\2\Phi(\alpha z)\\, which exceeds
one where \\\alpha z \> 0\\ and falls below it where \\\alpha z \< 0\\,
so a positive \\\alpha\\ moves mass to the right. At \\\alpha = 0\\ the
factor is one everywhere and the density is \\\phi(z)/\sigma\\.

The two logarithms are taken separately and added, `dnorm(log = TRUE)`
beside `pnorm(log.p = TRUE)`, so the light tail of the skewed side
returns a large negative number where forming \\\Phi(\alpha z)\\ first
and then taking its logarithm would return `-Inf`.

## Arguments

- distrib:

  A `SkewNormal1Distrib` object, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- y:

  A numeric vector of observations, anywhere on the real line.

- theta:

  A named list with components `mu`, `sigma` and `alpha`, each a numeric
  vector of length 1 or of the length of `y`; a component of length 1 is
  recycled. `sigma` must be strictly positive, and `mu` and `alpha` may
  take any finite value.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma), length(alpha))`, one value
per observation.

## Notation

\\\phi\\ and \\\Phi\\ are the standard Gaussian density and distribution
function, \\\mu\\ the location, \\\sigma \> 0\\ the scale and \\\alpha\\
the shape. \\\mu\\ is not the mean unless \\\alpha = 0\\.

## See also

[`distrib_cdf.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.SkewNormal1Distrib.md)
for Owen's T identity,
[`distrib_gradient.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal1Distrib.md)
for the score,
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
for the family, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- skewnormal1_distrib()
y <- c(-2, -0.5, 0.5, 2)
th <- list(mu = 0, sigma = 1, alpha = 3)

# The formula written out.
all.equal(distrib_pdf(d, y, th),
          2 * dnorm(y) * pnorm(3 * y))
#> [1] TRUE

# It integrates to one.
integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value
#> [1] 1

# At shape zero the tilting factor is one and the family is the Gaussian.
all.equal(distrib_pdf(d, y, list(mu = 0, sigma = 1, alpha = 0)), dnorm(y))
#> [1] TRUE

# The light tail stays a number on the log scale where the density itself
# has underflowed to zero.
rbind(density = distrib_pdf(d, -40, th),
      log_density = distrib_pdf(d, -40, th, log = TRUE))
#>                  [,1]
#> density         0.000
#> log_density -8005.932
```
