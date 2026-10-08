# Pseudo-Huber Density, Standard-Deviation Parametrization

Computes the density \$\$f(y; \mu, \sigma, \nu) =
\dfrac{\sqrt{R(\nu)}}{2 \sigma \sqrt{\nu}\\ K_1(\sqrt{\nu})}
\exp\left(-\sqrt{\nu + R(\nu)
\left(\dfrac{y-\mu}{\sigma}\right)^2}\right), \qquad R(\nu) =
\dfrac{\sqrt{\nu}\\ K_2(\sqrt{\nu})}{K_1(\sqrt{\nu})},\$\$ whose
variance is \\\sigma^2\\. The log-density is \\-D - \log 2 -
\log\sigma + h(\nu)\\ with \\D = \sqrt{\nu + R(\nu)(y -
\mu)^2/\sigma^2}\\ and \\h(\nu) = -\tfrac14 \log\nu + \tfrac12 \log
K_2(\sqrt\nu) - \tfrac32 \log K_1(\sqrt\nu)\\. \\R\\ and \\h\\ are
formed with the exponentially scaled Bessel functions, so they stay
finite where \\K_1(\sqrt\nu)\\ underflows.

## Arguments

- distrib:

  A `PseudoHuber2Distrib` object, from
  [`pseudohuber2_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber2_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each of length 1
  or of the length of `y`.

- log:

  Logical of length 1. When `TRUE` the log-density is returned.

- ...:

  Unused.

## Value

A numeric vector of densities, one per observation.

## See also

[`pseudohuber2_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber2_distrib.md);
[`distrib_pdf.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.PseudoHuberDistrib.md)
for the scale parametrization.

## Examples

``` r
d <- pseudohuber2_distrib()
th <- list(mu = 0.4, sigma = 1.5, nu = 2)
integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value
#> [1] 1

# the same density as pseudohuber_distrib() at sigma1 = sigma / sqrt(R)
R <- sqrt(2) * besselK(sqrt(2), 2) / besselK(sqrt(2), 1)
y <- c(-2, 0, 3)
all.equal(distrib_pdf(d, y, th),
          distrib_pdf(pseudohuber_distrib(), y,
                      list(mu = 0.4, sigma = 1.5 / sqrt(R), nu = 2)))
#> [1] TRUE
```
