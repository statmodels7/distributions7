# Generalized Gamma Density

Computes the generalized gamma density in Stacy's form, \$\$f(y; a, d,
p) = \dfrac{p}{a^{d}\\\Gamma(d/p)}\\ y^{d-1} e^{-(y/a)^{p}}, \qquad y \>
0,\$\$ and 0 at or below zero. The compiled kernel works on the log
scale throughout, so a large \\d\\ does not form \\a^d\\ or \\y^{d-1}\\
separately.

The three parameters do three separate things: \\a\\ sets the scale,
\\d\\ the behavior near the origin (the density vanishes there for \\d
\> 1\\ and diverges for \\d \< 1\\), and \\p\\ the weight of the upper
tail.

## Arguments

- distrib:

  A `GenGamma1Distrib` object, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- y:

  A numeric vector of observations. A non-positive value gives a density
  of 0.

- theta:

  A named list with components `a`, `d` and `p`, each a numeric vector
  of length 1 or of the length of `y`. All three must be strictly
  positive.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`. The result does not depend on the count.

## Value

A numeric vector of densities, of the length of the recycled inputs.

## Notation

\\a \> 0\\ is the scale, \\d \> 0\\ and \\p \> 0\\ the two shapes, and
\\\Gamma\\ the gamma function. \\a\\ is not the mean, which is
\\a\\\Gamma\\(d+1)/p\\/\Gamma(d/p)\\.

## See also

[`distrib_cdf.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GenGamma1Distrib.md)
for the distribution function,
[`distrib_gradient.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GenGamma1Distrib.md)
for the score, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- gengamma1_distrib()
y <- c(0.5, 1.5, 4)
th <- list(a = 2, d = 3, p = 1.5)

# The formula written out.
all.equal(distrib_pdf(d, y, th),
          1.5 / (2^3 * gamma(3 / 1.5)) * y^(3 - 1) * exp(-(y / 2)^1.5))
#> [1] TRUE

# It integrates to one.
integrate(function(v) distrib_pdf(d, v, th), 0, Inf)$value
#> [1] 1

# The four exact special cases.
c(gamma = max(abs(distrib_pdf(d, y, list(a = 2, d = 3, p = 1)) -
                  dgamma(y, shape = 3, scale = 2))),
  weibull = max(abs(distrib_pdf(d, y, list(a = 2, d = 1.5, p = 1.5)) -
                    dweibull(y, shape = 1.5, scale = 2))),
  exponential = max(abs(distrib_pdf(d, y, list(a = 2, d = 1, p = 1)) -
                        dexp(y, rate = 1 / 2))),
  half_normal = max(abs(distrib_pdf(d, y, list(a = sqrt(2), d = 1, p = 2)) -
                        2 * dnorm(y))))
#>        gamma      weibull  exponential  half_normal 
#> 1.387779e-17 5.551115e-17 1.387779e-17 1.110223e-16 

# d decides what happens at the origin.
distrib_pdf(d, 1e-8, list(a = 2, d = c(0.5, 1, 2), p = 1.5))
#> [1] 3959.255
```
