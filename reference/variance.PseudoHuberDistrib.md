# Variance of the Pseudo-Huber Distribution

Closed form, replacing the numerical default: \$\$\operatorname{Var}(Y)
= \sigma^2 \sqrt{\nu}\\ \frac{K_2(\sqrt{\nu})}{K_1(\sqrt{\nu})},\$\$
with \\K_r\\ the modified Bessel function of the second kind. The scale
parameter is therefore not the standard deviation: the shape multiplies
it, and the factor grows roughly like \\\sqrt{\nu}\\ at large \\\nu\\.

## Arguments

- x:

  A `PseudoHuberDistrib`, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- theta:

  A named list with components `mu` (the location), `sigma` (the scale,
  positive) and `nu` (the shape, positive), each a numeric vector of
  length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of the length the recycled parameters
imply. The scale and `nu` enter the value, so a setting that varies `mu`
alone repeats one number.

## Details

The two Bessel functions are evaluated exponentially scaled, through
`besselK(sqrt(nu), r, expon.scaled = TRUE)`. Both carry the same factor
\\e^{\sqrt{\nu}}\\, which cancels in the ratio, so the quotient is exact
where the unscaled functions would have overflowed; the terms are degree
homogeneous, and the scaling has been checked to \\\nu = 2000\\.

The quantity is also used inside the package: the bracket
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
searches over for this family is built from it, so an error here would
surface as a failed root search several frames away.

## Notation

\\\sigma \> 0\\ is the scale, \\\nu \> 0\\ the shape, and \\K_r\\ the
modified Bessel function of the second kind of order \\r\\.

## See also

[`kurtosis.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.PseudoHuberDistrib.md),
which is a ratio of the same Bessel functions;
[`mean.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.PseudoHuberDistrib.md);
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

## Examples

``` r
d <- pseudohuber_distrib()

# The Bessel ratio, written out against the method.
nu <- 4
all.equal(variance(d, list(mu = 0, sigma = 1, nu = nu)),
          sqrt(nu) * besselK(sqrt(nu), 2, expon.scaled = TRUE) /
                     besselK(sqrt(nu), 1, expon.scaled = TRUE))
#> [1] TRUE

# sigma is not the standard deviation: the shape scales it.
variance(d, list(mu = 0, sigma = 1, nu = c(1, 4, 100)))
#> [1]  2.699484  3.628616 11.534173
```
