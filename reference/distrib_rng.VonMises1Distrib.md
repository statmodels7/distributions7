# von Mises Random Generation

Draws `n` independent angles by the rejection algorithm of Best and
Fisher (1979), which proposes from a wrapped Cauchy envelope and accepts
with a probability that involves no Bessel function, so that the
normalizing constant, the costly part of the density, is never
evaluated. The constant \\b\\ of the envelope is computed as
\\2a\kappa/\\(s+1)(a+\sqrt{2a})\\\\ with \\s = \sqrt{1+4\kappa^2}\\ and
\\a = 1 + s\\, which keeps its digits at a small concentration; below
\\\kappa = 10^{-300}\\ the envelope overflows, the law is uniform to
double precision, and the draws are uniform on \\\[-\pi, \pi)\\.

The accepted angles are drawn about zero and then shifted by \\\mu\\ and
wrapped back into \\\[-\pi, \pi)\\, so every draw lies in the declared
support. With a scalar `kappa` the envelope's constants are built once
and the loop proposes more candidates than needed and repeats until `n`
draws have been accepted. With a `kappa` that varies by observation each
draw has its own envelope, and the loop proposes once for every draw not
yet accepted, so that draw `i` is generated at `kappa[i]` and shifted by
`mu[i]`. Either way the number of uniform variates consumed depends on
the draws.

## Arguments

- distrib:

  A `VonMises1Distrib` object, from
  [`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `kappa`, each a numeric vector
  of length 1 or `n`. As in the generic, a component of length 1 is
  recycled to `n`. `mu` must lie in \\(-\pi, \pi)\\ and `kappa` be
  strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` angles in \\\[-\pi, \pi)\\.

## References

Best, D. J. and Fisher, N. I. (1979). Efficient simulation of the von
Mises distribution. *Journal of the Royal Statistical Society, Series
C*, **28**(2), 152-157.

## See also

[`distrib_pdf.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.VonMises1Distrib.md)
for the density,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- vonmises1_distrib()
th <- list(mu = 0.5, kappa = 2)
set.seed(1)
z <- distrib_rng(d, 3e5, th)

# Every draw is in the declared support.
range(z)
#> [1] -3.141443  3.140308

# The circular mean recovers mu and the mean resultant length recovers
# A(kappa) = I_1 / I_0; the ordinary mean recovers neither.
c(circular_mean = atan2(mean(sin(z)), mean(cos(z))), mu = 0.5)
#> circular_mean            mu 
#>     0.4990003     0.5000000 
c(resultant = sqrt(mean(cos(z))^2 + mean(sin(z))^2),
  A = numericals7::bessel_i_ratio(2))
#> resultant         A 
#> 0.6982004 0.6977747 
```
