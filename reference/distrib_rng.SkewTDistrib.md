# Skew t Random Generation

Draws from the skew \\t\\ exactly, from its scale-mixture
representation: with \\Z\\ standard skew normal of shape \\\alpha\\ and
\\V \sim \chi^2\_\nu\\ independent of it, \$\$Y = \mu +
\sigma\\\dfrac{Z}{\sqrt{V/\nu}}.\$\$ \\Z\\ itself is drawn from
[`distrib_rng.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.SkewNormal1Distrib.md)'s
representation, so the whole draw is three `rnorm`/`rchisq` calls and no
inversion or rejection.

The representation reads off both limits. As \\\nu \to \infty\\ the
mixing factor tends to one and \\Y\\ is skew normal; at \\\alpha = 0\\
the numerator is Gaussian and \\Y\\ is Student \\t\\.

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`, each a
  numeric vector of length 1 or of length `n`; a component of length 1
  is recycled.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws.

## Notation

\\\nu\\ is the degrees of freedom of the mixing chi-squared and \\\delta
= \alpha/\sqrt{1+\alpha^2}\\ the weight of the half-normal component of
\\Z\\.

## See also

[`distrib_pdf.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.SkewTDistrib.md)
for the density the draws follow,
[`distrib_rng.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.SkewNormal1Distrib.md)
for the inner representation, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- skewt_distrib()
th <- list(mu = 1, sigma = 2, alpha = 3, nu = 8)

set.seed(21)
x <- distrib_rng(d, 1e5, th)

# Two moments against their closed forms. Both exist here; at nu <= 2 the
# variance does not.
rbind(sample = c(mean(x), var(x)),
      theory = c(mean(d, th), variance(d, th)))
#>            [,1]     [,2]
#> sample 2.677048 2.534376
#> theory 2.677051 2.520833

# The mixture, written out at the same seed.
set.seed(21)
delta <- 3 / sqrt(1 + 9)
z <- delta * abs(rnorm(1e5)) + sqrt(1 - delta^2) * rnorm(1e5)
all.equal(x, 1 + 2 * z / sqrt(rchisq(1e5, df = 8) / 8))
#> [1] TRUE
```
