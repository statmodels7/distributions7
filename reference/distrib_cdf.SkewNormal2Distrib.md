# Skew Normal Distribution Function in the Centered Parametrization

Computes the skew normal distribution function at the direct parameters
the centered triple implies, through Azzalini's Owen's T identity \\F(q)
= \Phi(z) - 2T(z, \alpha)\\ with \\z = (q-\xi)/\omega\\. The arithmetic
is
[`distrib_cdf.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.SkewNormal1Distrib.md)'s;
this method supplies the translated parameters.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object, from
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` and `gamma1`, each a
  numeric vector of length 1 or of the length of `q`.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, the value is \\P(Y \le
  q)\\; when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned. Defaults to `FALSE`.

- ...:

  Passed to
  [`distrib_cdf.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.SkewNormal1Distrib.md).

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
with `log.p = TRUE`.

## Notation

\\\gamma_1\\ is the skewness, \\(\xi, \omega, \alpha)\\ the implied
location, scale and shape, and \\T\\ Owen's T function.

## See also

[`distrib_cdf.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.SkewNormal1Distrib.md)
for the identity,
[`distrib_quantile.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.SkewNormal2Distrib.md)
for its inverse, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- skewnormal2_distrib()
q <- c(-2, -0.5, 0.5, 2)
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)

# Against a direct quadrature of the density.
rbind(owen = distrib_cdf(d, q, th),
      quadrature = vapply(q, function(u)
        integrate(function(v) distrib_pdf(d, v, th), -Inf, u)$value, 0))
#>                   [,1]      [,2]      [,3]      [,4]
#> owen       0.009070287 0.3311782 0.7156719 0.9645037
#> quadrature 0.009070287 0.3311782 0.7156719 0.9645037

# At zero skewness it is the Gaussian's.
all.equal(distrib_cdf(d, q, list(mu = 0, sigma = 1, gamma1 = 0)), pnorm(q))
#> [1] TRUE

# A positive skewness puts more than half the mass below the mean.
distrib_cdf(d, 0, th)
#> [1] 0.5374817
```
