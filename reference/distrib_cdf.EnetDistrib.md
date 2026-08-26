# Elastic-Net Distribution Function

Computes the distribution function in closed form. Each half of the
density is a truncated Gaussian, so with \\z = q-\mu\\ and \\x = a/\sqrt
c\\, \$\$F(q) = \dfrac{\Phi(\sqrt{c}\\z - x)}{2\Phi(-x)} \quad (z \le
0), \qquad F(q) = 1 - \dfrac{\Phi(-\sqrt{c}\\z - x)}{2\Phi(-x)} \quad (z
\> 0).\$\$

Both ratios are taken on the log scale and exponentiated. The
denominator \\\Phi(-x)\\ is exactly zero past \\x = 38\\ in double
precision, and \\x\\ reaches that at ordinary values of \\\alpha\\: at
\\\lambda = 20\\, \\\alpha = 0.995\\ it is 63.

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `lambda` and `alpha`, each a
  numeric vector of length 1 or of the length of `q`. `lambda` must be
  strictly positive and `alpha` strictly inside \\(0, 1)\\.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, the value is \\P(Y \le
  q)\\; when `FALSE` it is one minus that.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
with `log.p = TRUE`, of the length of the recycled inputs.

## Notation

\\a = \lambda\alpha\\, \\c = \lambda(1-\alpha)\\, \\x = a/\sqrt c\\, and
\\\Phi\\ the standard Gaussian distribution function.

## See also

[`distrib_quantile.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.EnetDistrib.md),
which inverts this in closed form,
[`distrib_pdf.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.EnetDistrib.md)
for the density, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- enet_distrib()
q <- c(-2, -0.5, 0, 0.5, 2)
th <- list(mu = 0, lambda = 2, alpha = 0.5)

# Against a direct quadrature of the density.
rbind(closed = distrib_cdf(d, q, th),
      quadrature = vapply(q, function(u)
        integrate(function(v) distrib_pdf(d, v, th), -Inf, u)$value, 0))
#>                   [,1]     [,2] [,3]     [,4]      [,5]
#> closed     0.004254186 0.210542  0.5 0.789458 0.9957458
#> quadrature 0.004254189 0.210542  0.5 0.789458 0.9957458

# The density is symmetric about the location, so the median is mu.
distrib_cdf(d, 0, th)
#> [1] 0.5

# An ordinary elastic net puts the Mills argument past 38, where the
# denominator of both halves has underflowed to zero on the natural scale.
hard <- list(mu = 0, lambda = 20, alpha = 0.995)
c(x = 20 * 0.995 / sqrt(20 * 0.005),
  phi_minus_x = pnorm(-20 * 0.995 / sqrt(20 * 0.005)),
  cdf = distrib_cdf(d, 0.1, hard))
#>           x phi_minus_x         cdf 
#>  62.9293254   0.0000000   0.9317207 
```
