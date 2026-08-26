# Elastic-Net Density

Computes the elastic-net density \$\$f(y; \mu, \lambda, \alpha) =
\frac{1}{Z} \exp\left\\-a\|y-\mu\| - \tfrac{c}{2}(y-\mu)^2\right\\,\$\$
with \\a = \lambda\alpha\\, \\c = \lambda(1-\alpha)\\ and \\Z =
2M(a/\sqrt{c})/\sqrt{c}\\, where \\M\\ is the Mills ratio. The constant
is finite at both ends of the mixing weight: \\2/a\\ as \\\alpha \to 1\\
and \\\sqrt{2\pi/c}\\ as \\\alpha \to 0\\.

The constant goes through \\\log M\\. Written directly, \\\log M(x)\\
adds \\x^2/2\\ to a log-probability of the same size and opposite sign,
so it loses a digit for every factor of ten in \\x\\; past \\x = 30\\
the asymptotic series is used instead. At \\\alpha = 1 - 10^{-12}\\ the
argument reaches \\10^6\\, which an ordinary elastic net does.

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- y:

  A numeric vector of observations, anywhere on the real line.

- theta:

  A named list with components `mu`, `lambda` and `alpha`, each a
  numeric vector of length 1 or of the length of `y`. `lambda` must be
  strictly positive and `alpha` strictly inside \\(0, 1)\\; at either
  endpoint the constant's argument is not defined and the value is
  `NaN`.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of the length of the recycled inputs.

## Notation

\\\mu\\ is the location, \\\lambda \> 0\\ the overall rate, \\\alpha \in
(0,1)\\ the mixing weight, \\M\\ the Mills ratio \\\Phi(-x)/\phi(x)\\,
and \\Z\\ the normalizing constant.

## See also

[`distrib_cdf.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.EnetDistrib.md)
for the distribution function,
[`distrib_gradient.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.EnetDistrib.md)
for the score,
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
and
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the two ends, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- enet_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, lambda = 2, alpha = 0.5)

# It integrates to one.
integrate(function(u) distrib_pdf(d, u, th), -Inf, Inf)$value
#> [1] 1

# The kernel written out, up to the constant.
r <- distrib_pdf(d, y, th) /
     exp(-1 * abs(y) - 1 * y^2 / 2)     # a = c = 1 at alpha = 0.5
all.equal(r, rep(r[1], length(r)))
#> [1] TRUE

# The two ends, approached to ten figures.
rbind(alpha = c(1 - 1e-10, 1e-10),
      enet = c(distrib_pdf(d, 0.7, list(mu = 0, lambda = 2,
                                        alpha = 1 - 1e-10)),
               distrib_pdf(d, 0.7, list(mu = 0, lambda = 2,
                                        alpha = 1e-10))),
      limit = c(distrib_pdf(laplace2_distrib(), 0.7,
                            list(mu = 0, lambda = 2)),
                dnorm(0.7, 0, 1 / sqrt(2))))
#>           [,1]         [,2]
#> alpha 1.000000 0.0000000001
#> enet  0.246597 0.3456374302
#> limit 0.246597 0.3456374302
```
