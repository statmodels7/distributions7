# von Mises Cumulative Distribution Function

Computes \\F(q) = P(Y \le q)\\ on \\\[-\pi, \pi)\\ from the Fourier
series of the density integrated term by term, through
[`vm_cdf()`](https://statmodels7.github.io/distributions7/reference/vm_cdf.md).
The density has no elementary antiderivative, and the quadrature the
base class would use costs one integration per observation; the series
replaces \\n\\ quadratures with a few dozen vectorized steps, using only
the Bessel **ratios** \\I_j/I_0\\ from a backward recurrence.

Below \\-\pi\\ the value is 0 and at or above \\\pi\\ it is 1, the
support being the declared interval. The result is clamped to \\\[0,
1\]\\: the series is exact and its rounding is not.

## Arguments

- distrib:

  A `VonMises1Distrib` object, from
  [`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md).

- q:

  A numeric vector of angles.

- theta:

  A named list with components `mu` and `kappa`, each a numeric vector
  of length 1 or of the length of `q`. A component of length 1 is
  recycled. `mu` must lie in \\(-\pi, \pi)\\ and `kappa` be strictly
  positive; a non-positive `kappa` gives `NA`.

- ...:

  Unused, and accepted so that the signature matches the generic's. This
  method takes **no** `lower.tail` or `log.p`: the upper tail is
  `1 - F(q)` and the logarithm is `log(F(q))`.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(q), length(mu), length(kappa))`.

## See also

[`vm_cdf()`](https://statmodels7.github.io/distributions7/reference/vm_cdf.md)
for the series and its measured term count,
[`distrib_pdf.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.VonMises1Distrib.md)
for the density,
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md),
which inverts this by root finding, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- vonmises1_distrib()
th <- list(mu = 0.5, kappa = 2)
y <- c(-1, 0, 0.5, 2)

# The series agrees with a direct quadrature of the density.
rbind(series = distrib_cdf(d, y, th),
      quadrature = vapply(y, function(v)
        integrate(function(u) distrib_pdf(d, u, th), -pi, v)$value,
        numeric(1)))
#>                  [,1]      [,2]      [,3]      [,4]
#> series     0.04797676 0.2669514 0.5051436 0.9623105
#> quadrature 0.04797676 0.2669514 0.5051436 0.9623105

# It runs from 0 to 1 across the declared support.
c(distrib_cdf(d, -pi, th), distrib_cdf(d, pi - 1e-12, th))
#> [1] 0 1

# And the quantile function, which the base class obtains by root finding
# on this, inverts it.
q <- distrib_quantile(d, c(0.25, 0.5, 0.75), th)
all.equal(distrib_cdf(d, q, th), c(0.25, 0.5, 0.75), tolerance = 1e-6)
#> [1] TRUE
```
