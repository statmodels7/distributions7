# von Mises Distribution Function in the Resultant Length

Computes \\F(q) = P(Y \le q)\\ on \\\[-\pi, \pi)\\ from the Fourier
series of the concentration parametrization, read at \\\kappa =
A^{-1}(\rho)\\. The map touches the second parameter only and the
response not at all, so the distribution function is the other family's
at that concentration.

What it replaces is the base class's quadrature, one integration per
observation. See
[`vm_cdf()`](https://statmodels7.github.io/distributions7/reference/vm_cdf.md)
for the series, its measured term count and the blocking that keeps the
intermediate small.

## Arguments

- distrib:

  A `VonMises2Distrib` object, from
  [`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md).

- q:

  A numeric vector of angles. Below \\-\pi\\ the value is 0 and at or
  above \\\pi\\ it is 1.

- theta:

  A named list with components `mu` and `rho`, each a numeric vector of
  length 1 or of the length of `q`. `mu` must lie in \\(-\pi, \pi)\\ and
  `rho` in \\(0, 1)\\.

- ...:

  Unused, and accepted so that the signature matches the generic's. This
  method takes **no** `lower.tail` or `log.p`: the upper tail is
  `1 - F(q)` and the logarithm is `log(F(q))`.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(q), length(mu), length(rho))`.

## See also

[`vm_cdf()`](https://statmodels7.github.io/distributions7/reference/vm_cdf.md)
for the series,
[`distrib_cdf.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.VonMises1Distrib.md)
for the same quantity in the concentration, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d2 <- vonmises2_distrib()
th <- list(mu = 0.5, rho = 0.7)
y <- c(-1, 0, 0.5, 2)

# The series agrees with a direct quadrature of the density.
rbind(series = distrib_cdf(d2, y, th),
      quadrature = vapply(y, function(v)
        integrate(function(u) distrib_pdf(d2, u, th), -pi, v)$value,
        numeric(1)))
#>                  [,1]      [,2]      [,3]      [,4]
#> series     0.04723942 0.2659828 0.5050289 0.9628184
#> quadrature 0.04723942 0.2659828 0.5050289 0.9628184

# And with the concentration parametrization at the implied concentration.
k <- numericals7::bessel_i_ratio_inverse(0.7)$kappa
all.equal(distrib_cdf(d2, y, th),
          distrib_cdf(vonmises1_distrib(), y, list(mu = 0.5, kappa = k)))
#> [1] TRUE
```
