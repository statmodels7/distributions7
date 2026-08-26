# Pseudo-Huber Cumulative Distribution Function

Computes \\F(q) = P(Y \le q)\\ by numerical integration of the density.
The family has no elementary distribution function, so there is nothing
closed form to call.

Two devices keep the quadrature honest. The law is symmetric about
\\\mu\\, so a quantile above the location is **reflected**, \\F(q) = 1 -
F(2\mu - q)\\, and only the lower tail is ever integrated, where the
integrand decays away from a finite endpoint. And every quantile is one
**row** of a single batched quadrature through
[`quad_rows()`](https://statmodels7.github.io/distributions7/reference/quad_rows.md),
so a vector of `q` is integrated in a single call.

A row that fails to reach the requested accuracy signals an error naming
the positions, instead of returning a plausible number.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `q`. A component of length 1 is
  recycled; a vector gives one integration per parameter setting.
  `sigma` and `nu` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are \\P(Y
  \le q)\\; when `FALSE` they are \\P(Y \> q)\\, formed as \\1 - F\\.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned, taken after the quadrature, so it carries the quadrature's
  own accuracy rather than improving on it in the far tail. Defaults to
  `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(q), length(mu), length(sigma), length(nu))`, clamped to that
range.

## See also

[`distrib_quantile.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PseudoHuberDistrib.md)
for the inverse,
[`distrib_pdf.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.PseudoHuberDistrib.md)
for the integrand,
[`distrib_grad_cdf.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PseudoHuberDistrib.md)
for the derivatives of this function, which are closed form in the
location and the scale, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
th <- list(mu = 0.4, sigma = 1.2, nu = 2)

# The quadrature, and the symmetry the method exploits.
distrib_cdf(d, c(-1, 0.4, 2), th)
#> [1] 0.2223070 0.5000000 0.8061021
distrib_cdf(d, 0.4 - 1.5, th) + distrib_cdf(d, 0.4 + 1.5, th)
#> [1] 1

# At the location it is one half exactly, the law being symmetric.
distrib_cdf(d, 0.4, th)
#> [1] 0.5

# It agrees with a direct integration of the density.
c(method = distrib_cdf(d, 2, th),
  integral = integrate(function(v) distrib_pdf(d, v, th), -Inf, 2)$value)
#>    method  integral 
#> 0.8061021 0.8061021 
```
