# Weibull Derivatives in the Mean and Shape

Return the derivatives of the log-density in \\(m, \sigma)\\ of orders
one to five, the expected information and its expected third and fourth
derivatives, and the first two derivatives of the expected information,
each from its own compiled kernel.

## Arguments

- distrib:

  A `Weibull3Distrib` object.

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mean` and `sigma`.

- scale:

  `"parameter"` or `"link"`; the link scale is applied by the generic.

- expected:

  Logical; for orders three and four, whether the expected derivative is
  returned.

- approx, nsim:

  Accepted for the generic's signature; the expectations are closed.

- ..., threads:

  Unused.

## Value

A named list with one numeric vector per component: 2, 3, 4, 5 and 6
components at orders one to five, the second order keyed `mean_mean`,
`sigma_sigma`, `mean_sigma`.

## Details

With \\w = \log(y/m) + \log\Gamma(1 + 1/\sigma) = \log(y/b)\\ the
log-density is \$\$\ell = \log\sigma - \log y + \sigma w - e^{\sigma
w}.\$\$ Every component is a closed form in \\\log(y/m)\\,
\\\log\Gamma(1 + 1/\sigma)\\ and the polygamma functions at \\1 +
1/\sigma\\, derived offline and written out per component, one kernel
per order. \\T = e^{\sigma w}\\ is a standard exponential, so an
expected derivative replaces each \\T^j (\log T)^k\\ by
\\\Gamma^{(k)}(j + 1)\\.

## See also

[`distrib_pdf.Weibull3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Weibull3Distrib.md).

## Examples

``` r
d <- weibull3_distrib()
th <- list(mean = 4, sigma = 1.7)
y <- c(0.7, 2.5, 6)
g <- distrib_gradient(d, y, th)
h <- 1e-6
up <- distrib_pdf(d, y, list(mean = 4, sigma = 1.7 + h), log = TRUE)
dn <- distrib_pdf(d, y, list(mean = 4, sigma = 1.7 - h), log = TRUE)
all.equal(g$sigma, (up - dn) / (2 * h), tolerance = 1e-6)
#> [1] TRUE
```
