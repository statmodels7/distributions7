# Generalized Gamma Derivatives in the Mean

Return the derivatives of the log-density in \\(m, d, p)\\ of orders one
to five, the expected information and its expected third and fourth
derivatives, and the first two derivatives of the expected information,
each from its own compiled kernel.

## Arguments

- distrib:

  A `GenGamma2Distrib` object.

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mean`, `d` and `p`.

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

A named list with one numeric vector per component: 3, 6, 10, 15 and 21
components at orders one to five, the second order keyed `mean_mean`,
`d_d`, `p_p`, `mean_d`, `mean_p`, `d_p`.

## Details

With \\k = d/p\\, \\k_1 = (d+1)/p\\ and \\w = \log(y/m) -
\log\Gamma(k) + \log\Gamma(k_1) = \log(y/a)\\, the log-density is
\$\$\ell = \log p - \log y + d\\w - \log\Gamma(k) - e^{p w}.\$\$ Every
component is a closed form derived offline and written out per
component, one kernel per order. With \\U = p\\w\\, the data enter
through \\X = U - \psi(k)\\ and \\Q = e^U/k - 1\\, and each component is
a polynomial in them whose coefficients depend on \\(d, p)\\ alone,
through \\\psi^{(n)}(k)\\ and the differences \\\psi^{(n)}(k_1) -
\psi^{(n)}(k)\\. As \\k\\ and \\d\\ grow the terms of these coefficients
agree to several orders, so above \\k = 10\\ or \\d = 10\\ the
coefficients are formed in double-double arithmetic (about 32 digits)
and rounded once; below, in double. Where \\\|X + \psi(k) - \log k\| \<
1\\ the polynomial is written in that variable, so that the cancellation
between \\X\\ and \\Q\\ is also carried out on the coefficients. \\e^U\\
is a gamma variable with shape \\k\\ and unit rate, so the expected
derivatives are polynomials in the central moments of its logarithm,
whose cumulants are the polygamma functions at \\k\\.

## See also

[`distrib_pdf.GenGamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GenGamma2Distrib.md).

## Examples

``` r
d <- gengamma2_distrib()
th <- list(mean = 5, d = 3, p = 1.5)
y <- c(1.2, 4, 9)
g <- distrib_gradient(d, y, th)
h <- 1e-6
up <- distrib_pdf(d, y, list(mean = 5, d = 3, p = 1.5 + h), log = TRUE)
dn <- distrib_pdf(d, y, list(mean = 5, d = 3, p = 1.5 - h), log = TRUE)
all.equal(g$p, (up - dn) / (2 * h), tolerance = 1e-6)
#> [1] TRUE
```
