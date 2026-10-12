# Lognormal Derivatives in the Mean and Variance of Y

Return the derivatives of the log-density in \\(m, v)\\ of orders one to
five, the expected information and its expected third and fourth
derivatives, and the first two derivatives of the expected information,
each from its own compiled kernel.

## Arguments

- distrib:

  A `Lognormal2Distrib` object.

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mean` and `var`.

- scale:

  `"parameter"` or `"link"`. For orders one to four the generic applies
  the link scale after dispatch; the method for
  [`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
  applies it itself.

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
`var_var`, `mean_var`.

## Details

With \\L = \log y\\, \\S = \log(1 + v/m^2)\\ and \\\mu_l = \log m -
S/2\\ the log-density is \$\$\ell = -L - \tfrac12 \log 2\pi - \tfrac12
\log S - \frac{(L - \mu_l)^2}{2S}.\$\$ The observed derivatives of
orders one and two, and every expected derivative, are closed forms in
\\(L, S, m, v)\\, derived offline and written out per component, one
kernel per order. \\\ell\\ is a polynomial of degree two in \\L\\, so
each expected derivative replaces \\L\\ and \\L^2\\ by \\E\[L\] =
\mu_l\\ and \\E\[L^2\] = \mu_l^2 + S\\.

The observed derivatives of orders three to five are computed from
\\\ell = h_1(S) + h_2(S)\\C^2 - C/2 - L - \tfrac12\log 2\pi\\, with \\C
= \log(y/m)\\, \\h_1(S) = -\tfrac12\log S - S/8\\ and \\h_2(S) =
-1/(2S)\\: the derivatives of \\S = \log(1 + v/m^2)\\ in \\(m, v)\\ and
of \\h_1\\ and \\h_2\\ of \\S\\ by the formula of Faa di Bruno over the
set partitions of the differentiation indices, and those of \\C^2\\ by
the Leibniz rule. Each factor keeps the size of the quantity that it
represents when \\v/m^2\\ is small, where the expanded closed forms
cancel.

## See also

[`distrib_pdf.Lognormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Lognormal2Distrib.md).

## Examples

``` r
d <- lognormal2_distrib()
th <- list(mean = 3, var = 2)
y <- c(0.7, 2.5, 6)
g <- distrib_gradient(d, y, th)
h <- 1e-6
up <- distrib_pdf(d, y, list(mean = 3 + h, var = 2), log = TRUE)
dn <- distrib_pdf(d, y, list(mean = 3 - h, var = 2), log = TRUE)
all.equal(g$mean, (up - dn) / (2 * h), tolerance = 1e-6)
#> [1] TRUE
```
