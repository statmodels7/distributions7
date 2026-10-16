# Pseudo-Huber Derivatives in the Parameters, Standard-Deviation Parametrization

Return the derivatives of the log-density in \\(\mu, \sigma, \nu)\\ of
orders one to five, in closed form; a compiled kernel of its own
computes each order.

## Arguments

- distrib:

  A `PseudoHuber2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- scale:

  `"parameter"` or `"link"`. The generic applies the link scale at
  orders one to four;
  [`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
  applies it in the method, through
  [`deriv5_scale()`](https://statmodels7.github.io/distributions7/reference/deriv5_scale.md).

- expected:

  Logical; for orders three and four, whether the expected derivative is
  returned.

- approx, nsim:

  The approximation and Monte Carlo size of
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md).

- ...:

  Unused.

## Value

A named list with one numeric vector per component, named by the
parameters joined with `_`: 3 components at order one, 6 at order two
(`mu_mu`, `sigma_sigma`, `nu_nu`, `mu_sigma`, `mu_nu`, `sigma_nu`), and
10, 15 and 21 at orders three to five in lexicographic order.

## Details

With \\r = y - \mu\\ the log-density is \\\ell = -D - \log 2 -
\log\sigma + h(\nu)\\, \\D = \sqrt{Q}\\, \\Q = \nu + R(\nu)
r^2/\sigma^2\\. \\Q\\ is a product of one-variable functions, so each of
its partial derivatives is closed: \$\$\partial\_\mu^i
\partial\_\sigma^j \partial\_\nu^k Q = \[i = j = 0,\\ k = 1\] +
R^{(k)}(\nu)\\ B_i\\ C_j,\$\$ with \\B_0 = r^2\\, \\B_1 = -2r\\, \\B_2 =
2\\, \\B_i = 0\\ for \\i \ge 3\\, and \\C_j = (-1)^j (j+1)!\\
\sigma^{-(j+2)}\\. A component of \\D\\ of order one and two is written
out per component. For orders three to five the data part \\-A/\sigma\\,
with \\A = \sqrt{\nu\sigma^2 + R(\nu) r^2}\\, is differentiated in
closed form, in expressions written in \\u = r/A\\ and powers of \\1/A\\
that are tabulated in the C++ sources. The derivatives of \\R\\ and \\h
= \log(R)/2 - B\\ come from those of \\B(\nu) = \log\\\sqrt\nu
K_1(\sqrt\nu)\\\\ and \\\log R\\, computed by three routes: the series
of \\K_1\\ and \\K_2\\ at the origin below \\\nu = 1\\, the cumulants of
the generalized inverse Gaussian law of their integral representation
for \\1 \le \nu \< 400\\, and the asymptotic series from \\\nu = 400\\.
A kernel of order \\m\\ computes \\R\\ and \\h\\ to order \\m\\ and
nothing above it.

The expected derivatives of orders three and four are those of
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md);
the expected information is registered through
[`register_loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/register_loc_scale_expected.md).

## See also

[`distrib_pdf.PseudoHuber2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.PseudoHuber2Distrib.md).

## Examples

``` r
d <- pseudohuber2_distrib()
th <- list(mu = 0.4, sigma = 1.5, nu = 2)
y <- c(-1, 0.5, 3)
g <- distrib_gradient(d, y, th)
# The score in sigma is compared with a central difference.
h <- 1e-6
up <- distrib_pdf(d, y, list(mu = 0.4, sigma = 1.5 + h, nu = 2), log = TRUE)
dn <- distrib_pdf(d, y, list(mu = 0.4, sigma = 1.5 - h, nu = 2), log = TRUE)
all.equal(g$sigma, (up - dn) / (2 * h), tolerance = 1e-6)
#> [1] TRUE
names(distrib_deriv5(d, y, th))
#>  [1] "mu_mu_mu_mu_mu"                "mu_mu_mu_mu_sigma"            
#>  [3] "mu_mu_mu_mu_nu"                "mu_mu_mu_sigma_sigma"         
#>  [5] "mu_mu_mu_sigma_nu"             "mu_mu_mu_nu_nu"               
#>  [7] "mu_mu_sigma_sigma_sigma"       "mu_mu_sigma_sigma_nu"         
#>  [9] "mu_mu_sigma_nu_nu"             "mu_mu_nu_nu_nu"               
#> [11] "mu_sigma_sigma_sigma_sigma"    "mu_sigma_sigma_sigma_nu"      
#> [13] "mu_sigma_sigma_nu_nu"          "mu_sigma_nu_nu_nu"            
#> [15] "mu_nu_nu_nu_nu"                "sigma_sigma_sigma_sigma_sigma"
#> [17] "sigma_sigma_sigma_sigma_nu"    "sigma_sigma_sigma_nu_nu"      
#> [19] "sigma_sigma_nu_nu_nu"          "sigma_nu_nu_nu_nu"            
#> [21] "nu_nu_nu_nu_nu"               
```
