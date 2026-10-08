# Pseudo-Huber Derivatives in the Parameters, Standard-Deviation Parametrization

Return the derivatives of the log-density in \\(\mu, \sigma, \nu)\\ of
orders one to five, in closed form, each order by its own compiled
kernel.

## Arguments

- distrib:

  A `PseudoHuber2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- scale:

  `"parameter"` or `"link"`; the link scale is applied by the generic.

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
\sigma^{-(j+2)}\\. A component of \\D\\ is Faa di Bruno's sum over the
set partitions of its indices, with \\(d/dQ)^m \sqrt{Q} = c_m Q^{1/2 -
m}\\. The derivatives of \\R\\ and \\h\\ are read off \\\log
K_n(\sqrt\nu)\\, whose derivatives in \\t = \sqrt\nu\\ follow from the
ratios \\K_n^{(j)}(t)/K_n(t) = (-1/2)^j \sum_i \binom{j}{i}
K\_{\|n-j+2i\|}(t)/ K_n(t)\\, formed with exponentially scaled Bessel
functions. A kernel of order \\m\\ computes \\R\\ and \\h\\ to order
\\m\\ and the partitions of order \\m\\ only.

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
# against a central difference in sigma
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
