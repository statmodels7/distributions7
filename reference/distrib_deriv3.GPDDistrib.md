# Generalized Pareto Derivatives of Orders Three to Five

Return the third, fourth and fifth derivatives of the log-density in
\\(\sigma, \xi)\\, and for orders three and four their expectations,
each from its own compiled kernel.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` only its
  length is read.

- theta:

  A named list with components `sigma` (positive) and `xi` (any real
  value), each a numeric vector of length 1 or of the length of `y`.

- expected:

  Logical of length 1; for orders three and four, whether the
  expectation under the model is returned. Defaults to `FALSE`.

- scale:

  `"parameter"` or `"link"`; the link scale is applied by the generic
  for orders three and four and by
  [`deriv5_scale()`](https://statmodels7.github.io/distributions7/reference/deriv5_scale.md)
  for order five.

- approx, nsim:

  Accepted for the generic's signature; the expectations are closed
  forms.

- ...:

  Unused.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`. The result does not depend on the count.

## Value

A named list with one numeric vector per component, keyed as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
keys them: 4, 5 and 6 components at orders three, four and five.

## Details

With \\z = y/\sigma\\, \\t = 1 + \xi z\\ and \\u = \xi z\\, the
log-density is \\-\log\sigma - \log t - W\\ with \\W = \log(t)/\xi =
z\\\phi(u)\\ and \\\phi(u) = \log(1+u)/u\\. Every component is a closed
form derived offline and written out, \\\sigma^{-r}\\ times a form in
\\z\\, \\\xi\\ and \\1/t\\ for \\r\\ derivatives in \\\sigma\\. Since
\\\partial W/\partial z = 1/t\\, only the pure \\\xi\\ component
contains \\\phi\\, through \\\partial^j W/\partial\xi^j =
z^{j+1}\phi^{(j)}(u)\\ with \$\$\phi^{(j)}(u) = \frac{(-1)^j
j!}{u^{j+1}} \sum\_{i \> j} \frac{v^i}{i}, \qquad v = \frac{u}{1+u},\$\$
summed as a series in \\v\\ for \\\|v\| \le 3/4\\ and as \\\log t\\
minus the first \\j\\ terms elsewhere, so that no form divides by
\\\xi\\ and \\\xi = 0\\ is an ordinary point. \\t\\ is formed from the
exact product \\\xi y\\, which keeps its relative accuracy near the
upper end of the support when \\\xi \< 0\\.

The expected derivatives follow from \\V = (1 + \xi
Y/\sigma)^{-1/\xi}\\, uniform on \\(0, 1)\\: each is a rational function
of \\\xi\\ over \\\sigma^r\\. They exist for \\\xi \> -1/3\\ at order
three and \\\xi \> -1/4\\ at order four, and are `NA` below; `approx`
and `nsim` are not read.

Outside the support, \\y \< 0\\ or \\1 + \xi y/\sigma \le 0\\, every
component is `NaN`.

## See also

[`distrib_hessian.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GPDDistrib.md)
for the order below,
[`gpd_components()`](https://statmodels7.github.io/distributions7/reference/gpd_components.md)
for an independent assembly of orders one to four, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- gpd_distrib()
y <- c(0.3, 1.5, 6)
th <- list(sigma = 1.5, xi = 0.3)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "sigma_sigma_sigma" "sigma_sigma_xi"    "sigma_xi_xi"      
#> [4] "xi_xi_xi"         

# A central difference of the Hessian reproduces the pure shape component.
eps <- 1e-5
up <- distrib_hessian(d, y, list(sigma = 1.5, xi = 0.3 + eps))$xi_xi
dn <- distrib_hessian(d, y, list(sigma = 1.5, xi = 0.3 - eps))$xi_xi
all.equal((up - dn) / (2 * eps), d3$xi_xi_xi, tolerance = 1e-6)
#> [1] TRUE

# The shape zero is an ordinary point.
distrib_deriv4(d, y, list(sigma = 1.5, xi = 0))$xi_xi_xi_xi
#> [1]  8.0640e-03  1.2000e+00 -3.3792e+03

# The expected third derivatives, NA where they do not exist.
distrib_deriv3(d, 1, list(sigma = 1.5, xi = c(0.2, -0.4)),
               expected = TRUE)$xi_xi_xi
#> [1] 8.928571       NA
```
