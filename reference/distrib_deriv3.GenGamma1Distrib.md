# Generalized Gamma Derivatives of Orders Three to Five

Return the third, fourth and fifth derivatives of the log-density in
\\(a, d, p)\\, and for orders three and four their expectations, each
from its own compiled kernel.

## Arguments

- distrib:

  A `GenGamma1Distrib` object, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- y:

  A numeric vector of strictly positive observations. With
  `expected = TRUE` only its length is read.

- theta:

  A named list with components `a`, `d` and `p`, each a numeric vector
  of length 1 or of the length of `y`, all strictly positive.

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
keys them: 10, 15 and 21 components at orders three, four and five.

## Details

With \\k = d/p\\ and \\U = p\log(y/a)\\, the log-density is \$\$\ell =
\log p - \log y + kU - e^{U} - \log\Gamma(k),\$\$ and \\e^{U} =
(y/a)^p\\ is a gamma variable with shape \\k\\ and unit rate. Every
component is a closed form derived offline and written out: a polynomial
in \\U\\, \\e^{U}\\ and \\a/y\\ whose coefficients depend on \\(d, p)\\
alone, multiplied by \\a^{-r}\\ for \\r\\ derivatives in \\a\\. The
coefficients contain the polygamma functions at \\k + 1\\, through
\\\psi^{(n)}(k) = \psi^{(n)}(k+1) + (-1)^{n+1} n!/k^{n+1}\\, so that the
poles in \\1/k\\ which the derivatives combine cancel in the closed form
and not in floating point. The expected derivatives are polynomials in
the central moments of \\U\\, whose cumulants are the polygamma
functions at \\k\\; `approx` and `nsim` are not read.

## See also

[`distrib_hessian.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GenGamma1Distrib.md)
for the order below,
[`gengamma_components()`](https://statmodels7.github.io/distributions7/reference/gengamma_components.md)
for an independent assembly of orders one to four, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- gengamma1_distrib()
y <- c(0.6, 1.4, 3.1)
th <- list(a = 2, d = 1.5, p = 1.3)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#>  [1] "a_a_a" "a_a_d" "a_a_p" "a_d_d" "a_d_p" "a_p_p" "d_d_d" "d_d_p" "d_p_p"
#> [10] "p_p_p"

# A central difference of the Hessian reproduces a mixed component.
eps <- 1e-5
up <- distrib_hessian(d, y, list(a = 2, d = 1.5, p = 1.3 + eps))$d_p
dn <- distrib_hessian(d, y, list(a = 2, d = 1.5, p = 1.3 - eps))$d_p
all.equal((up - dn) / (2 * eps), d3$d_p_p, tolerance = 1e-6)
#> [1] TRUE

# The expected fourth derivatives, closed form.
distrib_deriv4(d, 1, th, expected = TRUE)$p_p_p_p
#> [1] -13.64482

# Order five on the link scale.
distrib_deriv5(d, y, th, scale = "link")$p_p_p_p_p
#> [1]   6.920295   4.371609 -14.937992
```
