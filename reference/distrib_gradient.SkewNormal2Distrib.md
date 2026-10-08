# Skew Normal Derivatives in the Centered Parametrization

Return the derivatives of the log-density in \\(\mu, \sigma, \gamma_1)\\
of orders one to five, each from its own compiled kernel.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `gamma1`.

- scale:

  `"parameter"` or `"link"`; the link scale is applied by the generic
  (by
  [`deriv5_scale()`](https://statmodels7.github.io/distributions7/reference/deriv5_scale.md)
  at the fifth order).

- expected:

  Logical; for orders three and four, whether the expected derivative is
  returned, by
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md).

- approx, nsim:

  Passed to
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
  when `expected` is `TRUE`.

- ...:

  Unused.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list with one numeric vector per component: 3, 6, 10, 15 and 21
components at orders one to five.

## Details

With \\w = (y - \mu)/\sigma\\, \\c = (4 - \pi)/2\\ and \\r =
\mathrm{sign}(\gamma_1)(\|\gamma_1\|/c)^{1/3}\\, the log-density is
\$\$\ell = -\log\sigma - \tfrac12\log(1 + r^2) - \tfrac12 z^2 +
\log\Phi(x), \quad z = \frac{w + r}{\sqrt{1 + r^2}}, \quad x =
\frac{r\\z}{\sqrt{b^2 - (1 - b^2) r^2}},\$\$ with \\b = \sqrt{2/\pi}\\,
and a derivative in \\\gamma_1\\ is \\(3 c r^2)^{-1}\partial_r\\. Every
component is a combination of the derivatives of \\F = \ell +
\log\sigma\\ in \\w\\ and \\\gamma_1\\, derived offline. Their closed
forms in \\x\\, \\w\\, \\r\\ and \\\phi(x)/\Phi(x)\\ cancel terms of
order \\r^{-2k}\\ as \\\gamma_1 \to 0\\, so where \\\|x\| \< 0.4\\ and
\\\|r\| \< 0.4\\ the kernels take the derivatives in \\\gamma_1\\ from
the series \\F = \sum\_{n \ge 3} F_n(w) r^n\\, whose coefficients are
polynomials in \\w\\ computed offline at 60 digits.

The series has no \\r\\ and no \\r^2\\ term, so the score is finite at
\\\gamma_1 = 0\\; its \\r^4\\ and \\r^5\\ terms do not vanish, so the
derivatives of order two or more in \\\gamma_1\\ diverge there and the
methods of order two and more signal an error at zero skewness.

## See also

[`distrib_expected_hessian.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal2Distrib.md).

## Examples

``` r
d <- skewnormal2_distrib()
th <- list(mu = 0.2, sigma = 1.3, gamma1 = 0.4)
y <- c(-1.7, 0.3, 2.4)
distrib_gradient(d, y, th)
#> $mu
#> [1] -1.425500  0.243771  1.049694
#> 
#> $sigma
#> [1]  1.3141921 -0.7504792  1.0071750
#> 
#> $gamma1
#> [1]  0.04328132 -0.10028722 -0.02941041
#> 
# the score is finite at zero skewness
distrib_gradient(d, y, list(mu = 0.2, sigma = 1.3, gamma1 = 0))$gamma1
#> [1]  0.21043848 -0.03838568 -0.03838568
```
