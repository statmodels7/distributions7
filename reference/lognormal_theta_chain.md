# The Mixed Grid of the Lognormal

Computes the gaussian's components at \\t = \log y\\, carried by the
Jacobian of the transformation. It serves all three of the lognormal's
mixed methods,
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md),
with `order` selecting the response derivative and `second` the order in
\\\theta\\.

## Usage

``` r
lognormal_theta_chain(y, theta, order, second)
```

## Arguments

- y:

  A numeric vector of observations, strictly positive.

- theta:

  A named list containing `mu` and `sigma2`, read by the gaussian
  unchanged.

- order:

  `1` for the derivative of
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
  `2` for that of
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md).

- second:

  `FALSE` for one derivative in \\\theta\\, keyed by parameter, and
  `TRUE` for two, keyed by parameter pair.

## Value

A named list of numeric vectors, keyed by parameter when `second` is
`FALSE` and by parameter pair when it is `TRUE`.

## Why theta passes through untouched

With \\t = \log y\\ the log-density is the gaussian's in \\t\\ plus the
log-Jacobian, \\\ell(y) = \ell_N(t;\mu,\sigma^2) - \log y\\. The
transformation carries NO parameter, so \\t\\ does not move with
\\\theta\\ and every \\\theta\\-derivative is the gaussian's own at
\\t\\. The parent is `gaussian2`, which carries the same \\(\mu,
\sigma^2)\\ the lognormal does, so `theta` needs no map at all.

## What the Jacobian does carry

Writing \\g\\ for the gaussian's log-density in \\t\\, \$\$\ell^{(y)} =
\frac{g' - 1}{y}, \qquad \ell^{(yy)} = \frac{g'' - g' + 1}{y^2},\$\$ and
differentiating in \\\theta\\ kills the constants, leaving
\$\$\frac{\partial\ell^{(y)}}{\partial\theta} = \frac{1}{y}
\frac{\partial g'}{\partial\theta}, \qquad
\frac{\partial\ell^{(yy)}}{\partial\theta} = \frac{1}{y^2}\left(
\frac{\partial g''}{\partial\theta} - \frac{\partial
g'}{\partial\theta}\right),\$\$ with the same two lines at second order
in \\\theta\\.

## Notation

\\t = \log y\\, \\g\\ is the gaussian's log-density in \\t\\, and a
prime denotes a derivative in \\t\\.

## See also

[`distrib_cross2_y.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.Lognormal1Distrib.md),
the three methods it serves, and
[`mapped_theta2()`](https://statmodels7.github.io/distributions7/reference/mapped_theta2.md),
the route a family whose transformation DOES carry a parameter takes
instead.

## Examples

``` r
d <- lognormal1_distrib()
y <- c(0.4, 1.1, 2.3)
theta <- list(mu = 0.2, sigma2 = 0.6)

distributions7:::lognormal_theta_chain(y, theta, 2L, FALSE)
#> $mu
#> [1] -10.4166667  -1.3774105  -0.3150599
#> 
#> $sigma2
#> [1] 36.7411585  2.5360189  0.1927593
#> 

# Which is what the family's distrib_cross2_y reports.
distrib_cross2_y(d, y, theta)
#> $mu
#> [1] -10.4166667  -1.3774105  -0.3150599
#> 
#> $sigma2
#> [1] 36.7411585  2.5360189  0.1927593
#> 

# And it is the gaussian's own component at t = log y, over y^2.
g2 <- gaussian2_distrib()
t <- log(y)
first <- distrib_cross_y(g2, t, theta)
hi <- distrib_cross2_y(g2, t, theta)
(hi$mu - first$mu) / y^2
#> [1] -10.4166667  -1.3774105  -0.3150599
```
