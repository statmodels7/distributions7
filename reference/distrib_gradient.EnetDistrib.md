# Elastic-Net Score

Computes the three first derivatives of the log-density in closed form.
With \\z = y-\mu\\, \\a = \lambda\alpha\\ and \\c = \lambda(1-\alpha)\\,
\$\$\dfrac{\partial\ell}{\partial\mu} = a\\\mathrm{sgn}(z) + cz,\$\$ and
the two rate components are the data terms less the derivatives of
\\\log Z\\, which are \\G/\sqrt c\\ in \\a\\ and \\-(1+xG)/(2c)\\ in
\\c\\. The chain to \\(\lambda, \alpha)\\ is linear, the map being
bilinear.

The location component is **undefined at \\y = \mu\\**, where
\\\mathrm{sgn}\\ jumps; `sign(0)` is 0 in R, so the value returned there
is \\cz = 0\\, one point of the subdifferential. `params_smooth` records
`mu = FALSE` so that
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)'s
finite-difference guard knows.

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `lambda` and `alpha`, each a
  numeric vector of length 1 or of the length of `y`.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation to
  the link scale is applied in the generic's body, so this method always
  returns the parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `mu`, `lambda` and `alpha`, each
of the length of the recycled inputs.

## Notation

\\a = \lambda\alpha\\, \\c = \lambda(1-\alpha)\\, \\x = a/\sqrt c\\, \\G
= \mathrm{d}\log M/\mathrm{d}x\\ with \\M\\ the Mills ratio, and \\Z\\
the normalizing constant.

## See also

[`distrib_hessian.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.EnetDistrib.md)
for the second derivatives,
[`distrib_expected_hessian.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.EnetDistrib.md)
for the information, which does not agree with the observed curvature at
the kink, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- enet_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, lambda = 2, alpha = 0.5)
g <- distrib_gradient(d, y, th)

# The location component written out.
all.equal(g$mu, 1 * sign(y) + 1 * y)
#> [1] TRUE

# The two rate components against numerical differentiation.
f <- function(p) sum(distrib_pdf(d, y, list(mu = 0, lambda = p[1],
                                            alpha = p[2]), log = TRUE))
rbind(analytic = c(sum(g$lambda), sum(g$alpha)),
      numeric = numDeriv::grad(f, c(2, 0.5)))
#>               [,1]      [,2]
#> analytic -2.352365 0.6116233
#> numeric  -2.352365 0.6116233

# At the location the derivative does not exist; the value returned is
# one point of the subdifferential.
c(at_kink = distrib_gradient(d, 0, th)$mu,
  just_below = distrib_gradient(d, -1e-9, th)$mu,
  just_above = distrib_gradient(d, 1e-9, th)$mu)
#>    at_kink just_below just_above 
#>          0         -1          1 
```
