# von Mises Score in the Resultant Length

Computes the first derivatives of the log-density with respect to the
mean direction \\\mu\\ and the mean resultant length \\\rho\\, one value
per observation, in closed form: \$\$\dfrac{\partial\ell}{\partial\mu} =
\kappa\sin(y-\mu), \qquad \dfrac{\partial\ell}{\partial\rho} =
\left\\\cos(y-\mu) - A(\kappa)\right\\\kappa'(\rho),\$\$ with \\\kappa =
A^{-1}(\rho)\\ and \\\kappa'(\rho) = 1/A'(\kappa)\\ from the inverse
function rule.

The map touches the **second parameter only**, so the chain rule is the
one-variable one: the direction's component is unchanged from the
concentration parametrization, and the second is that family's
multiplied by a single factor. No multivariate expansion and no
cancellation are involved at any order.

With `scale = "link"` the generic applies the chain rule for the links
the family carries. This method always returns the parameter scale.

## Arguments

- distrib:

  A `VonMises2Distrib` object, from
  [`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md).

- y:

  A numeric vector of angles in \\\[-\pi, \pi)\\.

- theta:

  A named list with components `mu` and `rho`, each a numeric vector of
  length 1 or of the length of `y`. `mu` must lie in \\(-\pi, \pi)\\ and
  `rho` in \\(0, 1)\\.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, `mu` and `rho`, each of length
`max(length(y), length(mu), length(rho))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean
direction, \\\rho \in (0,1)\\ the mean resultant length, \\\kappa\\ the
concentration and \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\, so that
\\\rho = A(\kappa)\\.

## See also

[`distrib_hessian.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.VonMises2Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.VonMises2Distrib.md)
for their expectation,
[`distrib_gradient.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.VonMises1Distrib.md)
for the same quantity in the concentration, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d2 <- vonmises2_distrib()
y <- c(-1, 0, 0.5, 2)
th <- list(mu = 0.5, rho = 0.7)
g <- distrib_gradient(d2, y, th)

# numDeriv on the summed log-density reproduces the summed score.
fn <- function(v)
  sum(distrib_pdf(d2, y, list(mu = v[1], rho = v[2]), log = TRUE))
rbind(numeric = numDeriv::grad(fn, c(0.5, 0.7)),
      closed = vapply(g, sum, numeric(1)))
#>                 mu       rho
#> numeric -0.9653846 -4.809689
#> closed  -0.9653846 -4.809689

# The direction component is unchanged from the concentration
# parametrization, the map touching the second parameter only.
k <- numericals7::bessel_i_ratio_inverse(0.7)$kappa
all.equal(g$mu,
          distrib_gradient(vonmises1_distrib(), y,
                           list(mu = 0.5, kappa = k))$mu)
#> [1] TRUE
```
