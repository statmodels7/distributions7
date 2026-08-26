# von Mises Score

Computes the first derivatives of the von Mises log-density with respect
to the mean direction \\\mu\\ and the concentration \\\kappa\\, one
value per observation, in closed form: \$\$\dfrac{\partial
\ell}{\partial \mu} = \kappa \sin(y - \mu), \qquad \dfrac{\partial
\ell}{\partial \kappa} = \cos(y - \mu) - A(\kappa),\$\$ with \\A(\kappa)
= I_1(\kappa)/I_0(\kappa)\\ the derivative of \\\log I_0\\ and also the
**mean resultant length** of the family.

The ratio comes from
[`numericals7::bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.html),
which switches to an asymptotic expansion past \\\kappa = 10^4\\. R's
own scaled `besselI` underflows to an exact zero between \\10^5\\ and
\\10^6\\, so forming the ratio from two calls gives `NaN` over part of
that band.

With `scale = "link"` the generic applies the chain rule for the links
the family carries. This method always returns the parameter scale.

## Arguments

- distrib:

  A `VonMises1Distrib` object, from
  [`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md).

- y:

  A numeric vector of angles in \\\[-\pi, \pi)\\.

- theta:

  A named list with components `mu` and `kappa`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `mu` must lie in \\(-\pi, \pi)\\ and `kappa` be strictly
  positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, `mu` and `kappa`, each of length
`max(length(y), length(mu), length(kappa))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean
direction, \\\kappa \> 0\\ the concentration, \\I_m\\ the modified
Bessel function of the first kind of order \\m\\, and \\A(\kappa) =
I_1(\kappa)/I_0(\kappa)\\.

## See also

[`distrib_hessian.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.VonMises1Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.VonMises1Distrib.md)
for their expectation,
[`numericals7::bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.html)
for \\A\\, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- vonmises1_distrib()
y <- c(-1, 0, 0.5, 2)
th <- list(mu = 0.5, kappa = 2)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out.
all.equal(g$mu, 2 * sin(y - 0.5))
#> [1] TRUE
all.equal(g$kappa, cos(y - 0.5) - besselI(2, 1) / besselI(2, 0))
#> [1] TRUE

# numDeriv on the summed log-density reproduces the summed score.
fn <- function(p)
  sum(distrib_pdf(d, y, list(mu = p[1], kappa = p[2]), log = TRUE))
rbind(numeric = numDeriv::grad(fn, c(0.5, 2)),
      closed = vapply(g, sum, numeric(1)))
#>                 mu      kappa
#> numeric -0.9588511 -0.7720417
#> closed  -0.9588511 -0.7720417

# The summed score vanishes at the maximum likelihood estimate.
set.seed(7)
z <- distrib_rng(d, 2000, list(mu = 0.8, kappa = 3))
mle <- as.list(coef(fit_distrib(d, z)))
vapply(distrib_gradient(d, z, mle), sum, numeric(1))
#>            mu         kappa 
#> -5.080589e-13  2.680249e-08 
```
