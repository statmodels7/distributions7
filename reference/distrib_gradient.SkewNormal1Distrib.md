# Skew Normal Score

Computes the three first derivatives of the log-density in closed form.
With \\z = (y-\mu)/\sigma\\, \\t = \alpha z\\ and \\R(t) =
\phi(t)/\Phi(t)\\ the inverse Mills ratio, \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{z - \alpha R}{\sigma}, \qquad
\dfrac{\partial \ell}{\partial \sigma} = \dfrac{z^2 - 1 - \alpha z
R}{\sigma}, \qquad \dfrac{\partial \ell}{\partial \alpha} = z R.\$\$
Every derivative of \\\log\Phi(t)\\ is a polynomial in \\t\\ and \\R\\,
because \\R' = -R(t+R)\\ closes the recursion, so the whole derivative
surface of this family stays elementary.

The ratio comes from
[`numericals7::mills_ratio()`](https://statmodels7.github.io/numericals7/reference/mills_ratio.html),
which forms it on the log scale. Below about \\t = -38\\ both
\\\phi(t)\\ and \\\Phi(t)\\ underflow while their ratio is finite and
close to \\-t\\: measured at \\t = -400\\ the ratio is 400.0025 where
`dnorm(t)/pnorm(t)` is `NaN`.

At \\\alpha = 0\\ the shape score is \\z\sqrt{2/\pi}\\ and the location
score is \\z/\sigma\\, so the two are exactly proportional. The expected
information is singular at symmetry for that reason.

## Arguments

- distrib:

  A `SkewNormal1Distrib` object, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `alpha`, each a numeric
  vector of length 1 or of the length of `y`. `sigma` must be strictly
  positive.

- scale:

  Either `"parameter"`, the default, or `"link"`. On the link scale the
  derivatives are taken with respect to the unconstrained coordinates
  \\\eta\\ of the three link functions. The transformation is applied in
  the generic's body, so this method always returns the parameter scale
  and the argument is here to match the generic's signature.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `mu`, `sigma` and `alpha`, each
of the length of the recycled inputs.

## Notation

\\\ell\\ is the log-density of one observation, \\\phi\\ and \\\Phi\\
the standard Gaussian density and distribution function, and \\R(t) =
\phi(t)/\Phi(t)\\ the inverse Mills ratio.

## See also

[`distrib_hessian.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewNormal1Distrib.md)
for the second derivatives,
[`distrib_grad_y.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormal1Distrib.md)
for the derivative in the response,
[`numericals7::mills_ratio()`](https://statmodels7.github.io/numericals7/reference/mills_ratio.html)
for the ratio, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- skewnormal1_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3)
g <- distrib_gradient(d, y, th)

# The shape component written out.
all.equal(g$alpha, y * numericals7::mills_ratio(3 * y)$r)
#> [1] TRUE

# The score sums to nearly zero at the maximum likelihood estimate.
set.seed(3)
x <- distrib_rng(d, 4000, th)
fit <- fit_distrib(d, x)
vapply(distrib_gradient(d, x, as.list(coef(fit))), sum, 0) / 4000
#>            mu         sigma         alpha 
#> -5.679367e-08 -1.828155e-07  1.972956e-09 

# At symmetry the shape score is a fixed multiple of the location score,
# which is where this parametrization loses rank.
g0 <- distrib_gradient(d, y, list(mu = 0, sigma = 1.4, alpha = 0))
c(ratio = unique(round(g0$alpha / g0$mu, 12)), sigma_root_2_pi = 1.4 * sqrt(2 / pi))
#>           ratio sigma_root_2_pi 
#>        1.117038        1.117038 
```
