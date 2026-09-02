# Skew Normal Observed Hessian

Computes the six second derivatives of the log-density in closed form.
In the notation of
[`distrib_gradient.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal1Distrib.md),
with \\R' = -R(t+R)\\, \$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
\dfrac{\alpha^2 R' - 1}{\sigma^2}, \qquad \dfrac{\partial^2
\ell}{\partial \mu \\ \partial \sigma} = \dfrac{\alpha^2 z R' - 2z +
\alpha R}{\sigma^2}, \qquad \dfrac{\partial^2 \ell}{\partial \mu \\
\partial \alpha} = -\dfrac{R + \alpha z R'}{\sigma},\$\$
\$\$\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{1 - 3z^2 +
2\alpha z R + \alpha^2 z^2 R'}{\sigma^2}, \qquad \dfrac{\partial^2
\ell}{\partial \sigma \\ \partial \alpha} = -\dfrac{z R + \alpha z^2
R'}{\sigma}, \qquad \dfrac{\partial^2 \ell}{\partial \alpha^2} = z^2
R'.\$\$

This is the **observed** curvature at the data. The family has no
elementary expected information, so
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
falls to the base class and approximates it; see there for the
strategies and their cost.

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

  Either `"parameter"`, the default, or `"link"`. The transformation to
  the link scale is applied in the generic's body, so this method always
  returns the parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of six numeric vectors, in the order `mu_mu`,
`sigma_sigma`, `alpha_alpha`, `mu_sigma`, `mu_alpha`, `sigma_alpha`,
each of the length of the recycled inputs. The names are the diagonal
first, then the off-diagonal pairs, which is
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
ordering.

## Singularity at symmetry

At \\\alpha = 0\\ the expected information of this parametrization has
rank 2, not 3. The reason is in the score: the shape and location
components are exactly proportional there, so no data can separate them.
Measured on the approximated information at \\\mu = 0\\, \\\sigma = 1\\,
its eigenvalues are 2, 1.637 and \\-5.6\times10^{-17}\\, and the
smallest one grows like \\\alpha^4\\ as the shape moves off zero:
\\4.4\times10^{-10}\\ at \\\alpha = 0.01\\ and \\1.9\times10^{-3}\\ at
\\\alpha = 0.5\\.

The singularity belongs to the family as parametrized, not to the
implementation. Azzalini's centered parametrization removes it, and
lives in the separate object
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

## Notation

\\z = (y-\mu)/\sigma\\, \\t = \alpha z\\, \\R(t) = \phi(t)/\Phi(t)\\ and
\\R' = -R(t+R)\\.

## See also

[`distrib_gradient.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal1Distrib.md)
for the score,
[`distrib_deriv3.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewNormal1Distrib.md)
for the next order,
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the approximated expectation, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- skewnormal1_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3)
h <- distrib_hessian(d, y, th)

# The shape component written out.
all.equal(h$alpha_alpha, y^2 * numericals7::mills_ratio(3 * y)$dr)
#> [1] TRUE

# Against a central difference of the score.
eps <- 1e-5
rbind(analytic = h$mu_alpha,
      numeric = (distrib_gradient(d, y, list(mu = 0, sigma = 1, alpha = 3 + eps))$mu -
                 distrib_gradient(d, y, list(mu = 0, sigma = 1, alpha = 3 - eps))$mu) /
                (2 * eps))
#>               [,1]      [,2]     [,3]         [,4]
#> analytic -9.029656 -2.155568 0.154335 3.714795e-08
#> numeric  -9.029656 -2.155568 0.154335 3.714806e-08

# The information loses rank at symmetry, and recovers it as alpha^4.
rank_gap <- function(a) {
  e <- distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1, alpha = a))
  M <- matrix(c(e$mu_mu, e$mu_sigma, e$mu_alpha,
                e$mu_sigma, e$sigma_sigma, e$sigma_alpha,
                e$mu_alpha, e$sigma_alpha, e$alpha_alpha), 3, 3)
  min(abs(eigen(-M, only.values = TRUE)$values))
}
vapply(c(0, 0.01, 0.5), rank_gap, 0)
#> [1] 5.738220e-27 4.355679e-10 1.947533e-03
```
