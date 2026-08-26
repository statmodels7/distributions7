# Cauchy Score

Computes the first derivatives of the Cauchy log-density with respect to
\\\mu\\ and \\\sigma\\, one value per observation, in closed form.
Writing \\r = y - \mu\\ and \\d = \sigma^2 + r^2\\, \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{2r}{d}, \qquad \dfrac{\partial
\ell}{\partial \sigma} = \dfrac{r^2 - \sigma^2}{\sigma d}.\$\$

The score in \\\mu\\ is **redescending**: it rises to \\1/\sigma\\ at
\\r = \sigma\\ and falls back towards zero as \\\|r\|\\ grows, so an
observation far from the location contributes almost nothing to the
estimating equation. That is the mechanism behind the robustness of a
Cauchy likelihood, and it is also why the log-likelihood can have
several local maxima.

## Arguments

- distrib:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of two numeric vectors, `mu` and `sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location and
\\\sigma \> 0\\ the scale. Neither is a moment: no moment of this family
exists.

## See also

[`distrib_hessian.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.CauchyDistrib.md)
for the second derivatives,
[`distrib_expected_hessian.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.CauchyDistrib.md)
for their expectation,
[`distrib_gradient.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian1Distrib.md)
for the unbounded score of a light-tailed family, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- cauchy_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out.
r <- y - 0.4; dd <- 1.5^2 + r^2
all.equal(g$mu, 2 * r / dd)
#> [1] TRUE
all.equal(g$sigma, (r^2 - 1.5^2) / (1.5 * dd))
#> [1] TRUE

# The score redescends: it peaks at r = sigma, where it is 1/sigma, and
# decays afterwards. An outlier is discounted, not chased.
r <- c(0, 0.75, 1.5, 3, 6, 12)
round(distrib_gradient(d, 0.4 + r, th)$mu, 4)
#> [1] 0.0000 0.5333 0.6667 0.5333 0.3137 0.1641
1 / 1.5
#> [1] 0.6666667

# A Gaussian score at the same residuals grows without bound instead.
round(distrib_gradient(gaussian1_distrib(), 0.4 + r, th)$mu, 4)
#> [1] 0.0000 0.3333 0.6667 1.3333 2.6667 5.3333
```
