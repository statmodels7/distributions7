# Student t Score

Computes the first derivatives of the location-scale Student t
log-density with respect to \\\mu\\, \\\sigma\\ and \\\nu\\, one value
per observation, in closed form. Writing \\r = y - \mu\\ and \\D =
\nu\sigma^2 + r^2\\, \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{(\nu+1)r}{D}, \qquad \dfrac{\partial \ell}{\partial \sigma} =
\dfrac{\nu\left(r^2 - \sigma^2\right)}{\sigma D},\$\$
\$\$\dfrac{\partial \ell}{\partial \nu} = \dfrac{1}{2}\left\[
-\dfrac{1}{\nu} - \psi\\\left(\dfrac{\nu}{2}\right) +
\psi\\\left(\dfrac{\nu+1}{2}\right) + \dfrac{(\nu+1)r^2}{\nu D} -
\log\\\left(1 + \dfrac{r^2}{\nu\sigma^2}\right)\right\].\$\$

The location component **redescends**: it rises to
\\(\nu+1)/(2\sigma\sqrt{\nu})\\ at \\\|r\| = \sigma\sqrt{\nu}\\ and
falls back towards zero beyond it, so a gross outlier contributes almost
nothing to the estimating equation. A Gaussian's location score grows
without bound instead, which is the whole reason this family is used for
robust location estimation.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale; the transformation happens in the generic. The
arithmetic runs in a compiled kernel decomposed over the elements of the
output, so the result does not depend on the thread count.

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive.

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

A named list of three numeric vectors, `mu`, `sigma` and `nu`, each of
length `max(length(y), length(mu), length(sigma), length(nu))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location,
\\\sigma \> 0\\ the scale and \\\nu \> 0\\ the degrees of freedom.
\\\psi\\ is the digamma function,
[`digamma()`](https://rdrr.io/r/base/Special.html) in R.

## Large degrees of freedom

Every component is written as a ratio in \\z = r/\sigma\\ and \\u =
z^2/\nu\\ instead of as a quotient by powers of \\D\\, because
\\\nu\sigma^2\\ overflows well before the log link's own clamp is
reached. All three components stay finite at every \\\nu\\ the chart can
produce, up to `.Machine$double.xmax`.

## See also

[`distrib_hessian.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentT1Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.StudentT1Distrib.md)
for their expectation,
[`distrib_grad_y.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.StudentT1Distrib.md)
for the derivative in the response, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- student_t1_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 5)
g <- distrib_gradient(d, y, th)

# The location and scale components, written out.
r <- y - 0.4; D <- 5 * 1.2^2 + r^2
all.equal(g$mu, 6 * r / D)
#> [1] TRUE
all.equal(g$sigma, 5 * (r^2 - 1.2^2) / (1.2 * D))
#> [1] TRUE

# The location score redescends: it peaks at |r| = sigma * sqrt(nu) = 2.68
# and falls back, where a Gaussian's grows without bound.
rr <- c(0.5, 1, 2, 4, 8, 16)
rbind(residual = rr,
      t = 6 * rr / (5 * 1.2^2 + rr^2),
      gaussian = rr / 1.2^2)
#>               [,1]      [,2]     [,3]     [,4]      [,5]       [,6]
#> residual 0.5000000 1.0000000 2.000000 4.000000 8.0000000 16.0000000
#> t        0.4026846 0.7317073 1.071429 1.034483 0.6741573  0.3647416
#> gaussian 0.3472222 0.6944444 1.388889 2.777778 5.5555556 11.1111111

# The summed score vanishes at the maximum likelihood estimate.
set.seed(4)
z <- distrib_rng(d, 3000, list(mu = 1, sigma = 2, nu = 4))
mle <- as.list(coef(fit_distrib(d, z)))
vapply(distrib_gradient(d, z, mle), sum, numeric(1))
#>           mu        sigma           nu 
#> 1.787019e-04 3.952579e-04 9.730911e-05 
```
