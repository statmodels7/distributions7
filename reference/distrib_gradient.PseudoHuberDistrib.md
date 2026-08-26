# Pseudo-Huber Score

Computes the first derivatives of the pseudo-Huber log-density with
respect to \\\mu\\, \\\sigma\\ and \\\nu\\, one value per observation,
in closed form. With \\r = y - \mu\\ and \\D = \sqrt{\nu +
(r/\sigma)^2}\\, \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{r}{\sigma^2 D}, \qquad \dfrac{\partial \ell}{\partial \sigma} =
\dfrac{1}{\sigma}\left(\dfrac{r^2}{\sigma^2 D} - 1\right),\$\$
\$\$\dfrac{\partial \ell}{\partial \nu} = -\dfrac{1}{2}\left\[
\dfrac{1}{\nu} + \dfrac{1}{D} + \dfrac{K_1'(\sqrt{\nu})}{\sqrt{\nu}\\
K_1(\sqrt{\nu})}\right\].\$\$

The location component **redescends**: it grows like \\r/\sigma^2\\ near
the location and tends to \\\pm 1/\sigma\\ far from it, so a gross
outlier contributes a bounded amount to the estimating equation. That
bounded influence is the whole reason the pseudo-Huber loss exists, and
here it is the score of a genuine likelihood.

With `scale = "link"` the generic applies the chain rule for the links
the family carries. This method always returns the parameter scale. The
arithmetic runs in a compiled kernel; the Bessel ratio enters only the
\\\nu\\ component and is formed from the exponentially scaled functions.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

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

## Value

A named list of three numeric vectors, `mu`, `sigma` and `nu`, each of
length `max(length(y), length(mu), length(sigma), length(nu))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location,
\\\sigma \> 0\\ the scale, \\\nu \> 0\\ the shape, \\r = y - \mu\\ and
\\K_1\\ the modified Bessel function of the second kind of order one.

## See also

[`distrib_hessian.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PseudoHuberDistrib.md)
for the second derivatives,
[`distrib_expected_hessian.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PseudoHuberDistrib.md)
for their expectation,
[`distrib_grad_y.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.PseudoHuberDistrib.md)
for the derivative in the response, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 2)
g <- distrib_gradient(d, y, th)

# The location and scale components, written out.
r <- y - 0.4; D <- sqrt(2 + (r / 1.2)^2)
all.equal(g$mu, r / (1.2^2 * D))
#> [1] TRUE
all.equal(g$sigma, (r^2 / (1.2^2 * D) - 1) / 1.2)
#> [1] TRUE

# numDeriv on the summed log-density reproduces the summed score.
fn <- function(p)
  sum(distrib_pdf(d, y, list(mu = p[1], sigma = p[2], nu = p[3]), log = TRUE))
rbind(numeric = numDeriv::grad(fn, c(0.4, 1.2, 2)),
      closed = vapply(g, sum, numeric(1)))
#>                 mu     sigma          nu
#> numeric -0.2379499 -0.139081 0.003053376
#> closed  -0.2379499 -0.139081 0.003053376

# The location score is bounded by 1 / sigma, where a Gaussian's grows
# without bound.
rr <- c(1, 4, 16, 64)
rbind(residual = rr,
      pseudohuber = rr / (1.2^2 * sqrt(2 + (rr / 1.2)^2)),
      gaussian = rr / 1.2^2,
      bound = rep(1 / 1.2, 4))
#>                  [,1]      [,2]       [,3]       [,4]
#> residual    1.0000000 4.0000000 16.0000000 64.0000000
#> pseudohuber 0.4230609 0.7671455  0.8286850  0.8330405
#> gaussian    0.6944444 2.7777778 11.1111111 44.4444444
#> bound       0.8333333 0.8333333  0.8333333  0.8333333
```
