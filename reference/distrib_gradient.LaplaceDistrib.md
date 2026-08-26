# Laplace Score

Computes the first derivatives of the Laplace log-density with respect
to \\\mu\\ and \\\sigma\\, one value per observation, in closed form.
Writing \\r = y - \mu\\, \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{\mathrm{sign}(r)}{\sigma}, \qquad \dfrac{\partial \ell}{\partial
\sigma} = \dfrac{1}{\sigma}\left(\dfrac{\|r\|}{\sigma} - 1\right).\$\$

The score in \\\mu\\ takes only the three values \\-1/\sigma\\, 0 and
\\+1/\sigma\\: it carries the **sign** of the residual and nothing about
its size, which is why the maximum likelihood estimate of \\\mu\\ is the
sample median. At \\r = 0\\ exactly, `sign(0)` is 0 and the method
returns 0, the midpoint of the subdifferential \\\[-1/\sigma,
1/\sigma\]\\; the derivative does not exist there.

## Arguments

- distrib:

  A `LaplaceDistrib` object, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

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

## Value

A named list of two numeric vectors, `mu` and `sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location and
\\\sigma \> 0\\ the scale, with variance \\2\sigma^2\\. \\r = y - \mu\\
is the residual.

## See also

[`distrib_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LaplaceDistrib.md)
for the second derivatives, which vanish in \\\mu\\;
[`distrib_expected_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LaplaceDistrib.md)
for the information, which does not; and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- laplace_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out.
r <- y - 0.4
all.equal(g$mu, sign(r) / 1.5)
#> [1] TRUE
all.equal(g$sigma, (abs(r) / 1.5 - 1) / 1.5)
#> [1] TRUE

# The score in mu is a sign: three values, whatever the residual.
unique(distrib_gradient(d, 0.4 + c(-100, -1, 0, 1, 100), th)$mu)
#> [1] -0.6666667  0.0000000  0.6666667

# It vanishes summed at the sample median, which is the estimate of mu.
set.seed(12)
z <- distrib_rng(d, 401, list(mu = 3, sigma = 2))
sum(distrib_gradient(d, z, list(mu = median(z), sigma = 2))$mu)
#> [1] 0
```
