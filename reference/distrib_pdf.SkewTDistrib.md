# Skew t Density

Computes the skew \\t\\ density, with \\z = (y-\mu)/\sigma\\ and \\w =
\alpha z\sqrt{(\nu+1)/(\nu+z^2)}\\: \$\$f(y; \mu, \sigma, \alpha, \nu) =
\dfrac{2}{\sigma}\\ t\_\nu(z)\\T\_{\nu+1}(w),\$\$ with \\t\_\nu\\ the
standard Student \\t\\ density and \\T\_{\nu+1}\\ its distribution
function on **one more** degree of freedom. The extra degree of freedom
is deliberate: without it the tilted density would not integrate to one.

The two logarithms are taken separately and added, `dt(log = TRUE)`
beside `pt(log.p = TRUE)`, so the light tail of the skewed side returns
a large negative number instead of `-Inf`.

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- y:

  A numeric vector of observations, anywhere on the real line.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`, each a
  numeric vector of length 1 or of the length of `y`. `sigma` and `nu`
  must be strictly positive; `mu` and `alpha` may take any finite value.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of the length of the recycled inputs.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\alpha\\ the shape
and \\\nu \> 0\\ the degrees of freedom. Neither \\\mu\\ nor \\\sigma\\
is a moment.

## See also

[`distrib_gradient.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md)
for the score,
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
for the \\\nu \to \infty\\ limit,
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
for the \\\alpha = 0\\ case, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- skewt_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)

# The formula written out.
w <- 3 * y * sqrt(7 / (6 + y^2))
all.equal(distrib_pdf(d, y, th), 2 * dt(y, 6) * pt(w, 7))
#> [1] TRUE

# It integrates to one.
integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value
#> [1] 1

# Shape zero is the Student t.
all.equal(distrib_pdf(d, y, list(mu = 0, sigma = 1, alpha = 0, nu = 6)),
          dt(y, 6))
#> [1] TRUE

# Large degrees of freedom give the skew normal, at O(1/nu).
sn <- skewnormal1_distrib()
vapply(c(1e2, 1e4, 1e6), function(v)
  max(abs(distrib_pdf(d, y, list(mu = 0, sigma = 1, alpha = 3, nu = v)) -
          distrib_pdf(sn, y, list(mu = 0, sigma = 1, alpha = 3)))), 0)
#> [1] 2.425307e-03 2.435181e-05 2.435280e-07
```
