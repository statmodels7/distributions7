# Multivariate Student t Density

Computes the log-density \$\$\ell =
\log\Gamma\\\left(\tfrac{\nu+p}{2}\right) -
\log\Gamma\\\left(\tfrac{\nu}{2}\right) - \tfrac{p}{2}\log(\nu\pi) -
\tfrac{1}{2}\log\lvert\Sigma\rvert - \tfrac{\nu+p}{2}\log\\\left(1 +
\tfrac{q}{\nu}\right),\$\$ with \\q = (y-\mu)^\top \Sigma^{-1}(y-\mu)\\,
row by row, and exponentiates it unless `log = TRUE`. The last logarithm
is taken with [`base::log1p()`](https://rdrr.io/r/base/Log.html), which
is the difference between a number and a loss of every significant digit
when \\q/\nu\\ is small: near the center of a family with many degrees
of freedom, `log(1 + q/nu)` cancels its argument away.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations, one row each. A
  plain numeric vector of length \\p\\ is read as a single observation.
  Every point of \\\mathbb{R}^p\\ is in the support, so no row is
  rejected. A matrix with zero rows returns `numeric(0)`.

- theta:

  A named list of parameters, each component a single number: `mu1`,
  ..., `mup`, the matrix parametrization's prefixed free values, and
  `nu`, which must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-density is returned, which
  stays finite far into a tail. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length \\n\\, one density per row of `y`.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom, \\p\\ the dimension and \\q\\ the squared
Mahalanobis distance of an observation from the location.

## See also

[`distrib_gradient.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvStudentTDistrib.md)
for the derivatives,
[`distrib_rng.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.MvStudentTDistrib.md)
to draw from it,
[`distrib_pdf.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.MvGaussianDistrib.md)
for the limiting family, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- mvstudent_t_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
y <- rbind(c(0, 0), c(1, -1), c(0.5, -0.3))

distrib_pdf(d, y, theta)
#> [1] 0.12386429 0.07758092 0.17589341

# Against the formula written out.
S <- mv_sigma(d, theta)
q <- mahalanobis(y, c(0.5, -0.3), S)
ref <- lgamma(4) - lgamma(3) - log(6 * pi) - 0.5 * log(det(S)) -
  4 * log1p(q / 6)
all.equal(distrib_pdf(d, y, theta, log = TRUE), ref)
#> [1] TRUE

# As nu grows the density approaches the gaussian's at the same matrix.
g <- mvgaussian_distrib(2)
vapply(c(6, 60, 6000), function(nu) {
  t2 <- theta; t2$nu <- nu
  distrib_pdf(d, c(1, -1), t2, log = TRUE)
}, numeric(1))
#> [1] -2.556434 -2.433958 -2.419273
distrib_pdf(g, c(1, -1), theta[1:5], log = TRUE)
#> [1] -2.419124
```
