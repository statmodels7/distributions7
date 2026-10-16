# Gamma Third-Order Derivatives in Mean and Variance

Computes the four distinct third derivatives of the gamma log-density
with respect to \\\mu\\ and \\\sigma^2\\, in closed form. With \\\alpha
= \mu^2/\sigma^2\\ and \\z = y/\mu\\ the log-density is \\\ell =
G(\alpha) + \alpha e(\mu) - \log y\\, with \\G(\alpha) =
\alpha\log\alpha - \log\Gamma(\alpha) - \alpha\\ and \\e(\mu) = \log z -
(z - 1)\\. The derivatives of \\G(\alpha(\mu, \sigma^2))\\ come from the
formula of Faa di Bruno over the set partitions of the differentiation
indices, and those of \\\alpha e(\mu)\\ from the product rule. The
derivatives of \\G\\ are \\\log\alpha - \psi(\alpha)\\, \\1/\alpha -
\psi_1(\alpha)\\, \\-1/\alpha^2 - \psi_2(\alpha)\\ and \\2/\alpha^3 -
\psi_3(\alpha)\\. The kernel evaluates each of them from its asymptotic
series at a large shape, and \\e\\ from a series near \\z = 1\\, so the
components keep their digits when \\\alpha\\ is large. Every component
of this order carries \\\psi_2\\, the second derivative of the digamma
function.

With `expected = TRUE` the expectations are returned, obtained by
replacing \\Y\\ with \\\mu\\ and \\\log Y\\ with \\\psi(\alpha) -
\log\lambda\\. Both routes are closed form, so the method uses no
quadrature and ignores `approx` and `nsim`.

## Arguments

- distrib:

  A `Gamma2Distrib` object, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- y:

  A numeric vector of strictly positive observations. With
  `expected = TRUE` only its length is used.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the expectations under the model are
  returned in place of the observed values. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). The
  generic reads it, and this method does not.

- approx:

  Ignored, because a closed form is available for both the observed and
  the expected values.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  a fixed number of observations, set per kernel, the kernel runs
  sequentially. Defaults to `1L`.

## Value

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma2`,
`mu_sigma2_sigma2` and `sigma2_sigma2_sigma2`, each of length
`max(length(y), length(mu), length(sigma2))`. The names enumerate the
distinct multi-indices of order three in two parameters.

## Notation

\\\ell^{(ijk)}\\ is the third derivative of the log-density in
parameters \\i\\, \\j\\ and \\k\\. \\\psi\\ is the digamma function and
\\\psi_m\\ its \\m\\th derivative.

## See also

[`distrib_hessian.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gamma2Distrib.md)
for the order below and
[`distrib_deriv4.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gamma2Distrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic and for the numerical route that a family without a
closed form uses.

## Examples

``` r
d <- gamma2_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, sigma2 = 2)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "mu_mu_mu"             "mu_mu_sigma2"         "mu_sigma2_sigma2"    
#> [4] "sigma2_sigma2_sigma2"

# Every component varies with the observation.
d3$mu_mu_mu
#> [1] -0.2431584 -0.2431584 -0.2431584

# The expected values are constants at a fixed parameter setting.
lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu
#> [1] -0.2431584
#> 
#> $mu_mu_sigma2
#> [1] 0.3016318
#> 
#> $mu_sigma2_sigma2
#> [1] -0.1728947
#> 
#> $sigma2_sigma2_sigma2
#> [1] 0.2638418
#> 

# A central difference of the Hessian reproduces the same component.
eps <- 1e-6
up <- distrib_hessian(d, y, list(mu = 3 + eps, sigma2 = 2))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 3 - eps, sigma2 = 2))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-4)
#> [1] TRUE
```
