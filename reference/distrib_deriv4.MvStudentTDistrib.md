# Multivariate Student t Fourth Derivatives

Computes every fourth derivative of the log-density in the parameters,
in closed form, by the same construction as
[`distrib_deriv3.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MvStudentTDistrib.md)
one order up: the location and matrix part through
[`mvg_ptensors()`](https://statmodels7.github.io/distributions7/reference/mvg_ptensors.md)'s
array expansion, and the \\\nu\\ part by differentiating
\\(-1)^{m-1}(m-1)!/(\nu+q)^m\\ against the linear prefactor
\\(\nu+p)/2\\.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters, each component a single number.

- expected:

  Logical of length 1. When `TRUE` the expectation is returned, by
  SAMPLING: `nsim` draws are made and the observed components averaged.
  Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. The two differ here, `nu` carrying a log link by
  default.

- approx:

  Ignored: sampling is the only multivariate route to an expectation.
  Present so that the signature matches the generic's.

- nsim:

  The number of draws used when `expected = TRUE`. Defaults to `10000`.
  Ignored otherwise.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed and ordered as
`deriv_names(distrib@params, 4)`. At \\p = 2\\ on an unstructured scale
matrix there are 126 components. With `expected = TRUE` every vector is
constant.

## Details

The license for this order is that the SAME assembly run at orders one
and two reproduces the hand-written score and information, which are
derived separately and are already under
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md).

With `expected = TRUE` the expectation is taken by sampling and carries
Monte Carlo error of order `nsim^(-1/2)`.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom, \\q\\ the squared Mahalanobis distance, \\p\\ the
dimension and \\\ell\\ the log-density of one observation.

## See also

[`distrib_deriv3.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MvStudentTDistrib.md)
for the order below,
[`mvt_higher()`](https://statmodels7.github.io/distributions7/reference/mvt_higher.md)
for the shared engine, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- mvstudent_t_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(3)
y <- distrib_rng(d, 4, theta)

d4 <- distrib_deriv4(d, y, theta)
length(d4)
#> [1] 126

# Against one stencil on the analytic third order.
h <- 1e-4
tp <- theta; tp$nu <- tp$nu + h
tm <- theta; tm$nu <- tm$nu - h
c(exact = sum(d4[["mu1_mu2_nu_nu"]]),
  stencil = (sum(distrib_deriv3(d, y, tp)[["mu1_mu2_nu"]]) -
             sum(distrib_deriv3(d, y, tm)[["mu1_mu2_nu"]])) / (2 * h))
#>        exact      stencil 
#> 0.0009689403 0.0009689403 
```
