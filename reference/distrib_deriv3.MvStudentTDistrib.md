# Multivariate Student t Third Derivatives

Computes every third derivative of the log-density in the parameters, in
closed form. Nothing is obstructed here: the log-density carries
`lgamma`, a logarithm and a quadratic form, each elementary in \\\nu\\,
so there is no distribution function to differentiate in its degrees of
freedom. The component splits into a location-and-matrix part, handled
by the same array expansion the gaussian uses through
[`mvg_ptensors()`](https://statmodels7.github.io/distributions7/reference/mvg_ptensors.md),
and a \\\nu\\ part obtained by differentiating
\\(-1)^{m-1}(m-1)!/(\nu+q)^m\\ against the linear prefactor
\\(\nu+p)/2\\.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

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
`deriv_names(distrib@params, 3)`. At \\p = 2\\ on an unstructured scale
matrix there are 56 components. With `expected = TRUE` every vector is
constant.

## Details

The license for this order is that the SAME assembly run at orders one
and two reproduces the hand-written score and information, which are
derived separately and are already under
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md).
Against one stencil on the analytic Hessian, every component here agrees
to a relative \\1.2\times10^{-8}\\ or better.

Unlike the gaussian's, no component vanishes. There a tuple with three
location indices is exactly zero, the quadratic form being quadratic in
\\\mu\\; here the weight depends on \\q\\, so
\\\ell^{(\mu_1\mu_1\mu_1)}\\ is an ordinary non-zero number.

With `expected = TRUE` the expectation is taken by sampling and carries
Monte Carlo error of order `nsim^(-1/2)`.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom, \\q\\ the squared Mahalanobis distance, \\p\\ the
dimension and \\\ell\\ the log-density of one observation.

## See also

[`distrib_deriv4.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.MvStudentTDistrib.md)
for the next order,
[`distrib_hessian.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvStudentTDistrib.md)
for the second,
[`mvt_higher()`](https://statmodels7.github.io/distributions7/reference/mvt_higher.md)
for the shared engine,
[`distrib_deriv3.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MvGaussianDistrib.md)
for the limiting family, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(3)
y <- distrib_rng(d, 4, theta)

d3 <- distrib_deriv3(d, y, theta)
length(d3)
#> [1] 56

# A gaussian's three-location component is exactly zero and this one is not.
c(t = sum(d3[["mu1_mu1_mu1"]]),
  gaussian = sum(distrib_deriv3(mvgaussian1_distrib(2), y,
                                theta[1:5])[["mu1_mu1_mu1"]]))
#>        t gaussian 
#>  1.36229  0.00000 

# Against one stencil on the analytic Hessian.
h <- 1e-4
tp <- theta; tp$nu <- tp$nu + h
tm <- theta; tm$nu <- tm$nu - h
c(exact = sum(d3[["mu1_mu2_nu"]]),
  stencil = (sum(distrib_hessian(d, y, tp)[["mu1_mu2"]]) -
             sum(distrib_hessian(d, y, tm)[["mu1_mu2"]])) / (2 * h))
#>     exact   stencil 
#> 0.0165313 0.0165313 
```
