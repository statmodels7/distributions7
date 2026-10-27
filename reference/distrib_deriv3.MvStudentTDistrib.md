# Multivariate Student t Third Derivatives

Computes every third derivative of the log-density in the parameters, in
closed form. The log-density carries `lgamma`, a logarithm and a
quadratic form, and each of them is elementary in \\\nu\\. The
location-and-matrix part of a component is handled by the array
expansion that the gaussian uses through
[`mvg_ptensors()`](https://statmodels7.github.io/distributions7/reference/mvg_ptensors.md),
and the \\\nu\\ part by differentiating \\(-1)^{m-1}(m-1)!/(\nu+q)^m\\
against the linear prefactor \\(\nu+p)/2\\. The components in \\\nu\\
alone come from
[`mvt_cp_d()`](https://statmodels7.github.io/distributions7/reference/mvt_cp_d.md)
and
[`mvt_h_d()`](https://statmodels7.github.io/distributions7/reference/mvt_h_d.md).

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
  or
  [`mvstudent_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters, where each component is a single number.

- expected:

  Logical of length 1. When `TRUE` the method returns the expectation,
  which it takes by sampling: it makes the number of draws that `method`
  carries and averages the observed components. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`. The method returns the
  parameter scale and the generic applies the link scale to its result.
  The two differ here, because `nu` carries a log link by default.

- method:

  `NULL`: the generic applies an object before dispatching on the
  family, so this method receives `NULL`. With `expected = TRUE` the
  method takes the expectation by sampling 10000 draws; an object such
  as
  [`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
  with another `nsim` is applied by the generic. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed and ordered as
`deriv_names(distrib@params, 3)`. At \\p = 2\\ on an unstructured scale
matrix there are 56 components. With `expected = TRUE` every vector is
constant.

## Details

Run at orders one and two, the same assembly reproduces
[`distrib_gradient.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvStudentTDistrib.md)
and
[`distrib_hessian.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvStudentTDistrib.md),
which are written separately, for the scale matrix and for its inverse.

For the gaussian, every component with three location indices is zero,
because the quadratic form is quadratic in \\\mu\\. Here the weight
depends on \\q\\, and \\\ell^{(\mu_1\mu_1\mu_1)}\\ is an ordinary
non-zero number.

With `expected = TRUE` the expectation is taken by sampling and carries
Monte Carlo error that falls as the inverse square root of the number of
draws.

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

# The result is compared with a central difference of the Hessian in nu.
h <- 1e-4
tp <- theta; tp$nu <- tp$nu + h
tm <- theta; tm$nu <- tm$nu - h
c(exact = sum(d3[["mu1_mu2_nu"]]),
  stencil = (sum(distrib_hessian(d, y, tp)[["mu1_mu2"]]) -
             sum(distrib_hessian(d, y, tm)[["mu1_mu2"]])) / (2 * h))
#>     exact   stencil 
#> 0.0165313 0.0165313 
```
