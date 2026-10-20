# The Closed-Form Higher Derivatives of a Multivariate Gaussian

Computes every third- or fourth-order derivative of the log-density in
the parameters. It enumerates the parameter tuples the way
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
does, splits each into mean indices and matrix indices, and reads the
surviving cases off the gaussian's algebra. Writing \\r = y - \mu\\ and
\\P_t\\ for the precision's derivative array over the matrix indices
\\t\\, the cases are as follows.

- A tuple with three or more mean indices is exactly zero, because the
  quadratic form is quadratic in \\\mu\\.

- A tuple with two mean indices \\(i, j)\\ gives \\-P_t\[i, j\]\\.

- A tuple with one mean index \\i\\ gives \\(r P_t)\_i\\, one value per
  observation.

- A tuple without mean indices gives \\\pm\tfrac{1}{2}\\\partial^{\lvert
  t\rvert}\log\lvert M\rvert - \tfrac{1}{2}\\ r^\top P_t\\ r\\, where
  the sign depends on which side the parametrization carries.

## Usage

``` r
mvg_higher(distrib, y, theta, order, expected = FALSE)
```

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters in which each component is a single number.

- order:

  The derivative order, `3L` or `4L`.

- expected:

  Logical of length 1. When `TRUE` the expectation of each component is
  returned. Defaults to `FALSE`.

## Value

A named list of numeric vectors of length \\n\\, one per derivative
component, keyed and ordered as `deriv_names(distrib@params, order)`.
There are 35 components at order 3 and 70 at order 4 for a
two-dimensional gaussian on an unstructured covariance.

## Details

The log-determinant derivatives come from
[`parameters7::param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.html)
and `param_d4logdet()`, and the precision arrays from
[`mvg_ptensors()`](https://statmodels7.github.io/distributions7/reference/mvg_ptensors.md).
With `expected = TRUE` the response enters through \\r\\ and \\r^\top
P_t r\\, whose expectations are 0 and \\\operatorname{tr}(P_t\Sigma)\\,
so the expectation is exact.

## Notation

\\\mu\\ is the mean, \\M\\ the matrix that the parametrization carries,
\\P = \Sigma^{-1}\\ the precision, \\\eta\\ the free vector, \\r = y -
\mu\\ the centered response and \\P_t\\ the precision's derivative array
over the multiset of free values \\t\\.

## See also

[`distrib_deriv3.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MvGaussianDistrib.md)
and
[`distrib_deriv4.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.MvGaussianDistrib.md),
the two methods that it serves, and
[`mvg_ptensors()`](https://statmodels7.github.io/distributions7/reference/mvg_ptensors.md)
for the arrays that it reads.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
set.seed(1)
y <- distrib_rng(d, 4, theta)

d3 <- distributions7:::mvg_higher(d, y, theta, 3L)
length(d3)
#> [1] 35
length(distributions7:::mvg_higher(d, y, theta, 4L))
#> [1] 70

# Three mean indices vanish exactly.
vapply(d3[grep("^mu[0-9]+_mu[0-9]+_mu[0-9]+$", names(d3))],
       function(z) max(abs(z)), numeric(1))
#> mu1_mu1_mu1 mu1_mu1_mu2 mu1_mu2_mu2 mu2_mu2_mu2 
#>           0           0           0           0 

# The result is compared with one stencil on the analytic Hessian, which
# shares no algebra with the expansion.
h <- 1e-4
tp <- theta; tp$sigma_L2.1 <- tp$sigma_L2.1 + h
tm <- theta; tm$sigma_L2.1 <- tm$sigma_L2.1 - h
c(exact = sum(d3[["sigma_log_L1_sigma_log_L1_sigma_L2.1"]]),
  stencil = (sum(distrib_hessian(d, y, tp)[["sigma_log_L1_sigma_log_L1"]]) -
             sum(distrib_hessian(d, y, tm)[["sigma_log_L1_sigma_log_L1"]])) /
            (2 * h))
#>     exact   stencil 
#> -6.063894 -6.063894 
```
