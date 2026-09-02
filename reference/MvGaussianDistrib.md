# Multivariate Gaussian Distribution Class

The S7 class of the gaussian family on \\\mathbb{R}^p\\, with density
\$\$f(y) = (2\pi)^{-p/2}\lvert\Sigma\rvert^{-1/2}
\exp\\\left\\-\tfrac{1}{2}(y-\mu)^\top\Sigma^{-1}(y-\mu)\right\\.\$\$
The mean \\\mu\\ contributes \\p\\ scalar parameters and the matrix is
carried by a parameters7 parametrization, whose free values become
scalar parameters in their turn. It inherits from
`multivariate_distrib`, so the response is an \\n \times p\\ matrix and
the distribution function and the quantile function are refused.

Build one with
[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md),
which fills the properties in and checks that the matrix parametrization
has full rank. This page documents the raw S7 constructor, which
validates neither the rank nor the agreement between `n_dim` and the
parametrization's dimension.

## Usage

``` r
MvGaussianDistrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  n_dim = integer(0),
  param = parameters7::parameter(),
  inverted = logical(0)
)
```

## Arguments

- distrib_name:

  A single character string specifying the name of the distribution
  (e.g., `"student t"`).

- dimension:

  A character string indicating the dimensionality (`"univariate"` or
  `"multivariate"`).

- bounds:

  A numeric vector of length 2 defining the overall support of the
  distribution `c(lower, upper)`.

- params:

  A character vector containing the names of the distribution parameters
  (e.g., `c("mu", "sigma")`).

- params_interpretation:

  A character vector (typically named) providing the statistical
  interpretation of each parameter (e.g., `c(mu = "location")`).

- n_params:

  A numeric value specifying the total number of parameters.

- params_bounds:

  A list of numeric vectors of length 2, specifying the valid
  mathematical domain for each individual parameter.

- link_params:

  A list of link function objects corresponding to each parameter,
  primarily used to map parameters to the unconstrained real line for
  optimization algorithms.

- params_smooth:

  An optional named logical vector flagging, for each parameter, whether
  the log-likelihood is differentiable with respect to it. Defaults to
  all `TRUE` (leave empty). Set an entry to `FALSE` for parameters at
  which the log-likelihood has a kink (e.g. the location of a Laplace
  distribution): the observed Hessian is then degenerate and the
  expected information must be obtained from the score variance rather
  than from \\-\mathbb{E}\[H\]\\ (see
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)).

- n_dim:

  The dimension \\p\\ of one observation. A single positive integer.

- param:

  A parameters7 parametrization of the matrix, inheriting from
  [`parameters7::parameter`](https://statmodels7.github.io/parameters7/reference/parameter.html).
  Its `n_free` free values are flattened into scalar parameters of the
  distribution.

- inverted:

  Logical of length 1. `TRUE` when `param` carries the precision
  \\\Omega = \Sigma^{-1}\\ and `FALSE` when it carries the covariance
  \\\Sigma\\. Nothing but the sign of the log-determinant and one matrix
  inversion depends on it; the law is the same either way.

## Value

An S7 object of class `MvGaussianDistrib`, inheriting from
`multivariate_distrib` and from `distrib`. Beyond the parent's
`distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params`,
`params_smooth` and `n_dim`, it carries `param` and `inverted` as
supplied. For an object built by
[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
at \\p = 2\\ on an unstructured covariance, `params` is
`c("mu1", "mu2", "sigma_log_L1", "sigma_log_L2", "sigma_L2.1")`, every
`params_bounds` entry is \\(-\infty, \infty)\\ and every link is the
identity.

## Methods

Registered on this class:
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.MvGaussianDistrib.md),
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.MvGaussianDistrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MvGaussianDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.MvGaussianDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvGaussianDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.MvGaussianDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvGaussianDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvGaussianDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvGaussianDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.MvGaussianDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.MvGaussianDistrib.md),
[`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.MvGaussianDistrib.md),
[`mean()`](https://statmodels7.github.io/distributions7/reference/mean.MvGaussianDistrib.md),
[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.MvGaussianDistrib.md),
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvGaussianDistrib.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.MvGaussianDistrib.md)

Everything else comes from
[`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md),
including the two refusals.

## See also

[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
to build one,
[`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
for the heavy-tailed sibling,
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
for the covariance a parameter vector describes, and
[`distrib_gradient.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvGaussianDistrib.md)
for the closed-form score.

## Examples

``` r
d <- mvgaussian1_distrib(2)
S7::S7_inherits(d, multivariate_distrib)
#> [1] TRUE

# Five scalar parameters: two means and the three free values of the
# log-Cholesky covariance, prefixed by the matrix they describe.
d@params
#> [1] "mu1"          "mu2"          "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
d@n_dim
#> [1] 2
d@param@free_names
#> [1] "log_L1" "log_L2" "L2.1"  

# Every parameter is already unconstrained, so every link is the identity
# and the positive definiteness lives inside the matrix parametrization.
vapply(d@link_params, function(l) l@link_name, character(1))
#>          mu1          mu2 sigma_log_L1 sigma_log_L2   sigma_L2.1 
#>   "identity"   "identity"   "identity"   "identity"   "identity" 

# The precision form differs in one property and in the prefix.
o <- mvgaussian2_distrib(2, parameters7::log_cholesky(2))
c(covariance = d@inverted, precision = o@inverted)
#> covariance  precision 
#>      FALSE       TRUE 
o@params
#> [1] "mu1"          "mu2"          "omega_log_L1" "omega_log_L2" "omega_L2.1"  
```
