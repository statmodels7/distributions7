# Multivariate Student t Distribution Class

The S7 class of the elliptical Student t family on \\\mathbb{R}^p\\,
with density \$\$f(y) \propto \lvert\Sigma\rvert^{-1/2} \left(1 +
\frac{(y-\mu)^\top \Sigma^{-1} (y-\mu)}{\nu}\right)^{-(\nu+p)/2}.\$\$
The location \\\mu\\ contributes \\p\\ scalar parameters, the scale
matrix \\\Sigma\\ is carried by a parameters7 parametrization whose free
values become scalar parameters, and the degrees of freedom \\\nu\\ is
added last. It inherits from `multivariate_distrib`, so the response is
an \\n \times p\\ matrix and the distribution function and the quantile
function are refused.

Build one with
[`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md),
which fills the properties in and checks the rank of the
parametrization. This page documents the raw S7 constructor, which
validates nothing.

## Usage

``` r
MvStudentTDistrib(
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

  A parameters7 parametrization of the scale matrix or of its inverse,
  according to `inverted`, inheriting from
  [`parameters7::parameter`](https://statmodels7.github.io/parameters7/reference/parameter.html).
  Its `n_free` free values are flattened into scalar parameters of the
  distribution.

- inverted:

  Logical of length 1. `TRUE` when `param` carries the inverse scale
  matrix \\\Sigma^{-1}\\ and `FALSE` when it carries the scale matrix
  \\\Sigma\\. Nothing but the sign of the log-determinant and one matrix
  inversion depends on it; the law is the same either way, and neither
  matrix is a moment of the response.

## Value

An S7 object of class `MvStudentTDistrib`, inheriting from
`multivariate_distrib` and from `distrib`. Beyond the parent's
`distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params`,
`params_smooth` and `n_dim`, it carries `param` and `inverted` as
supplied.

## Methods

Registered on this class:
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.MvStudentTDistrib.md),
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.MvStudentTDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvStudentTDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.MvStudentTDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvStudentTDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvStudentTDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvStudentTDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.MvStudentTDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.MvStudentTDistrib.md),
[`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.MvStudentTDistrib.md),
[`mean()`](https://statmodels7.github.io/distributions7/reference/mean.MvStudentTDistrib.md),
[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.MvStudentTDistrib.md),
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MvStudentTDistrib.md),
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvStudentTDistrib.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.MvStudentTDistrib.md)

Everything else comes from
[`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).
The third and fourth derivatives are registered in `mv_higher.R`, which
the two multivariate families share.

## See also

[`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
to build one,
[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
for the limit \\\nu \to \infty\\,
[`mv_sigma.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvStudentTDistrib.md)
for the scale matrix and
[`variance.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.MvStudentTDistrib.md)
for the covariance, which are different matrices here.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
S7::S7_inherits(d, multivariate_distrib)
#> [1] TRUE

# Six parameters: two locations, three free values of the scale matrix, and
# the degrees of freedom last.
d@params
#> [1] "mu1"          "mu2"          "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
#> [6] "nu"          
d@params_interpretation
#>                  mu1                  mu2         sigma_log_L1 
#>           "location"           "location"              "scale" 
#>         sigma_log_L2           sigma_L2.1                   nu 
#>              "scale"              "scale" "degrees of freedom" 

# nu is the one parameter with a bound, and the one with a link that is not
# the identity, so this family's link scale is not its parameter scale.
d@params_bounds$nu
#> [1]   0 Inf
vapply(d@link_params, function(l) l@link_name, character(1))
#>          mu1          mu2 sigma_log_L1 sigma_log_L2   sigma_L2.1           nu 
#>   "identity"   "identity"   "identity"   "identity"   "identity"        "log" 
```
