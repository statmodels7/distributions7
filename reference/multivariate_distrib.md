# S7 Class for Multivariate Distributions

Defines the subclass of `distrib` for distributions whose observations
are vectors. The response is an \\n \times p\\ matrix, one row per
observation, and `n_dim` records \\p\\. The class is a sibling of
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
and
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md)
and not a subclass of either, so it inherits none of their
one-dimensional defaults. The methods registered on this class signal an
error for the distribution function, the quantile function and the
response derivatives; the gaussian and the Student t register response
derivatives of their own.

## Usage

``` r
multivariate_distrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  n_dim = integer(0)
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

  The dimension \\p\\ of one observation, a single positive integer. The
  validator rejects any other value, and also a `dimension` property
  other than `"multivariate"`.

## Value

An S7 object of class `multivariate_distrib`, inheriting from `distrib`.
Beyond the parent's `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`, it carries `n_dim`.

## The parameters are scalars

The parameters of a multivariate distribution are scalars, as for every
other distribution of the package. A mean vector contributes \\p\\ of
them and a matrix contributes the free values of the parameters7
parametrization that carries it, so `theta` is a named list of numbers
and
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
the Hessian keys, the link scale and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
apply without a special case.

The constraint on a matrix or on a probability vector is carried by its
parametrization, so the free values of a matrix and of a simplex take
the identity link. A scalar parameter keeps a link of its own: the
degrees of freedom of the Student t and the concentration of the
Dirichlet take the log link by default.

## The one-dimensional defaults

The defaults registered on
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
are one-dimensional: a distribution function by quadrature over an
interval, a quantile by root finding and a generator for a scalar
density. They have no counterpart in \\p\\ dimensions, so the class does
not inherit them.

## Parameters that vary by observation

A parametrization describes one matrix for the whole sample, and
[`mv_flat_theta()`](https://statmodels7.github.io/distributions7/reference/mv_flat_theta.md)
rejects a `theta` component longer than one. A distribution whose
parameters depend on covariates is a model, which belongs to the layer
above this package.

## See also

[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
and
[`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
for the elliptical families,
[`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
and
[`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
for the simplex-valued ones,
[`n_obs()`](https://statmodels7.github.io/distributions7/reference/n_obs.md)
for the observation count, and
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for the quantities that a fit reports.

## Examples

``` r
d <- mvgaussian1_distrib(2)
S7::S7_inherits(d, multivariate_distrib)
#> [1] TRUE
c(n_dim = d@n_dim, n_params = d@n_params)
#>    n_dim n_params 
#>        2        5 

# The class is a sibling of continuous_distrib and not a subclass.
c(multivariate = S7::S7_inherits(d, multivariate_distrib),
  continuous = S7::S7_inherits(d, continuous_distrib))
#> multivariate   continuous 
#>         TRUE        FALSE 

# The free values of the matrix take the identity link; the degrees of
# freedom of a Student t take the log link.
vapply(d@link_params, function(l) l@link_name, character(1))
#>          mu1          mu2 sigma_log_L1 sigma_log_L2   sigma_L2.1 
#>   "identity"   "identity"   "identity"   "identity"   "identity" 
vapply(mvstudent_t1_distrib(2)@link_params, function(l) l@link_name,
       character(1))
#>          mu1          mu2 sigma_log_L1 sigma_log_L2   sigma_L2.1           nu 
#>   "identity"   "identity"   "identity"   "identity"   "identity"        "log" 

# The shipped multivariate families all inherit from the class.
vapply(list(mvgaussian1_distrib(2), mvstudent_t1_distrib(2),
            dirichlet_distrib(3), multinomial_distrib(3, size = 5)),
       function(x) S7::S7_inherits(x, multivariate_distrib), TRUE)
#> [1] TRUE TRUE TRUE TRUE
```
