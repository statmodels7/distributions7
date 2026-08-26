# Gamma Distribution Class, Mean and Dispersion

The S7 class of the gamma family parametrized by its mean \\\mu \> 0\\
and a dispersion \\\phi \> 0\\, so that \\\operatorname{Var}(Y) =
\phi\mu^2\\. The shape is \\1/\phi\\ and the rate \\1/(\phi\mu)\\. It
inherits from `continuous_distrib`, so it answers every generic of the
`distrib` contract; the eleven methods listed below are registered on it
directly and everything else comes from the parent.

Build one with
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md),
which supplies the two link functions and fills the properties in. This
page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
Gamma1Distrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0)
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

## Value

An S7 object of class `Gamma1Distrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
they hold `"gamma1"`, `"univariate"`, `c(0, Inf)`, `c("mu", "phi")`, the
interpretations `c(mu = "mean", phi = "dispersion")`, `2`, the domain
\\(0, \infty)\\ for both parameters, and the two links.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gamma1Distrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gamma1Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gamma1Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gamma1Distrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Gamma1Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gamma1Distrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Gamma1Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gamma1Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gamma1Distrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Gamma1Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Gamma1Distrib.md)

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
to build one;
[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md)
for the same law in mean and variance;
[`distrib_pdf.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gamma1Distrib.md)
and
[`distrib_gradient.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gamma1Distrib.md)
for the closed forms this class supplies.

## Examples

``` r
d <- gamma1_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The properties a consumer reads to drive the family without knowing it.
d@params
#> [1] "mu"  "phi"
d@params_interpretation
#>           mu          phi 
#>       "mean" "dispersion" 
d@bounds
#> [1]   0 Inf

# The dispersion multiplies the variance function V(mu) = mu^2, so the
# coefficient of variation is sqrt(phi) whatever the mean.
sqrt(variance(d, list(mu = 3, phi = 0.5))) / 3
#> [1] 0.7071068
sqrt(variance(d, list(mu = 300, phi = 0.5))) / 300
#> [1] 0.7071068
```
