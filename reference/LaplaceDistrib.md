# Laplace Distribution Class, Location and Scale

This is the S7 class of the Laplace (double exponential) family with
location \\\mu\\ and scale \\\sigma \> 0\\, with density \\f(y) =
(2\sigma)^{-1}\exp(-\|y-\mu\|/\sigma)\\ on the whole real line. It
inherits from `continuous_distrib`, so every generic of the `distrib`
class applies to it; the methods listed below are registered on it
directly and everything else comes from the parent.

The density has a kink at \\y = \mu\\, where \\\|y - \mu\|\\ is not
differentiable. This makes the family non-regular in its location: the
observed second derivative in \\\mu\\ is 0 almost everywhere while the
information is \\1/\sigma^2\\. The class records this by setting
`params_smooth` to `c(mu = FALSE, sigma = TRUE)`, a property that
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
and the finite-difference guards read.

Build one with
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).
This page documents the raw S7 constructor, which takes the parent's
properties; the validator of `distrib` checks them, and the constructor
adds no check of its own.

## Usage

``` r
LaplaceDistrib(
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

An S7 object of class `LaplaceDistrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
they hold `"laplace"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "sigma")`, the interpretations
`c(mu = "location", sigma = "scale")`, `2`, the domains \\(-\infty,
\infty)\\ and \\(0, \infty)\\, the two links, and
`c(mu = FALSE, sigma = TRUE)` for `params_smooth`.

## Methods

The following methods are registered on this class, and each is linked
to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.LaplaceDistrib.md),
  [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.LaplaceDistrib.md),
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.LaplaceDistrib.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.LaplaceDistrib.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LaplaceDistrib.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LaplaceDistrib.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LaplaceDistrib.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LaplaceDistrib.md),
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.LaplaceDistrib.md)
  and
  [`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.LaplaceDistrib.md);

- the derivatives of the expected information:
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.elementary.md)
  and
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.elementary.md);

- the derivatives in the response:
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.LaplaceDistrib.md),
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.LaplaceDistrib.md),
  [`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md)
  and
  [`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md);

- the mixed derivatives in the response and the parameters:
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.LaplaceDistrib.md),
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross2_y.md),
  [`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/loc_cross3_y.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_y_hess.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_y_hess.md);

- the derivatives of the distribution function:
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.LaplaceDistrib.md),
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.LaplaceDistrib.md),
  [`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_deriv_cdf_k.md)
  and
  [`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_deriv_cdf_k.md);

- the moments:
  [`mean()`](https://statmodels7.github.io/distributions7/reference/mean.LaplaceDistrib.md),
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.LaplaceDistrib.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.LaplaceDistrib.md)
  and
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LaplaceDistrib.md);

- the decomposition at the kink:
  [`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.LaplaceDistrib.md).

Every other generic is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
to build one;
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
for the same law written by its rate \\\lambda = 1/\sigma\\;
[`distrib_expected_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LaplaceDistrib.md)
for the information at the kink;
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
for a smooth relative.

## Examples

``` r
d <- laplace_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The class declares that the log-likelihood has a kink in mu.
d@params_smooth
#>    mu sigma 
#> FALSE  TRUE 

# sigma is a scale, not a standard deviation: the variance is 2 sigma^2.
th <- list(mu = 0.4, sigma = 1.5)
c(variance = variance(d, th), two_sigma_sq = 2 * 1.5^2)
#>     variance two_sigma_sq 
#>          4.5          4.5 

# The tails are heavier than a Gaussian's: the excess kurtosis is 3.
kurtosis(d, th)
#> [1] 3
```
