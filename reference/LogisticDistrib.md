# Logistic Distribution Class

This is the S7 class of the logistic family with mean \\\mu\\ and scale
\\\sigma \> 0\\, whose distribution function is the logistic sigmoid
\\F(q) = \[1 + e^{-(q-\mu)/\sigma}\]^{-1}\\ on the whole real line. It
inherits from `continuous_distrib`, so every generic of the `distrib`
class applies to it; the methods listed below are registered on it
directly and everything else comes from the parent.

The family is symmetric about \\\mu\\, which is therefore both the mean
and the median. Its variance is \\\pi^2\sigma^2/3\\, so \\\sigma\\ is
not the standard deviation, and its excess kurtosis is \\6/5\\, slightly
heavier than a Gaussian's.

Build one with
[`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).
This page documents the raw S7 constructor, which takes the parent's
properties; the validator of `distrib` checks them, and the constructor
adds no check of its own.

## Usage

``` r
LogisticDistrib(
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

An S7 object of class `LogisticDistrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
they hold `"logistic"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "sigma")`, the interpretations
`c(mu = "mean", sigma = "scale")`, `2`, the domains \\(-\infty,
\infty)\\ and \\(0, \infty)\\, and the two links.

## Methods

The following methods are registered on this class, and each is linked
to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.LogisticDistrib.md),
  [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.LogisticDistrib.md),
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.LogisticDistrib.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.LogisticDistrib.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LogisticDistrib.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LogisticDistrib.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LogisticDistrib.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LogisticDistrib.md),
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.LogisticDistrib.md)
  and
  [`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.LogisticDistrib.md);

- the derivatives of the expected information:
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.elementary.md)
  and
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.elementary.md);

- the derivatives in the response:
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.LogisticDistrib.md),
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.LogisticDistrib.md),
  [`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md)
  and
  [`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md);

- the mixed derivatives in the response and the parameters:
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.LogisticDistrib.md),
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross2_y.md),
  [`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/loc_cross3_y.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_y_hess.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_y_hess.md);

- the derivatives of the distribution function:
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.LogisticDistrib.md),
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.LogisticDistrib.md),
  [`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_deriv_cdf_k.md)
  and
  [`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_deriv_cdf_k.md);

- the moments:
  [`mean()`](https://statmodels7.github.io/distributions7/reference/mean.LogisticDistrib.md),
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.LogisticDistrib.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.LogisticDistrib.md)
  and
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LogisticDistrib.md).

Every other generic is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
to build one;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md),
which it resembles with slightly heavier tails;
[`distrib_expected_hessian.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LogisticDistrib.md)
for the information.

## Examples

``` r
d <- logistic_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE
d@params_interpretation
#>      mu   sigma 
#>  "mean" "scale" 

# sigma is a scale, not a standard deviation: the variance is pi^2 sigma^2/3.
th <- list(mu = 0.4, sigma = 1.5)
c(variance = variance(d, th), pi_sq_over_3 = pi^2 * 1.5^2 / 3)
#>     variance pi_sq_over_3 
#>     7.402203     7.402203 

# Symmetric, so the mean is the median; the excess kurtosis is 6/5.
c(mean = mean(d, th), median = distrib_quantile(d, 0.5, th),
  excess_kurtosis = kurtosis(d, th))
#>            mean          median excess_kurtosis 
#>             0.4             0.4             1.2 
```
