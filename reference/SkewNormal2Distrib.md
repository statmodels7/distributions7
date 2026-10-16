# Skew Normal Distribution Class, Centered Parametrization

This is the S7 class of the skew normal written in its first three
moments: the mean \\\mu\\, the standard deviation \\\sigma\\ and the
skewness \\\gamma_1\\. It is the same law as
[SkewNormal1Distrib](https://statmodels7.github.io/distributions7/reference/SkewNormal1Distrib.md),
reached through the map of
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md);
the two classes differ only in which three numbers name a member of the
family.

The expected information of this parametrization is non-singular at
\\\gamma_1 = 0\\, where the expected information of the direct
parametrization loses a rank; at \\\sigma = 1\\ its eigenvalues tend to
2, 1 and \\1/6\\ as the skewness goes to zero.

Build one with
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md),
which supplies the three link functions and bounds the skewness at \\\pm
0.9952717\\. This page documents the raw S7 constructor, which takes the
parent's properties; the validator of `distrib` checks them, and the
constructor adds no check of its own.

## Usage

``` r
SkewNormal2Distrib(
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

An S7 object of class `SkewNormal2Distrib`, inheriting from
`continuous_distrib` and from `distrib`. For an object built by
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
the properties hold `"skew normal2"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "sigma", "gamma1")`, the interpretations
`c(mu = "mean", sigma = "standard deviation", gamma1 = "skewness")`,
`3`, and the domains \\(-\infty,\infty)\\, \\(0,\infty)\\ and
\\(-0.9952717, 0.9952717)\\.

## Methods

The following methods are registered on this class, and each is linked
to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.SkewNormal2Distrib.md),
  [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.SkewNormal2Distrib.md),
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.SkewNormal2Distrib.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.SkewNormal2Distrib.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal2Distrib.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal2Distrib.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal2Distrib.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal2Distrib.md),
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal2Distrib.md)
  and
  [`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal2Distrib.md);

- the derivatives of the expected information:
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal2Distrib.md)
  and
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal2Distrib.md);

- the derivatives in the response:
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormal2Distrib.md),
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormal2Distrib.md),
  [`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md)
  and
  [`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md);

- the mixed derivatives in the response and the parameters:
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.SkewNormal2Distrib.md),
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormal2Distrib.md),
  [`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/loc_cross3_y.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormal2Distrib.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormal2Distrib.md);

- the derivatives of the distribution function:
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/register_mapped_cdf_k.md),
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/register_mapped_cdf_k.md),
  [`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/register_mapped_cdf_k.md)
  and
  [`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/register_mapped_cdf_k.md);

- the moments:
  [`mean()`](https://statmodels7.github.io/distributions7/reference/mean.SkewNormal2Distrib.md),
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.SkewNormal2Distrib.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewNormal2Distrib.md)
  and
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewNormal2Distrib.md);

- the fitting helpers:
  [`expected_hessian_costly()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_costly.md).

Every other generic is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## The point at zero skewness

The map runs through a cube root, so \\\partial\alpha/\partial\gamma_1\\
is unbounded at \\\gamma_1 = 0\\. The score has a finite limit there and
is returned, and so is the expected information. A derivative of the
log-density with \\k \ge 2\\ differentiations in \\\gamma_1\\ diverges
like \\\gamma_1^{4/3-k}\\ (\\\gamma_1^{-2/3}\\ for \\k = 2\\). The
methods of order two and more, and the derivatives of the expected
information, signal an error at \\\gamma_1 = 0\\ exactly. The density
and the distribution function equal the Gaussian's there.

## See also

[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
to build one;
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
for the direct parametrization;
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md)
for the map between them.

## Examples

``` r
d <- skewnormal2_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

d@params
#> [1] "mu"     "sigma"  "gamma1"
d@params_interpretation
#>                   mu                sigma               gamma1 
#>               "mean" "standard deviation"           "skewness" 

# The skewness is bounded, which the direct parametrization's shape is not.
d@params_bounds$gamma1
#> [1] -0.9952717  0.9952717

# All three parameters are moments, which is why it is called "centered".
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
c(mean = mean(d, th), sd = sqrt(variance(d, th)), skewness = skewness(d, th))
#>     mean       sd skewness 
#>      0.0      1.0      0.5 
```
