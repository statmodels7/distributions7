# Skew Normal Distribution Class

This is the S7 class of Azzalini's skew normal family, a Gaussian
carrying a shape parameter that tilts it. With \\z = (y-\mu)/\sigma\\
the density is \\f(y) = 2\phi(z)\Phi(\alpha z)/\sigma\\, so the factor
\\2\Phi(\alpha z)\\ is above one on the side where \\\alpha z \> 0\\ and
below one on the other. At \\\alpha = 0\\ that factor is identically one
and the family is the Gaussian.

The class inherits from `continuous_distrib`. Its parametrization is
\\(\mu, \sigma, \alpha)\\ with \\\mu\\ a location and not the mean,
since \\E\[Y\] = \mu + \sigma\delta\sqrt{2/\pi}\\ with \\\delta =
\alpha/\sqrt{1+\alpha^2}\\.

Build one with
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
which supplies the three link functions and fills the properties in.
This page documents the raw S7 constructor, which takes the parent's
properties; the validator of `distrib` checks them, and the constructor
adds no check of its own.

## Usage

``` r
SkewNormal1Distrib(
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

An S7 object of class `SkewNormal1Distrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
they hold `"skew normal1"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "sigma", "alpha")`, the interpretations
`c(mu = "location", sigma = "scale", alpha = "shape")`, `3`, and the
domains \\(-\infty,\infty)\\, \\(0,\infty)\\ and \\(-\infty,\infty)\\.

## Methods

The following methods are registered on this class, and each is linked
to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.SkewNormal1Distrib.md),
  [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.SkewNormal1Distrib.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.SkewNormal1Distrib.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal1Distrib.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewNormal1Distrib.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal1Distrib.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewNormal1Distrib.md)
  and
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.SkewNormal1Distrib.md);

- the derivatives of the expected information:
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal1Distrib.md)
  and
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal1Distrib.md);

- the derivatives in the response:
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormal1Distrib.md),
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.SkewNormal1Distrib.md),
  [`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md)
  and
  [`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md);

- the mixed derivatives in the response and the parameters:
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.SkewNormal1Distrib.md),
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross2_y.md),
  [`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/loc_cross3_y.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_y_hess.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_y_hess.md);

- the derivatives of the distribution function:
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.SkewNormal1Distrib.md),
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.SkewNormal1Distrib.md),
  [`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.SkewNormal1Distrib.md)
  and
  [`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.SkewNormal1Distrib.md);

- the moments:
  [`mean()`](https://statmodels7.github.io/distributions7/reference/mean.SkewNormal1Distrib.md),
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.SkewNormal1Distrib.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewNormal1Distrib.md)
  and
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewNormal1Distrib.md);

- the fitting helpers:
  [`distrib_intercept_start()`](https://statmodels7.github.io/distributions7/reference/distrib_intercept_start.SkewNormal1Distrib.md)
  and
  [`expected_hessian_costly()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_costly.md).

Every other generic is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
to build one;
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
for the same law in Azzalini's centered parametrization, whose
information is not singular at symmetry;
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md),
which adds degrees of freedom and reaches a skewness that this family
cannot reach;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the case \\\alpha = 0\\.

## Examples

``` r
d <- skewnormal1_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

d@params
#> [1] "mu"    "sigma" "alpha"
d@params_interpretation
#>         mu      sigma      alpha 
#> "location"    "scale"    "shape" 
d@params_bounds
#> $mu
#> [1] -Inf  Inf
#> 
#> $sigma
#> [1]   0 Inf
#> 
#> $alpha
#> [1] -Inf  Inf
#> 

# The scale rides a log link; the location and the shape are unconstrained.
vapply(d@link_params, function(l) l@link_name, character(1))
#>         mu      sigma      alpha 
#> "identity"      "log" "identity" 

# mu is a location and not the mean.
th <- list(mu = 0, sigma = 1, alpha = 3)
c(mu = th$mu, mean = mean(d, th))
#>        mu      mean 
#> 0.0000000 0.7569398 
```
