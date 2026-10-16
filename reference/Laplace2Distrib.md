# Laplace Distribution Class, Location and Rate

This is the S7 class of the Laplace (double exponential) family written
by its rate \\\lambda \> 0\\ in place of its scale, with density \\f(y)
= (\lambda/2)\exp(-\lambda\|y-\mu\|)\\ on the whole real line. It is the
same law as
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
at \\\lambda = 1/\sigma\\; the parametrization differs and the
derivatives with it. It inherits from `continuous_distrib`, so every
generic of the `distrib` class applies to it.

The rate form is the one in which a lasso penalty is written: with \\\mu
= 0\\ the negative log-density is \\\lambda\|y\| - \log(\lambda/2)\\, so
\\\lambda\\ is the penalty's own tuning parameter and larger values
shrink harder. `penalties7::lasso_penalty()` is this family when the
location is held at zero.

Like
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md),
this family has a kink at \\y = \mu\\ and is non-regular in its
location; `params_smooth` is `c(mu = FALSE, lambda = TRUE)`.

Build one with
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).
This page documents the raw S7 constructor, which takes the parent's
properties; the validator of `distrib` checks them, and the constructor
adds no check of its own.

## Usage

``` r
Laplace2Distrib(
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

An S7 object of class `Laplace2Distrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
they hold `"laplace2"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "lambda")`, the interpretations
`c(mu = "location", lambda = "rate")`, `2`, the domains \\(-\infty,
\infty)\\ and \\(0, \infty)\\, the two links, and
`c(mu = FALSE, lambda = TRUE)` for `params_smooth`.

## Methods

The following methods are registered on this class, and each is linked
to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Laplace2Distrib.md),
  [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Laplace2Distrib.md),
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Laplace2Distrib.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Laplace2Distrib.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Laplace2Distrib.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Laplace2Distrib.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Laplace2Distrib.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Laplace2Distrib.md),
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Laplace2Distrib.md)
  and
  [`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.Laplace2Distrib.md);

- the derivatives of the expected information:
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.elementary.md)
  and
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.elementary.md);

- the derivatives in the response:
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Laplace2Distrib.md),
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Laplace2Distrib.md),
  [`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md)
  and
  [`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md);

- the mixed derivatives in the response and the parameters:
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Laplace2Distrib.md),
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/mapped_theta2_methods.md),
  [`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/loc_cross3_y.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/mapped_theta2_methods.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/mapped_theta2_methods.md);

- the derivatives of the distribution function:
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Laplace2Distrib.md),
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.Laplace2Distrib.md),
  [`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/register_mapped_cdf_k.md)
  and
  [`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/register_mapped_cdf_k.md);

- the moments:
  [`mean()`](https://statmodels7.github.io/distributions7/reference/mean.Laplace2Distrib.md),
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.Laplace2Distrib.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.Laplace2Distrib.md)
  and
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Laplace2Distrib.md);

- the decomposition at the kink:
  [`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.Laplace2Distrib.md).

Every other generic is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
to build one;
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
for the same law in its scale;
[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md),
which mixes this family with a Gaussian.

## Examples

``` r
d <- laplace2_distrib()
d@params
#> [1] "mu"     "lambda"
d@params_interpretation
#>         mu     lambda 
#> "location"     "rate" 
d@params_smooth
#>     mu lambda 
#>  FALSE   TRUE 

# This is the same law as laplace_distrib() at lambda = 1/sigma.
y <- c(-1.2, 0.3, 2.5)
all.equal(distrib_pdf(d, y, list(mu = 0.4, lambda = 1 / 1.5)),
          distrib_pdf(laplace_distrib(), y, list(mu = 0.4, sigma = 1.5)))
#> [1] TRUE

# The variance is 2/lambda^2, so a larger rate is a tighter distribution.
vapply(c(0.5, 1, 2), function(l) variance(d, list(mu = 0, lambda = l)),
       numeric(1))
#> [1] 8.0 2.0 0.5
```
