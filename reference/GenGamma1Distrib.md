# Generalized Gamma Distribution Class

This S7 class represents the generalized gamma in Stacy's
three-parameter form: a scale \\a \> 0\\ and two shapes \\d \> 0\\ and
\\p \> 0\\, with density \\f(y) =
p\\y^{d-1}e^{-(y/a)^p}/\\a^d\Gamma(d/p)\\\\ on \\y \> 0\\.

It is the flexible family for a positive response, and it nests four
others exactly: the gamma at \\p = 1\\, the Weibull at \\d = p\\, the
exponential at \\d = p = 1\\ and the half-normal at \\a = \sqrt2, d = 1,
p = 2\\. Estimates of \\d\\ and \\p\\ therefore show which of these
special cases the data support.

Build one with
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md),
which supplies the three link functions. This page documents the raw S7
constructor; the validator of `distrib` checks that the properties are
consistent (the parameter count, a bound and a link for every parameter,
ordered support bounds).

## Usage

``` r
GenGamma1Distrib(
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

An S7 object of class `GenGamma1Distrib`, inheriting from
`continuous_distrib` and from `distrib`. For an object built by
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
the properties hold `"gengamma1"`, `"univariate"`, `c(0, Inf)`,
`c("a", "d", "p")`, the interpretations
`c(a = "scale", d = "shape", p = "power")`, `3`, and three domains \\(0,
\infty)\\.

## Methods

The following methods are registered on this class, and each is linked
to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GenGamma1Distrib.md),
  [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GenGamma1Distrib.md),
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GenGamma1Distrib.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GenGamma1Distrib.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GenGamma1Distrib.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GenGamma1Distrib.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GenGamma1Distrib.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GenGamma1Distrib.md),
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GenGamma1Distrib.md)
  and
  [`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GenGamma1Distrib.md);

- the derivatives of the expected information:
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.StudentT1Distrib.md)
  and
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.StudentT1Distrib.md);

- the derivatives in the response:
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GenGamma1Distrib.md),
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GenGamma1Distrib.md),
  [`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GenGamma1Distrib.md)
  and
  [`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GenGamma1Distrib.md);

- the mixed derivatives in the response and the parameters:
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GenGamma1Distrib.md),
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GenGamma1Distrib.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GenGamma1Distrib.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GenGamma1Distrib.md);

- the derivatives of the distribution function:
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md),
  [`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_deriv_cdf_k.md)
  and
  [`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_deriv_cdf_k.md);

- the moments:
  [`mean()`](https://statmodels7.github.io/distributions7/reference/mean.GenGamma1Distrib.md),
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.GenGamma1Distrib.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.GenGamma1Distrib.md)
  and
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.GenGamma1Distrib.md).

Every other generic is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## Representation as a gamma variable

\\u = (Y/a)^p\\ is Gamma with shape \\k = d/p\\ and unit rate, exactly.
The distribution function is that Gamma's, the quantile function inverts
it, the generator raises a Gamma draw to the power \\1/p\\, and every
expectation that the information requires is a moment of \\u\\. The
Weibull and the Gumbel use the same device, with a standard exponential
variable in place of the gamma.

## See also

[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
to build one;
[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md),
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
and
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
for the families that it nests;
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
for a limiting case.

## Examples

``` r
d <- gengamma1_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

d@params
#> [1] "a" "d" "p"
d@params_interpretation
#>       a       d       p 
#> "scale" "shape" "power" 

# All three parameters are positive, so all three ride a log link.
vapply(d@link_params, function(l) l@link_name, character(1))
#>     a     d     p 
#> "log" "log" "log" 

# Three of the four exact special cases are compared at one point each.
y <- c(0.5, 1.5, 4)
rbind(gamma = distrib_pdf(d, y, list(a = 2, d = 3, p = 1)) -
              dgamma(y, shape = 3, scale = 2),
      weibull = distrib_pdf(d, y, list(a = 2, d = 1.5, p = 1.5)) -
                dweibull(y, shape = 1.5, scale = 2),
      half_normal = distrib_pdf(d, y, list(a = sqrt(2), d = 1, p = 2)) -
                    2 * dnorm(y))
#>                      [,1]          [,2]         [,3]
#> gamma        1.734723e-18  1.387779e-17 0.000000e+00
#> weibull     -5.551115e-17  5.551115e-17 4.163336e-17
#> half_normal -1.110223e-16 -5.551115e-17 5.963112e-19
```
