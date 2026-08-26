# Skew t Distribution Class

The S7 class of Azzalini's skew \\t\\: a Student \\t\\ carrying a shape
parameter that tilts it, so that the tail weight and the asymmetry are
modeled by two parameters instead of one. With \\z = (y-\mu)/\sigma\\
and \\w = \alpha z\sqrt{(\nu+1)/(\nu+z^2)}\\ the density is \\2
t\_\nu(z)T\_{\nu+1}(w)/\sigma\\.

It contains three families as limits. At \\\alpha = 0\\ it is the
Student \\t\\; as \\\nu \to \infty\\ it is the skew normal, approached
at \\O(1/\nu)\\; and with both it is the Gaussian. Its reason to exist
is the skewness: the skew normal cannot pass 0.9953, and this family
reaches 2.05 at \\\nu = 6\\ and 4.00 at \\\nu = 4\\.

Build one with
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md),
which supplies the four link functions. This page documents the raw S7
constructor, which validates none of the relationships between its
properties.

## Usage

``` r
SkewTDistrib(
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

An S7 object of class `SkewTDistrib`, inheriting from
`continuous_distrib` and from `distrib`. For an object built by
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
the properties hold `"skew t"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "sigma", "alpha", "nu")`, the interpretations
`c(mu = "location", sigma = "scale", alpha = "shape", nu = "degrees of freedom")`,
`4`, and the domains \\(-\infty,\infty)\\, \\(0,\infty)\\,
\\(-\infty,\infty)\\ and \\(0,\infty)\\.

## Methods

Registered in this file:
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.SkewTDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.SkewTDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewTDistrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewTDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.SkewTDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewTDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.SkewTDistrib.md).

Registered elsewhere: all four moments in `moments.R`
([`mean()`](https://statmodels7.github.io/distributions7/reference/mean.SkewTDistrib.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.SkewTDistrib.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewTDistrib.md),
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewTDistrib.md));
the mixed response-parameter derivative
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
in `cross_derivatives_families.R`; and
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
in `cdf_derivatives_families.R`.

The **distribution function** and the **quantile function** come from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md),
by quadrature and by root finding on it. So does the **expected
information**: this family has none in elementary form, so
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
approximates it, and `method = "newton"` is much the cheaper way to fit
it.

## What is closed form and what is not

Every derivative in \\\mu\\, \\\sigma\\ and \\\alpha\\ is closed form.
Every derivative involving \\\nu\\ is not, and the obstruction is
mathematical: the density carries \\T\_{\nu+1}\\, and the derivative of
a Student \\t\\ distribution function in its degrees of freedom has no
elementary expression. Those components come from **one** stencil
applied to an analytic quantity, never from a difference of a
difference.

## See also

[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
to build one;
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
for the \\\nu \to \infty\\ limit;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
for the \\\alpha = 0\\ case;
[`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md)
for the scalar functions the derivatives are built from.

## Examples

``` r
d <- skewt_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

d@params
#> [1] "mu"    "sigma" "alpha" "nu"   
d@params_interpretation
#>                   mu                sigma                alpha 
#>           "location"              "scale"              "shape" 
#>                   nu 
#> "degrees of freedom" 

# Two parameters for two departures from the Gaussian, which is why the
# family exists.
vapply(d@link_params, function(l) l@link_name, character(1))
#>         mu      sigma      alpha         nu 
#> "identity"      "log" "identity"      "log" 

# It passes the skew normal's ceiling of 0.9953.
vapply(c(4, 6, 20, 1e6),
       function(v) skewness(d, list(mu = 0, sigma = 1, alpha = 50, nu = v)), 0)
#> [1] 3.9976019 2.0500317 1.1904904 0.9936341
```
