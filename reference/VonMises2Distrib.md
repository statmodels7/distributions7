# von Mises Distribution Class, Mean Resultant Length

The S7 class of the von Mises family written in its mean direction
\\\mu\\ and its **mean resultant length** \\\rho \in (0, 1)\\, the
quantity circular statistics reports and one minus the circular
variance. It inherits from `continuous_distrib`; the nine methods listed
below are registered on it directly.

The concentration of
[`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
is recovered as \\\kappa = A^{-1}(\rho)\\ with \\A(\kappa) =
I_1(\kappa)/I_0(\kappa)\\, a strictly increasing bijection from \\(0,
\infty)\\ onto \\(0, 1)\\ whose inverse has no closed form. That is why
the family is written out here instead of through
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

Build one with
[`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md),
which supplies the two link functions and fills the properties in. This
page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
VonMises2Distrib(
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

An S7 object of class `VonMises2Distrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md)
they hold `"von mises2"`, `"univariate"`, `c(-pi, pi)`,
`c("mu", "rho")`, the interpretations
`c(mu = "mean direction", rho = "mean resultant length")`, `2`, and the
domains \\(-\pi, \pi)\\ and \\(0, 1)\\.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.VonMises2Distrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.VonMises2Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.VonMises2Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.VonMises2Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.VonMises2Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.VonMises2Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.VonMises2Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.VonMises2Distrib.md),
[`mean()`](https://statmodels7.github.io/distributions7/reference/mean.VonMises2Distrib.md).

The **quantile** and the response derivatives come from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md),
the quantile by root finding on this class's own distribution function.
The remaining moments come from
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
and its siblings, numerically, and are the ordinary moments of \\Y\\ as
a number.

## See also

[`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md)
to build one;
[`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
for the same law in the concentration;
[`numericals7::bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.html)
for the map;
[`distrib_expected_hessian.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.VonMises2Distrib.md)
for the information.

## Examples

``` r
d <- vonmises2_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The resultant length is bounded, which is what makes it readable.
d@params
#> [1] "mu"  "rho"
d@params_bounds$rho
#> [1] 0 1
d@params_interpretation
#>                      mu                     rho 
#>        "mean direction" "mean resultant length" 

# The direction rides a bounded link and the resultant length a logit.
vapply(d@link_params, function(l) l@link_name, character(1))
#>                                                     mu 
#> "bounded(lwr=-3.14159265358979, upr=3.14159265358979)" 
#>                                                    rho 
#>                                                "logit" 
```
