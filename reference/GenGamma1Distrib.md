# Generalized Gamma Distribution Class

The S7 class of the generalized gamma in Stacy's three-parameter form: a
scale \\a \> 0\\ and two shapes \\d \> 0\\ and \\p \> 0\\, with density
\\f(y) = p\\y^{d-1}e^{-(y/a)^p}/\\a^d\Gamma(d/p)\\\\ on \\y \> 0\\.

It is the flexible family for a positive response, and it nests four
others exactly: the gamma at \\p = 1\\, the Weibull at \\d = p\\, the
exponential at \\d = p = 1\\ and the half-normal at \\a = \sqrt2, d = 1,
p = 2\\. Choosing between them becomes something to estimate.

Build one with
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md),
which supplies the three link functions. This page documents the raw S7
constructor, which validates none of the relationships between its
properties.

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

Registered in this file, the last three compiled:
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GenGamma1Distrib.md),
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GenGamma1Distrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GenGamma1Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GenGamma1Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GenGamma1Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GenGamma1Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GenGamma1Distrib.md).

Registered elsewhere: the third and fourth orders in
`gengamma1_higher.R`; the response derivatives in
`cross_derivatives_families.R`; the mixed one in
`cross_derivatives_simple.R`; and the four moments in `moments.R`.

## The one representation everything rests on

\\u = (Y/a)^p\\ is Gamma with shape \\k = d/p\\ and unit rate, exactly.
The distribution function is that Gamma's, the quantile function inverts
it, the generator raises a Gamma draw to the power \\1/p\\, and every
expectation the information needs is a moment of \\u\\. It is the same
device the Weibull and the Gumbel use, where the corresponding variable
is standard exponential.

## See also

[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
to build one;
[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md),
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
and
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
for the families it nests;
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
for the limit it approaches.

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

# Three of the four exact special cases, at one point each.
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
