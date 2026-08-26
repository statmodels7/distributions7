# Elastic-Net Distribution Class

The S7 class of the density whose negative logarithm is the elastic-net
penalty: a Laplace and a Gaussian at the same location, multiplied
together and normalized, \$\$f(y) \propto
\exp\\-\lambda\alpha\|y-\mu\| - \lambda(1-\alpha)(y-\mu)^2/2\\.\$\$ At
\\\alpha \to 1\\ it is the Laplace of
[Laplace2Distrib](https://statmodels7.github.io/distributions7/reference/Laplace2Distrib.md)
and at \\\alpha \to 0\\ the Gaussian; \\\alpha\\ is confined to the open
interval, as every bounded parameter in this package is, so neither end
is a member. Both remain families of their own.

Like the Laplace, its log-likelihood is **not differentiable in
\\\mu\\**: the absolute value carries a kink at the location, and
`params_smooth` records `mu = FALSE`.

Build one with
[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md),
which supplies the three link functions. This page documents the raw S7
constructor, which validates none of the relationships between its
properties.

## Usage

``` r
EnetDistrib(
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

An S7 object of class `EnetDistrib`, inheriting from
`continuous_distrib` and from `distrib`. For an object built by
[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
the properties hold `"enet"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "lambda", "alpha")`, the interpretations
`c(mu = "location", lambda = "rate", alpha = "mixing weight")`, `3`, the
domains \\(-\infty,\infty)\\, \\(0,\infty)\\ and \\(0,1)\\, and
`params_smooth = c(mu = FALSE, lambda = TRUE, alpha = TRUE)`.

## Methods

Registered in this file:
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.EnetDistrib.md),
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.EnetDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.EnetDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.EnetDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.EnetDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.EnetDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.EnetDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.EnetDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.EnetDistrib.md),
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.EnetDistrib.md),
[`mean()`](https://statmodels7.github.io/distributions7/reference/mean.EnetDistrib.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.EnetDistrib.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.EnetDistrib.md).

Registered in `enet_higher.R`, both closed form:
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.EnetDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.EnetDistrib.md).

The **kurtosis** has no closed form here and comes from the base class
by quadrature. So does the fourth response derivative and beyond, which
are zero: the log-density is quadratic in \\y\\ away from the location.

## The two quantities everything is written in

With \\a = \lambda\alpha\\, \\c = \lambda(1-\alpha)\\ and \\x = a/\sqrt
c\\, every derivative in the two rates is a polynomial in \\x\\ and \\G
= \mathrm{d}\log M/\mathrm{d}x\\, where \\M\\ is the Mills ratio. The
map \\(\lambda, \alpha) \mapsto (a, c)\\ is bilinear, so the chain onto
the reported parameters adds one cross term and nothing else.

## See also

[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
to build one;
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
and
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the two ends; `penalties7::elasticnet_penalty()`, the consumer this
family exists for.

## Examples

``` r
d <- enet_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

d@params
#> [1] "mu"     "lambda" "alpha" 
d@params_bounds
#> $mu
#> [1] -Inf  Inf
#> 
#> $lambda
#> [1]   0 Inf
#> 
#> $alpha
#> [1] 0 1
#> 

# The location carries a kink, so it is declared non-smooth. That is what
# switches off the finite-difference guard in check_distrib().
d@params_smooth
#>     mu lambda  alpha 
#>  FALSE   TRUE   TRUE 

# The two ends, approached but not reached.
th <- list(mu = 0, lambda = 2, alpha = 0.5)
c(laplace = distrib_pdf(d, 0.7, list(mu = 0, lambda = 2, alpha = 1 - 1e-10)),
  laplace_exact = distrib_pdf(laplace2_distrib(), 0.7,
                              list(mu = 0, lambda = 2)),
  gaussian = distrib_pdf(d, 0.7, list(mu = 0, lambda = 2, alpha = 1e-10)),
  gaussian_exact = dnorm(0.7, 0, 1 / sqrt(2)))
#>        laplace  laplace_exact       gaussian gaussian_exact 
#>      0.2465970      0.2465970      0.3456374      0.3456374 
```
