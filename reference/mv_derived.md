# Interpretable Quantities of a Multivariate Distribution

Returns the quantities a reader of a multivariate fit wants: standard
deviations, correlations, and whatever else the family's matrix
parametrization is about. It returns with them the Jacobian that carries
standard errors onto those quantities, and the scale each one's
confidence interval should be built on.

## Usage

``` r
mv_derived(distrib, theta, ...)
```

## Arguments

- distrib:

  An object inheriting from
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

- theta:

  A named list or vector of parameters, each component a single number.
  Aligned by the generic before dispatch, so any order is accepted.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A named list with

- `value`:

  a named numeric vector of the quantities;

- `jacobian`:

  a numeric matrix with one row per quantity and one column per
  parameter of `distrib`, in `distrib@params` order;

- `transform`:

  a character vector, one of `"identity"`, `"log"`, `"atanh"` or
  `"logit"` per quantity, naming the scale its interval is built on;

- `block`:

  a character vector labeling the group each quantity belongs to, by
  which the printed summary is laid out.

All four are named after the quantities and are the same length.

## Why coordinates are not quantities

The free values of a parameters7 parametrization are coordinates chosen
so that an optimizer can move freely. The logarithm of the second
diagonal entry of a Cholesky factor has an estimate and a standard
error, and neither answers a question anybody asked. This generic names
the quantities that do and supplies \\\partial g/\partial\theta\\, so
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
can apply the delta method to them.

## The interval scale each quantity declares

Each quantity carries the scale its interval should be built on, exactly
as
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
builds a univariate interval on the link scale and maps it back. A
standard deviation is intervalled on the log scale, so the interval
cannot reach zero; a correlation on Fisher's \\z =
\mathrm{artanh}(\rho)\\, so it cannot leave \\(-1, 1)\\; an
unconstrained quantity on its own scale. On the raw scale a
correlation's interval routinely exceeds one.

## What the default method reports

The method registered on
[`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
returns the distinct entries of the matrix
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
produces, named `sigma_v1_v1` and so on, with a Jacobian from one
central difference per parameter. A family whose matrix is not a
covariance still reports it on its own scale, which beats reporting a
Cholesky coordinate. The two shipped families override with closed-form
Jacobians.

## Notation

\\\theta\\ is the parameter vector, \\g\\ the map to the reported
quantities, \\\Sigma\\ the matrix the family carries and \\\rho\\ a
correlation.

## See also

[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md),
which applies the delta method to this,
[`mv_derived.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_derived.MvGaussianDistrib.md)
and
[`mv_derived.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_derived.MvStudentTDistrib.md)
for the two closed-form methods, and
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
for the matrix.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
der <- mv_derived(d, theta)
der$value
#>     sd_v1     sd_v2 cor_v1_v2 
#> 1.0000000 1.1180340 0.4472136 

# Each quantity declares the scale its interval is built on.
der$transform
#>     sd_v1     sd_v2 cor_v1_v2 
#>     "log"     "log"   "atanh" 
der$block
#>                 sd_v1                 sd_v2             cor_v1_v2 
#> "Standard deviations" "Standard deviations"        "Correlations" 

# The Jacobian is exact for this family, against a numerical one.
g <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  mv_derived(d, t2)$value
}
max(abs(der$jacobian - numDeriv::jacobian(g, unlist(theta))))
#> [1] 8.833823e-12

# The mean parameters do not enter any of them, so those columns are zero.
der$jacobian[, c("mu1", "mu2")]
#>           mu1 mu2
#> sd_v1       0   0
#> sd_v2       0   0
#> cor_v1_v2   0   0

# A structured matrix reports its own quantities as a further block.
a <- mvgaussian_distrib(3, sigma = parameters7::ar1(3))
mv_derived(a, as.list(stats::setNames(c(0, 0, 0, 0.1, 0.3), a@params)))$value
#>      sd_v1      sd_v2      sd_v3  cor_v1_v2  cor_v1_v3  cor_v2_v3      scale 
#> 1.05127110 1.05127110 1.05127110 0.29131261 0.08486304 0.29131261 1.10517092 
#>        rho 
#> 0.29131261 
```
