# Skew t Third Derivatives

Computes the twenty third derivatives of the log-density, assembled so
that no stencil is ever applied to another stencil's output. Ten of the
twenty involve \\\nu\\.

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` only its
  length matters.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`.

- expected:

  Logical of length 1. When `FALSE`, the default, the observed
  derivatives at `y` are returned.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- approx:

  One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read only when
  `expected = TRUE`.

- nsim:

  A single positive integer, the Monte Carlo sample size used when
  `approx = "mc"`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of twenty numeric vectors, one per distinct third-order
component, from `mu_mu_mu` to `nu_nu_nu` as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
names them.

## How the twenty are obtained

A component whose Hessian entry is closed form goes through the generic
construction of
[`numerical_deriv3()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv3.md),
which is one stencil on an analytic quantity. That covers both indices
in \\(\mu, \sigma, \alpha)\\, and also one index equal to \\\nu\\ where
the stencil runs along a different variable.

The components the generic construction would nest are replaced: \\(i,
\nu, \nu)\\ for \\i\\ in \\(\mu, \sigma, \alpha)\\ is one five-point
second difference of the **closed-form** score component \\i\\, through
[`fd5_second()`](https://statmodels7.github.io/distributions7/reference/fd5_second.md);
and \\(\nu, \nu, \nu)\\ is one five-point third difference of the
log-density itself, through
[`fd5_third()`](https://statmodels7.github.io/distributions7/reference/fd5_third.md).

## Accuracy

The pure-\\\nu\\ component is the loosest at this order. Measured at
\\\mu = 0\\, \\\sigma = 1\\, \\\alpha = 3\\, \\\nu = 6\\ on four
observations, `nu_nu_nu` is \\-0.0061280\\ against \\-0.0061270\\ from
an independent single stencil on the log-density at \\h = 0.05\\, so
about four significant digits. The package's own thirteen-check battery
reports order 3 against finite differences at \\2\times10^{-6}\\ for
this family.

With `expected = TRUE` the whole order is an expectation and comes from
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md);
the family has no closed-form expected information, so `approx` and
`nsim` are read.

## See also

[`distrib_hessian.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewTDistrib.md)
for the order below,
[`distrib_deriv4.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.SkewTDistrib.md)
for the order above,
[`fd5_third()`](https://statmodels7.github.io/distributions7/reference/fd5_third.md)
for the stencil, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- skewt_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
d3 <- distrib_deriv3(d, y, th)
c(components = length(d3), involving_nu = sum(grepl("nu", names(d3))))
#>   components involving_nu 
#>           20           10 

# A closed-form-block component against a difference of the Hessian.
eps <- 1e-5
rbind(analytic = d3$mu_mu_alpha,
      numeric = (distrib_hessian(d, y, list(mu = 0, sigma = 1,
                                            alpha = 3 + eps, nu = 6))$mu_mu -
                 distrib_hessian(d, y, list(mu = 0, sigma = 1,
                                            alpha = 3 - eps, nu = 6))$mu_mu) /
                (2 * eps))
#>               [,1]      [,2]       [,3]        [,4]
#> analytic 0.6111686 -1.082016 -0.4214539 0.004008904
#> numeric  0.6111686 -1.082016 -0.4214539 0.004008904

# The pure-nu component against an independent single stencil on the
# log-density, which shares no arithmetic with the route above.
ld <- function(v) sum(distrib_pdf(d, y, list(mu = 0, sigma = 1,
                                             alpha = 3, nu = v), log = TRUE))
c(ours = sum(d3$nu_nu_nu),
  stencil = numericals7::fd_derivative(ld, 6, 3L, h = 0.05))
#>         ours      stencil 
#> -0.006128010 -0.006127044 
```
