# Skew t Third Derivatives

Computes the twenty third derivatives of the log-density with the
compiled kernel described on
[`distrib_gradient.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md).

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

- threads:

  A single positive integer, how many threads the kernel may use.

## Value

A named list of twenty numeric vectors, one per distinct third-order
component, from `mu_mu_mu` to `nu_nu_nu` as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
names them.

## Details

With `expected = TRUE` the whole order is an expectation and comes from
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md);
the family has no closed-form expected information, so `approx` and
`nsim` are read.

## See also

[`distrib_hessian.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewTDistrib.md)
for the order below,
[`distrib_deriv4.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.SkewTDistrib.md)
for the order above, and
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

# The pure-nu component against a single stencil on the log-density.
ld <- function(v) sum(distrib_pdf(d, y, list(mu = 0, sigma = 1,
                                             alpha = 3, nu = v), log = TRUE))
c(ours = sum(d3$nu_nu_nu),
  stencil = numericals7::fd_derivative(ld, 6, 3L, h = 0.05))
#>         ours      stencil 
#> -0.006128039 -0.006127044 
```
